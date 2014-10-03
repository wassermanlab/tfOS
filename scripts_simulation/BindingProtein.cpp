//
// BindingProtein.cpp
//
// Author:
//	David Arenillas (dave@cmmt.ubc.ca)
//	Rebecca Worsley Hunt 
//
// Copywrite:
//	Wasserman Lab
//	Centre For Molecular Medicine and Therapeutics
//	Child & Family Research Institute
//	University of British Columbia
//
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <utility>
#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/BindingProtein.h"

/* fallOff: uniform = 1; gaussian = 2; XX = 3; etc. Thus all proteins have same fallOff distr.
      but could have proteins with different slide lengths one day... which would change size of distr */
const unsigned int  BindingProtein::BP_MAX_HOP = 45; // limiting distance to hop from where a protein leaves the DNA
const double  		BindingProtein::BP_HOP_PROB = 0.70; // 7/10 times the protein will hop, thus 3/10 times it will jump 
const unsigned int	BindingProtein::BP_MAX_SLIDE = 42;
// use kJ not J to roughly relate matrix scores to the same scale as energy 
const bp_energy_t	BindingProtein::BP_R = 8.314472/1000; // kJ/mol*K
//AVOGADRO =  6.02214179*10^23 mol^-1
const bp_energy_t	BindingProtein::BP_k = BP_R/(6.02214179*pow(10.0, 23.0)); // kJ/molecule*K
const bp_energy_t	BindingProtein::BP_T = 300;  // Kelvin

using namespace std;

// default constructor - initialize to empty values
BindingProtein::BindingProtein()
{
    myName = "";
	myUniqueID = ""; 		// first protein has id=0 e.g. Myf-0, Myf-1, Myf-2 
	myPioneerStatus = false;
	myLength = 0;
    myMaxSlide = BP_MAX_SLIDE;  //XXX perhaps should be in user input file
	myMaxHop = BP_MAX_HOP;
	myHopProb = BP_HOP_PROB;
	myMaxEnergy = 0;
	myThresholdBindingEnergy = 0;
	myState = BP_FREE;
	myPosition = -1; // position is on DNA molecule; position1 is for 1st nt 5' +ve and 3' -ve strand
    myFreeOrigin = -1; // DNA location which the now free protein originated from
    myDirection = 0;   // -1 is left, +1 is right. Irrespective of orientation.
    mySlidCount = 0;
	myPalindromic = false; // when palindromic=T orientation of protein doesnt matter
	myInteractingProt = NULL;
}

// constructor - initialize to given values
// should probably take out 'optional' things that user does not have to set
BindingProtein::BindingProtein(const string &name, const string &uniqueID, unsigned int length, unsigned int distance, unsigned int maxHop,
	double hopProb, double energy, double thresholdBinding, const vector<bp_energy_t> &fep, const vector<bp_energy_t> &rep ):
	 myName(name), myUniqueID(uniqueID), myLength(length), myMaxSlide(distance), myMaxHop(maxHop), myHopProb(hopProb), myMaxEnergy(energy), myThresholdBindingEnergy(thresholdBinding), 
	  myFEP(fep), myREP(rep)
{
	myPioneerStatus = false;
    myState = BP_FREE;
    myPosition = -1;
    myFreeOrigin = -1;	// DNA location the now free protein originated from
    myDirection = 0;	// -1 == left, +1 == right. Irrespective of orientation.
    mySlidCount = 0;
    myPalindromic = false;
	myInteractingProt = NULL;
}

BindingProtein::~BindingProtein()
{
}

//copy constructor. Only used when making multiple copies of a protein at beginning of simulation
BindingProtein::BindingProtein(const BindingProtein& orig_prot)
{
	myName = orig_prot.myName;
	myUniqueID = "";
	myPioneerStatus = orig_prot.myPioneerStatus;
	myLength = orig_prot.myLength;
	myState = orig_prot.myState;
	myPosition = orig_prot.myPosition;
	myFreeOrigin = orig_prot.myFreeOrigin;
	myDirection = orig_prot.myDirection;
	mySlidCount = orig_prot.mySlidCount;
	myMaxSlide = orig_prot.myMaxSlide;
	myMaxHop = orig_prot.myMaxHop;
	myHopProb = orig_prot.myHopProb;
	myMaxEnergy = orig_prot.myMaxEnergy;
	myThresholdBindingEnergy = orig_prot.myThresholdBindingEnergy;
	myPalindromic = orig_prot.myPalindromic;
	myFEP = orig_prot.myFEP;
	myREP = orig_prot.myREP;
	myInteractingProt = orig_prot.myInteractingProt;
}

// = operator overload (rhs is right-hand-side of the operator)
BindingProtein& BindingProtein::operator= (const BindingProtein& rhs_prot)
{
	if (this == &rhs_prot) return *this;
	
	myName = rhs_prot.myName;
	myUniqueID = "=";
	myPioneerStatus = rhs_prot.myPioneerStatus;
	myLength = rhs_prot.myLength;
	myState = rhs_prot.myState;
	myPosition = rhs_prot.myPosition;
	myFreeOrigin = rhs_prot.myFreeOrigin;
	myDirection = rhs_prot.myDirection;
	mySlidCount = rhs_prot.mySlidCount;
	myMaxSlide = rhs_prot.myMaxSlide;
	myMaxHop = rhs_prot.myMaxHop;
	myHopProb = rhs_prot.myHopProb;
	myMaxEnergy = rhs_prot.myMaxEnergy;
	myThresholdBindingEnergy = rhs_prot.myThresholdBindingEnergy;
	myPalindromic = rhs_prot.myPalindromic;
	myFEP = rhs_prot.myFEP; 				// Forward Energy Positions
	myREP = rhs_prot.myREP;
	myInteractingProt = rhs_prot.myInteractingProt;

	return *this;
}


// Protein name setter
void BindingProtein::name(const string &name)
{
    myName = name;
}

// Protein ID setter
void BindingProtein::uniqueID(const string &uniqueID)
{
	myUniqueID = uniqueID;
}

void BindingProtein::isPioneer( bool pioneer)
{
	myPioneerStatus = pioneer;
}

// Protein length setter (width of TFBS profile matrix) setter
void BindingProtein::length(unsigned int length)
{
    myLength = length;
}

// Protein position on DNA
void BindingProtein::position(int pos)
{
    myPosition = pos;
}

// Protein position prior to falling off the DNA
void BindingProtein::freeOrigin(int pos)
{
    myFreeOrigin = pos;
}

// Direction protein is moving on DNA. 0 = none; 1 = right; -1 = left.
void BindingProtein::direction(int dir)
{
    myDirection = dir;
}

// Maximum number of nucleotides before falling off
void BindingProtein::maxSlide(unsigned int distance) 
{
	myMaxSlide = distance;
}

// maximum dist protein is allowed to hop
void BindingProtein::maxHopDistance(unsigned int maxHopDist )
{
	myMaxHop = maxHopDist;
}

// the chance of hopping (1-hopProb = chance of jumping)
void BindingProtein::hopnotjumpProb(float hopProb)
{
	myHopProb = hopProb;
}

// Is this protein palindromic?
void BindingProtein::isPalindromic(bool palindromic)
{
    myPalindromic = palindromic;
}

/* Which strand is this protein oriented to?  i.e. looking down imaginary axis and assigning left and 
	right side to protein. 1 = left-side oriented to +ve strand, -1 = left-side oriented to -ve strand.
	Won't matter for palindromic */
void BindingProtein::strandOrientation (int strandOrientation)
{
    myStrandOrientation = strandOrientation;
}

// State of protein (free, tethered or bound
void BindingProtein::state(bp_state_t state)
{
    myState = state;
}

// Maximun binding energy of protein (maximum possible score of energy matrix)
void BindingProtein::maxEnergy(bp_energy_t energy)
{
    myMaxEnergy = energy;
}

// Treshold binding energy
void BindingProtein::thresholdBindingEnergy(bp_energy_t thresholdBinding)
{
	myThresholdBindingEnergy = thresholdBinding;
}

// Forward binding energy profile
void BindingProtein::forwardEnergyProfile(const vector<bp_energy_t> &fep)
{
    myFEP = fep;
}

// Reverse binding energy profile
void BindingProtein::reverseEnergyProfile(const vector<bp_energy_t> &rep)
{
    myREP = rep;
}

void BindingProtein::addForwardEnergy(bp_energy_t val)
{
    myFEP.push_back(val);
}

void BindingProtein::addReverseEnergy(bp_energy_t val)
{
    myREP.push_back(val);
}

// Functions so that protein knows where it is located on the sequence, and its state etc
// Set protein tethered pos; optionally set direction (if no strandOrientation given then default=0).
void BindingProtein::tetherAt(int pos, int dir, int strandOrientation = 0)
{
    myState = BP_TETHERED;
    myPosition = pos;
    myDirection = dir;
    myStrandOrientation = strandOrientation;
    mySlidCount = 0;
}

// Set protein bound pos
void BindingProtein::bindAt(int pos)
{
    myState = BP_BOUND;
    myPosition = pos;
    myDirection = 0;
    mySlidCount = 0;
}

// Set protein bound pos
void BindingProtein::bindAt(int pos, int strandOrientation = 0) //if no orientation given it'll be set to 0
{
    myState = BP_BOUND;
    myPosition = pos;
    myDirection = 0;
    myStrandOrientation = strandOrientation;
    mySlidCount = 0;
}

// Protein falls off DNA
void BindingProtein::fallOff()
{
    myState = BP_FREE;
    myFreeOrigin = myPosition; // remember relative location to DNA
    myPosition = -1;
    myDirection = 0;
    myStrandOrientation = 0;
    mySlidCount = 0;
	if(myInteractingProt != NULL){
		myInteractingProt->setInteractingProt(NULL);
	}
	myInteractingProt = NULL;
}

// Protein slides right along DNA
void BindingProtein::slideRight()
{
    myDirection = 1;
    myPosition++;	// XXX bounds check?
    mySlidCount++;
}

// Protein slides left along DNA
void BindingProtein::slideLeft()
{
    myDirection = -1;
    myPosition--;	// XXX bounds check?
    mySlidCount++;
}

// Protein slides in current direction
void BindingProtein::slideCurrent()
{
    myPosition += myDirection;	// XXX bounds check?
    mySlidCount++;
}

// Return probability of protein transitioning from a free state to a tethered
// state
double BindingProtein::freeToTetheredProb(int pos) const
{
    if (myState != BP_FREE) return 0;

    return 1;
}

// Return probability of protein transitioning from a free state to a bound
// state.  Currently no chance
double BindingProtein::freeToBoundProb(int pos) const
{
	if (myState != BP_FREE) return 0;
    
    //return bindingProbability();
    return 0;
}

// Protein is free, will it Hop or will it Jump
double BindingProtein::hopnotjumpProb() const
{
    if (myState != BP_FREE) return 0;
    
    return 0.70; // this is returned to testProb and thus decides likelihood
				// of rand picking a value < 0.70
	// around 30 of 100 decisions is a jump, rest of it is a hop
	// 50% of 0.3 = 0.15.  1-0.15 converts it to hop i.e. rand<(1-0.15) = hop
}					

// Return probability of protein transitioning from a tethered state to a
// bound state
double BindingProtein::tetheredToBoundProb(double seqAccess, pair<pair<float,string>, pair<double,double> > weightEnergy ) 
{
    if (myState != BP_TETHERED) return 0;
	
	return bindingProbability(seqAccess, weightEnergy );
}


// Return probability of protein transitioning from a tethered state to a free state
// uniform = 1, gaussian = 2 
double BindingProtein::tetheredToFreeProb(unsigned int fallOffDistr) const
{
	// fallOff_distr is from ProteinBindingSimulator right now, since other simulator defaults  
	if (myState != BP_TETHERED) return 0;
    
	double prob=0; /* the "test tetheredToFreeProb" comes before a protein slides, 
						thus position 0 has prob=0 of disengagement */  
	double pi = 3.14159265358979323846264338;

	if(mySlidCount != 0) { // position 0, is where protein sits before sliding
		double gaus_stdev = 0.25*myMaxSlide; // 0.15 works to keep distr bounded
		double gaus_mean = 1*myMaxSlide; // again 0.6 keeps the distr bounded. Best for when using prob*log(position*10)
							// may not be needed when not using log(position*10).
							// However, not using log drastically reduces number of proteins falling off. If mean is 30,
							// then expect 2000 proteins fall off over ~60,000 slides. Get that with the log, but 10fold 
							//reduced without it.
		switch (fallOffDistr) {
			case 1:	 // uniform
				// 2014.08.21 I thought it wasnt working... but 'mySlideCount' on first tethering got up to 52bp for 1 of 4 proteins and 32bp for others
				prob = 1/((double)myMaxSlide - (double)(mySlidCount-1)); 
				//cout << "uniform prob: "<< myMaxSlide<<"-"<<mySlidCount<<" "<< prob<< " slideposition\t"<< mySlidCount <<endl; //TTT
				break;
			case 2:  // 2014.08.21 below I said equation doesnt work... however the 'cout' comment shows it is working
			         //  I wish the results were skewed... in R they are, but not how it works here 
			         /* gaussian distribution = 
						   1/(stdev*sqrt(2*pi)) * e^(-(x-mean)^2/(2*stdev^2)), per Wolfram where x is the x-axis value */
						   //  1/(stdev*sqrt(2*pi)) * e^(-0.5*((x-mean)/stdev)^2), per Wikipedia where x is the x-axis value 
				prob = ( 1/(gaus_stdev*sqrt(2*pi)) )*exp( -pow((double)mySlidCount-gaus_mean,2)/(2*pow(gaus_stdev,2)) );	
				prob *= log(double(mySlidCount*10)); // this adjustment increases the height of the gaussian curve of probabilities
				//cout << "gaussian prob: "<< prob<< " slideposition\t"<< mySlidCount << endl; //TTT
				break;
			default: // gaussian
				prob = ( 1/(gaus_stdev*sqrt(2*pi)) )*exp( -pow((double)mySlidCount-(gaus_mean),2)/(2*pow(gaus_stdev,2)) );
				prob *= log(double(mySlidCount*10)); // this adjustment keeps right tail under control
		}
	}	
	return prob;
}				  


// Return probability of protein remaining in bound state
double BindingProtein::remainBoundProb( double seqAccess, pair<pair<float,string>, pair<double, double> > weightEnergy ) 
{
    if (myState != BP_BOUND) return 0;

    return bindingProbability(seqAccess, weightEnergy);
}

// Return probability of protein transitioning from a bound state to a
// tethered state. We don't allow bound to Free, so bound to tethered p=1
double BindingProtein::boundToTetheredProb() const
{
    if (myState != BP_BOUND) return 0;

    return 0; //XXX really, this should be allowed if binding site is weak.
			  //  as weak disengagement is likely similar to tether sliding?
}

// Return probability of protein transitioning from a bound state to a
// free state
double BindingProtein::boundToFreeProb() const
{
    if (myState != BP_BOUND) return 0;
	
	return 1; //XXX always release, assuming some kinetic energy involved. Should allow transition to tether
	          //  particularly for weaker sites
    //return 1 - bindingProbability(); // this is actually the prob. of binding no sites along entire DNA
}

double BindingProtein::energyAt(int pos) const
{
	bp_energy_t fEnrg = myFEP[pos]; 
	bp_energy_t rEnrg = myREP[pos];
	if (myStrandOrientation == 1)
		return fEnrg;
	else if (myStrandOrientation == -1)
		return rEnrg;
	else{
		return max(fEnrg, rEnrg);
	}
}

void BindingProtein::setInteractingProt( BindingProtein* interactingProtein)
{
	myInteractingProt = interactingProtein;
}

// Compute probability of binding at the current position
/* not using Boltzmann or Granek & Clarke.  Simulation is taking care of probability of landing on a
   segment of DNA.  What we want to know is how likely the TF is to bind there for a while... which
   is based on the affinity for the site. */
double BindingProtein::bindingProbability() const 
{
    bp_energy_t fEner = myFEP[myPosition];
    bp_energy_t rEner = myREP[myPosition];
    bp_energy_t dG = 0; // dG = sum(RTlnKa) sum positions in a site, then to get site Ka, use e^(dG/RT)

	/* XXX do something: the maximum energy needs to more than the maximum from the DNA seq... otherwise the
	   probability gets stuck at p=1 */
	
    if (myPalindromic) {  /* because palindromic, orientation doesn't matter... take orientation with
    					     highest energy */
    	dG = max(fEner, rEner);
    }else {
    	if (myStrandOrientation == 1)
	    	dG = fEner;
    	else if (myStrandOrientation == -1)
	    	dG = rEner;
    }
	//double kT = BP_k * BP_T;  // kJ molecule-1
    double RT = BP_R * BP_T; // kJ mol-1

    //XXX setting thresholdBindingEnergy as roughly the 50% probability of binding
    double p = 1 / (1 + exp( (double) -1*(dG-thresholdBindingEnergy() )/RT )); //sigmoidal of dG for sequence versus binding threshold
    
    // if ( p > 0.5) 
	//	cout << "  "<< myName <<"\tBind prob. = " << p << "\t  dG= " << dG << "\tthreshold= "<<thresholdBindingEnergy() << endl;

	return p;
}

// not const because might change this->myInteractingProt
double BindingProtein::bindingProbability(double seqaccess, pair<pair<float,string>, pair<double,double> > weightNrg) 
{
    bp_energy_t fEner = myFEP[myPosition];
    bp_energy_t rEner = myREP[myPosition];
    bp_energy_t dG = 0; // dG = sum(RTlnKa) sum positions in a site, then to get site Ka, use e^(dG/RT)
	double prob = 0;
	float weight = (weightNrg.first).first; // weight of interaction
	// coopbindingstate value is "onebound" if only one protein needs to be at a binding site; "bothbound" or both proteins must be at a binding site
	//         "onebound" is most useful for when had a matrix that represents two binding proteins, which can mean one or both proteins able to bind there
	//         although ultimately one wonders if any site for two proteins "requires" two proteins and won't function with only one
	string coopbindingstate = (weightNrg.first).second; // interaction requirement
	double prot2Nrg = (weightNrg.second).first;   // 2nd protein's energy for dna
	double prot2NrgThr = (weightNrg.second).second; // 2nd protein's bindng energy threshold
	int accessThr = 0; // sequence isn't accessible to binding
	int weightThr = 0; // no interaction weight
	
	//XXX seqaccess: one day may want to use the value to decide how likely it is the TF can access seq, rather than the binary system used at the moment
	//if sequence not accessible: return, unless protein is a pioneer
	//if( seqaccess < accessThr ){ 
	if( seqaccess < accessThr && myPioneerStatus == false){ 
		return prob=0; // not accessible
	}

    if (myPalindromic) {  /* because palindromic, orientation doesn't matter... take orientation with
    					     highest energy */
    	dG = max(fEner, rEner);
    }else {
    	if (myStrandOrientation == 1)
	    	dG = fEner;
    	else if (myStrandOrientation == -1)
	    	dG = rEner;
    }
	//===== prob equation ===============================

	//double kT = BP_k * BP_T;  // kJ molecule-1
    double RT = BP_R * BP_T; // kJ mol-1
    //XXX setting thresholdBindingEnergy as roughly the 50% probability of binding
    // (dG-thresholdBindingEnergy())/RT will be +ve for dG above threshold,  which then becomes a -ve exponential: 1+exp^-X,  thus denominator gets tiny with strong dG and prob appraoches 1
	prob = 1 / (1 + exp( (double) -1*(dG-thresholdBindingEnergy())/RT )); //sigmoidal of dG for sequence versus binding threshold
	
	//====================================================
	if ( prob > 0.5) //TTT
		cout << "BP line 471: "<< myUniqueID <<"\tp= " << prob << "  dG= " << dG << " thr= "<<thresholdBindingEnergy() << endl; //TTT
	/* Interactions requirements:
	       1) both proteins prob. needs to be >0.50 for interaction to happen (because able to bind specifically regardless of result from testProb())
	       2) the strongest DNA affinity of prot1 or prot2 becomes the base binding energy for the two
	       3) then the interaction weight increases that energy */ 
	if( weight > weightThr ){
		double prot2prob = 1 / (1 + exp( (double) -1*( prot2Nrg - prot2NrgThr )/RT ));
		double bindingthreshold = 0.50; //when prob >0.5 then dG >thresholdBindingEnergy(),
		
		cout <<"BP line 481: energies "<< dG << " (thr "<< thresholdBindingEnergy()<<") prot2 "<<prot2Nrg <<" (thr "<<prot2NrgThr<<")" <<endl;//TTT
		if( prot2prob > prob){ // keep the stronger affinity for co-binding DNA affinity  
			prob = prot2prob;	
		}
		cout << "BP: line 485:  p="<<prob <<" ";//TTT
		if(coopbindingstate == "onebound" && (prob >= bindingthreshold || prot2prob >= bindingthreshold) ){ //only one protein needs to be above bindingthreshold
			prob = prob*(double) weight;	
		 	cout << "BP: line 488: one bound, weight added; p becomes p="<<prob << "(" << myUniqueID<<" +partner)" <<endl;//TTT
		}else if(coopbindingstate == "bothbound" && prob >= bindingthreshold && prot2prob >= bindingthreshold ){ // require both proteins to be at least 50% able to bind 
			cout << "BP: tested both binding > 0.50: "<<dG << " ("<< prob<< ") vs "<< prot2Nrg << " ("<< prot2prob<< ")"<<endl;
			//cout << "BP: both pass 0.50: "<<dG << " ("<< prob<< ") vs "<< prot2Nrg << " ("<< prot2prob<< ")" <<endl; //TTT
			prob = prob*(double) weight;	
		 	cout << "BP: line 493: weight added to binding "<< weight << " becomes p="<<prob << "(" << myUniqueID<<" +partner)"<<endl;//TTT
		}else{
			// ProteinBindingSimulator set putative interactions, but they just failed so erase them
			if(myInteractingProt !=NULL){
				BindingProtein* prot2 = myInteractingProt;
				prot2->setInteractingProt(NULL);
				this->myInteractingProt = NULL;
			}
		cout << endl; //TTT
		}
	}
	/* XXX do something: the maximum energy needs to more than the maximum from the DNA seq... otherwise the
	   probability gets stuck at prob=1 */
	if( prob >=1 ){
		prob = 0.99; // had at 0.999 (or higher) but could take up a lot of steps: maybe 1/3 of all the bound counts?
	}
	return prob;
}


