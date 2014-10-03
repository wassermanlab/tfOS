//
// ProteinBindingSimulator.cpp
//
// Author:
//	David Arenillas (dave@cmmt.ubc.ca)
//
// Copywrite:
//	Wasserman Lab
//	Centre For Molecular Medicine and Therapeutics
//	Child & Family Research Institute
//	University of British Columbia
// 
#include <iostream>
#include <fstream>
#include <cstdlib>  // std::rand; std::srand
#include <ctime>
#include <math.h>
#include <string>
#include <algorithm>  // std::random_shuffle
#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/ProteinBindingSimulator.h"

const unsigned int  ProteinBindingSimulator::DFLT_FALLOFF_DISTR = 2; // uniform = 1; gaussian = 2;
const unsigned int  ProteinBindingSimulator::DFLT_SEQ_ACCESSIBILITY = 1; // anything >0 is considered open to binding  
const pair<pair<float,string>, pair<double,double> > ProteinBindingSimulator::DFLT_NOINTERACTION = make_pair(make_pair(0,""), make_pair(0,0)); // anything >0 is considered open to binding 


// Constructor
ProteinBindingSimulator::ProteinBindingSimulator( const string &seq)
	: mySeq(seq)
{
    // initialize the occupied array to size of sequence
    myOccupied.assign(seq.size(), NULL);
	myFallOffDistr = DFLT_FALLOFF_DISTR; // default is uniform distribution (applies to all proteins)
	mySeqAccess = NULL;
	myInteractions = NULL;
	myMaxCoopDists = 0;
}

// Destructor
ProteinBindingSimulator::~ProteinBindingSimulator()
{
}

// Set sequence
void ProteinBindingSimulator::sequence(const string &seq)
{
    if (!&seq) return;

    mySeq = seq;
    // initialize the occupied array to size of sequence
    myOccupied.assign(seq.size(), NULL);
}

// Set fall off distribution for sliding protein
void ProteinBindingSimulator::fallOffDistribution(unsigned int distributionID)
{
	myFallOffDistr = distributionID;
}

// Set the max of all cooperative binding distances (used as a first pass threshold for interactions)
void ProteinBindingSimulator::maxAllCoopDistances(unsigned int coopDist)
{
	myMaxCoopDists = coopDist;
}


//	Add a new protein. We can have mulitple instances of the same protein as
//	well as different proteins.
bool ProteinBindingSimulator::addBindingProtein(BindingProtein *prot)
{
    if (!prot) return false;

    if (prot->forwardEnergyProfile().size() != (mySeq.size() - prot->length() + 1) ) {
    	cerr << "Error (ProteinBindingSimulator.cpp):  energy profile size ("
			 << prot->forwardEnergyProfile().size() << ", "<< prot->reverseEnergyProfile().size()<< ") does not concur with sequence length (" 
			 << mySeq.size() - prot->length() + 1 << ")" << endl;
		return false;
    }
	
    myProteins.push_back(prot);
	
	// initialize the occupancy counts for this protein
    // XXX counts for each individual protein or counts for each type of
    // protein? XXX
    if (!myTetheredCounts.count(prot->uniqueID()) ) {
		vector<unsigned int> tcounts;
		tcounts.assign(mySeq.size(), 0);
		myTetheredCounts[prot->uniqueID()] = tcounts;
	}

    if (!myBoundCounts.count(prot->uniqueID()) ) {
		vector<unsigned int> bcounts;
		bcounts.assign(mySeq.size(), 0);
		myBoundCounts[prot->uniqueID()] = bcounts;
    }

    return true;
}

bool ProteinBindingSimulator::addAccessibilityObj(SequenceAccessibility *access)
{
	if(!access) return false;
	mySeqAccess = access; // so myInteractions needs defining as a pointer in .h file 
	return true;
}

bool ProteinBindingSimulator::addInteractionsObj(ProteinInteractions *interac)
{
	if(!interac) return false;
	myInteractions = interac; // so myInteractions needs defining as a pointer in .h file 
	return true;
}

//
// Perform simulation with the given number of iterations.
// Return true on success, false on failure.
//
bool ProteinBindingSimulator::simulate(unsigned int iters ) //iters is from sim_initiate.cpp
{
	SequenceAccessibility *seqAccess = mySeqAccess;
	ProteinInteractions *interact = myInteractions;
    const unsigned int seqLen = mySeq.length();
    unsigned int numProteins = myProteins.size();
	int randIndex;
	vector<int> randOrderProt; // flags proteins as processed or not
    // XXX should use better seed (finer granularity than 1 sec)
    srand((unsigned) time(NULL));  // used for order of proteins and for protein position on DNA

	cout << "Beginning simulation..."<<endl;
	//cout << "fallOff distr: "<< myFallOffDistr <<endl;
	for (register unsigned int it = 0; it < iters; it++) { // iters defined in sim_initiate()
		//finiProteins = myProteins;
		randOrderProt.clear();
		for(unsigned int j=0; j<numProteins; j++){
			randOrderProt.push_back( j ); // get 0 through to numProteins-1
		}
		random_shuffle( randOrderProt.begin(), randOrderProt.end() );

		// cycle through all proteins in random order for each iteration
		for (register unsigned int pcount = 0; pcount < numProteins; pcount++) {	
			if( randOrderProt[ pcount ] == -1 ){//when collision occured a -1 removed a secondary protein ahead of time
				continue;
			}
			BindingProtein *prot = myProteins[ randOrderProt[pcount] ];
	    	const unsigned int protLen = prot->length(); // no idea why of all prot values I stored this one
			// to determine if seq is accessible to being bound (mean of values at protein footprint) 
			double seqaccess = DFLT_SEQ_ACCESSIBILITY; // assume all sequence accessible unless seqAccessibility data provided 
			if( mySeqAccess != NULL){	
				if( ! prot->isFree() ) {
					seqaccess = seqAccess->getAccessibility( prot->position(), prot->length());
				}
			}
			// interaction weight is used to determine contribution of prot::prot interaction to bound state stability
			// energy is to compare against the interacting protein for the strongest binding prot of the pair
			pair<pair<float,string>, pair<double,double> > weightEnergy = DFLT_NOINTERACTION;
			if( !prot->isFree() && myInteractions != NULL ) {
				weightEnergy  = isInteraction( prot, interact );
			}

	    	// ***************** FREE **************************************
	    	if (prot->isFree()) {
	    		// Protein is not bound or tethered to the DNA.
	    		// hop locally, or jump distally ?
				/* XXX may need to set hop "distance" to a range... or need to set the range 
						threshold as a variable somewhere */
				int newPos = translocate( prot, seqLen); // hop p=70%, jump p=30%
				
				//int newPos = rand() % (seqLen - protLen); // original choice of random position
				
				if (newPos != -1 && isOccupiedRegion( newPos, newPos + protLen - 1) == NULL) {
					int newStrand = rand() % 2; // value 0 or 1
					if (newStrand == 0) 
						newStrand = -1;
					// test if bound... if not, then it will tether
//					cout << "\t\t-- test "<<prot->uniqueID() <<" freeToBound\t";
					if (testProb(prot->freeToBoundProb(newPos))) {
						prot->bindAt(newPos, newStrand);
//						cout << prot->uniqueID() << " bound at " << newPos << endl;
						setOccupiedRegion(newPos, newPos + protLen - 1, prot);
						// XXX update bound counts here???
						incBoundCounts(prot->uniqueID(), newPos, newPos + protLen - 1);
		    		} 
		    		else {     // it is tethered (testProb( prot->freeToTetheredProb( newPos))) 
						int newDir = rand() % 2;
						if (newDir == 0) 
							newDir = -1;
						prot->tetherAt(newPos, newDir, newStrand);
//						cout << prot->uniqueID() << " tethered at " << newPos << endl;
						setOccupiedRegion(newPos, newPos + protLen - 1, prot);
						// XXX update tethered counts here???
						incTetheredCounts(prot->uniqueID(), newPos, newPos + protLen - 1);
		    			// XXX must convert prob. from multi-nomial to binomial
					} 
				} else {     // remain free when pos. is occupied
					prot->freeOrigin(newPos); // either -1 or the place that was occupied by another
				//	if(newPos > -1)
//	 	    			cout << " region " << newPos << "-" << newPos + protLen - 1 << " occupied" << endl;
//	 	    		cout << prot->uniqueID() << " remains free " << endl;
	 	    	}
	 	    // ***************** TETHERED ***********************************
			} else if (prot->isTethered()) {
	    		// Protein is tethered to (i.e. sliding along) the DNA.
//	    		cout << "\t\t-- test "<<prot->uniqueID() <<"  tetheredToFree\t";
				if ( testProb(prot->tetheredToFreeProb( myFallOffDistr )) ) { //distribution from user via simulate.cpp 
					// protein's next step is to "fall off"
					cout << prot->uniqueID() << " fell off\t" << prot->slidCount() << endl; //TTT
					setUnoccupiedRegion(prot->position(), prot->position() + protLen - 1);
		    		prot->fallOff();
				} else {
					// protein slides (right or left) to next DNA position and assesses state
					int curPos = prot->position();
		  	  		int curDir = prot->direction();
					int newStartPos = curPos + curDir;
					int newEndPos = newStartPos + prot->length() - 1;
					
					// move to the right
		    		if (prot->direction() == 1) {
						if (newEndPos >= (int) mySeq.size()) {
			    			// slide off end of DNA
			    			setUnoccupiedRegion(curPos, newEndPos - 1);//why is this not 'til mySeq.size() ?
			    			prot->fallOff();
//			    			cout << prot->uniqueID() << " slid off end" << endl;
						} else if (isOccupiedPosition(newEndPos) != NULL) {
							// collision
							BindingProtein* prot2 = isOccupiedPosition(newEndPos);
//			    			cout << prot->uniqueID() <<" "<<curPos+protLen-1<< " collision! with: " <<
			    //					prot2->uniqueID() <<" "<<prot2->position() << endl;
					    	setUnoccupiedRegion(curPos, newEndPos - 1);
					    	prot->fallOff();
					    	if(prot2->isTethered()) {
								setUnoccupiedRegion(prot2->position(), prot2->position()+prot2->length()-1);
								prot2->fallOff();
								/* locate the element in randOrderProt[] with the index for prot2 and 
								     insert -1 so that if it is still ahead in the list it wont be used.
								   Even if the 2nd protein has already been "done" in this iteration
								   it has to unoccupy DNA now, otherwise in next iteration it may be
								   interfere with the other two if it is ordered after them */
								/*for( int kount=0; kount < (int) numProteins; kount++) {
									if(prot2->uniqueID().compare(myProteins[kount]->uniqueID()) == 0) {
//										cout << prot2->uniqueID() << " fell off - collision" << endl;
										// replace the element containing the index 'kount' with -1
										std::replace(randOrderProt.begin(), randOrderProt.end(), kount, -1);
										break;
									}
								}*/
			    			}
						} else {
			    			// next position is clear, slide there
			    			prot->slideCurrent(); // sets only the new position, and increments count
//							cout << prot->uniqueID() << " slid right" << endl;
							setOccupiedPosition(newEndPos, prot); // these two lines do the work of
			    			setUnoccupiedPosition(curPos);       // re-setting the occupied region
			    			
			    			/* XXX set the STATE of this new position i.e. tethered vs bound
			    			   (we already know it doesn't fall or collide) */
//			   				cout << "\t\t-- test "<<prot->uniqueID() <<"  tetheredToBound\t";
			    			if (testProb(prot->tetheredToBoundProb(seqaccess, weightEnergy))) {
		    					// State -> Bound
								//  interaction with another protein could cause transition to binding
		    					prot->bindAt(prot->position());
		    					setOccupiedRegion(prot->position(), prot->position()+protLen-1, prot);
		    					cout << "    " << prot->uniqueID() << " t-bound at " << prot->position() << endl; //TTT
							cout <<"T------------------------------"<<endl;//TTT
		    					// XXX update bound counts here???
		    					incBoundCounts(prot->uniqueID(), prot->position(), prot->position() + protLen - 1);
								// XXX must convert prob. from multi-nomial to binomial
								// if an interaction existed then give 2nd protein 'bound' status
								if( prot->isInteractingProt() != NULL ){
									BindingProtein* prot2 = prot->isInteractingProt();
									prot2->bindAt(prot2->position());  //interacting protein has acquired a bound state too
									incBoundCounts(prot2->uniqueID(), prot2->position(), prot2->position()+prot2->length()-1);
								}
							} else {
								// State = tethered and the strand and direction stay the same
								incTetheredCounts(prot->uniqueID(), newStartPos, newEndPos);
								// remove tentative interactions that just failed testProb()
								if( prot->isInteractingProt() != NULL ){
									BindingProtein* prot2 = prot->isInteractingProt();
									prot2->setInteractingProt(NULL);
									prot->setInteractingProt(NULL);
								}
			    			}
						}
		    		} 
		    		else { // move to the left
						if (newStartPos < 0) {
			    			// slide off beginning of DNA
			    			setUnoccupiedRegion(0, newEndPos + 1);
			    			prot->fallOff();
//			    			cout << prot->uniqueID() << " slid off beginning" << endl;
						} else if (isOccupiedPosition(newStartPos) != NULL) {
							// collision
							BindingProtein* prot2 = isOccupiedPosition(newStartPos);
//							cout << prot->uniqueID()<<" "<<curPos<<" collision! with: "<<prot2->uniqueID()
//								<<" "<<prot2->position()+prot2->length()-1<<endl;
							setUnoccupiedRegion(curPos, newEndPos + 1);
							prot->fallOff();
							if(prot2->isTethered()) {
								setUnoccupiedRegion(prot2->position(), prot2->position()+prot2->length()-1);
								prot2->fallOff();
								/* locate the protein in randOrderProt and mark "done"
								   Even if the 2nd protein has already been "done" in this iteration
								   it has to unoccupy DNA now, otherwise in next iteration it may be
								   interfere with the other two if it is ordered after them */
								/*for( int kount=0; kount < (int) numProteins; kount++) {
									if(prot2->uniqueID().compare(myProteins[kount]->uniqueID()) == 0) {
//										cout << prot2->uniqueID() << " fell off - collision" << endl;
										std::replace(randOrderProt.begin(), randOrderProt.end(), kount, -1);
										break;
									}
								}*/
			    			}
						} else {
			    			// next position is clear, slide there
			    			prot->slideCurrent();
//			    			cout << prot->uniqueID() << " slid left" << endl;
			    			setOccupiedPosition(newStartPos, prot);
			    			setUnoccupiedPosition(newEndPos + 1);
			    			
			    			/* XXX set the State of this new position i.e. tethered vs bound
			    			   (we already know it doesn't fall or collide) */
//			    			cout << "\t\t-- test "<<prot->uniqueID() <<"  tetheredToBound\t";
							
							if (testProb(prot->tetheredToBoundProb(seqaccess, weightEnergy))) {
		    					// protein binds DNA at current tethered position on current strand
								//  interaction with another protein could cause transition to binding
								prot->bindAt(prot->position());
		    					setOccupiedRegion(prot->position(), prot->position()+protLen-1, prot);
		    					cout << "  " << prot->uniqueID() << " t-bound at " << prot->position() << endl; //TTT
							cout <<"T------------------------------"<<endl;//TTT
		    					incBoundCounts(prot->uniqueID(), prot->position(), prot->position()+protLen-1);
								// XXX must convert prob. from multi-nomial to binomial
								// if an interaction existed then give 2nd protein 'bound' status
								if( prot->isInteractingProt() != NULL ){
									BindingProtein* prot2 = prot->isInteractingProt();
									prot2->bindAt(prot2->position()); //interacting protein was presumably already bound, but making sure
									incBoundCounts(prot2->uniqueID(), prot2->position(), prot2->position()+prot2->length()-1);
								}
							} else {
								// State = tethered and the strand and direction stay the same
								// XXX tethered counts
			    				incTetheredCounts(prot->uniqueID(), newStartPos, newEndPos);
								// remove tentative interactions that just failed testProb()
								if( prot->isInteractingProt() != NULL ){
									BindingProtein* prot2 = prot->isInteractingProt();
									prot2->setInteractingProt(NULL);
									prot->setInteractingProt(NULL);
								}
			    			}						
						}
					}
				}
			// ******************* BOUND *****************************************
	    	} else if (prot->isBound()) {
//	    		cout << "\t\t-- test "<<prot->uniqueID() <<"  remainBound, then test boundToFree\t";
                // XXX must convert prob. from multi-nomial to binomial
				if(testProb(prot->remainBoundProb(seqaccess, weightEnergy)) ){
					cout << prot->uniqueID()<<"remain Bound "<<  endl; //TTT				
					cout << "B------------------------------"<<endl;//TTT
//					cout << prot->uniqueID() << " remains bound at " << prot->position() << endl;
					incBoundCounts(prot->uniqueID(), prot->position(), prot->position()+protLen-1);
					if( prot->isInteractingProt() != NULL ){
						BindingProtein* prot2 = prot->isInteractingProt();
						prot2->bindAt(prot2->position());  //interacting protein has acquired a bound state too
						incBoundCounts(prot2->uniqueID(), prot2->position(), prot2->position()+prot2->length()-1);
					}
				} else if( testProb(prot->boundToFreeProb()) ) { // isSeqAccess(pos, prot->length(), seqAccess) 
					setUnoccupiedRegion(prot->position(), prot->position() + protLen - 1);
					// if one protein in an interaction breaks off, both become free
					if( prot->isInteractingProt() != NULL ){
						BindingProtein* prot2 = prot->isInteractingProt();
						prot2->fallOff();
					}
		    		prot->fallOff();
//		    		cout << prot->uniqueID() << " fell off" << endl;
				} else { // enter tethered state. Currently never happens after being bound
					// won't tether if was in an interaction
					if( prot->isInteractingProt() != NULL ){
						BindingProtein* prot2 = prot->isInteractingProt();
						prot2->fallOff();
						prot->fallOff();	
						continue; // back to protein loop
					}
					// tether if was not interacting with another protein
					int newDir = rand() % 2;
					if (newDir == 0)
						newDir = -1;
		    		int newStrand = rand() % 2;
		    		if (newStrand == 0) 
		    			newStrand = -1;
		    		prot->tetherAt(prot->position(), newDir, newStrand);
//		    		cout << prot->uniqueID() << " tethered at " << prot->position() << endl;
		    		// XXX tethered counts
		    		incTetheredCounts(prot->uniqueID(), prot->position(), prot->position()+protLen-1);
				} 
			// ******************** ERROR ***************************************
	    	} else {
	    		cerr << "ERROR: protein in undefined state" << endl; 
	    	}
		}
	}
	cout << "End of simulation."<<endl;
    return true;
}

//******************************************************************************
//
// Private methods follow
//
//******************************************************************************

inline bool ProteinBindingSimulator::testProb(float prob) const
{	
	// random seed is initialized to a value representing the second in which the program is executed
    double testStat = (double) rand() / RAND_MAX; /* rand() is number from 0 to RAND_MAX
    						  			                RAND_MAX is a C++ constant >= 32767 */
    
	//XXX need to check why this is called twice for a TF sometimes
	if(prob>0.8){ //TTT
		cout << "PBS: line 418 testProb: " << testStat << " prob: "<<prob<<endl; //TTT
	} //TTT
	return testStat < prob; 
}

/*// get accessibility value of sequence (used to decide if TF can bind)
double ProteinBindingSimulator::isSeqAccess(BindingProtein *prot, SequenceAccessibility *seqAccess) const
{
	if( seqAccess->accessDataExists() ){
	cout << "PBS: line 366"<<endl; //TTT
		return seqAccess->getAccessibility(prot->position(), prot->length());
	}else{
		return seqAccess->noCoverageVal(); // no coverage is likely set as -1
	}
}*/

// get the weight of interaction and 2nd protein's energy of binding 
pair<pair<float,string>, pair<double, double> > ProteinBindingSimulator::isInteraction(BindingProtein *prot1, ProteinInteractions *inter ) const
{
	BindingProtein *prot2;
	pair <pair<float, int>, string> interWghtDistState;
	// weight=0 signals no interaction, other values not considered
	pair<pair<float, string>, pair<double, double> > wghtNrg = DFLT_NOINTERACTION;
	int pos1, pos2;
	int leftEnd, rightEnd;

	pos1 = prot1->position();
	rightEnd = pos1+prot1->length()+myMaxCoopDists;
	leftEnd =  pos1-myMaxCoopDists;
	if( rightEnd > (int) mySeq.size()-1 ){
		rightEnd = (int) mySeq.size()-1; // last position
	}
	if(leftEnd <1 ){
		leftEnd=1; //first position	
	}
	prot2 = isOccupiedRegion( pos1+prot1->length(), rightEnd ); // is there a protein somewhere close to the right?
	if( prot2 == NULL){ // check other direction
		prot2 = isOccupiedRegion( leftEnd, pos1-1 ); // is there a protein somewhere close to the left?
	}
	//XXX should check whether either protein is already in an interaction 'cause at the moment we won't allow
	if(prot2 !=NULL){
		// check 1) proteins can interact
		//       2) close enough to interact
		//       3) at least one is bound (bindingProbability() in BindingProtein.cpp will require both to be above bindingThreshold) 
		pos2 = prot2->position(); // position of the near-by protein
		interWghtDistState = inter->getInteractionValues(prot1->name(), prot2->name() );
		float weight = (interWghtDistState.first).first;
		int dist = (interWghtDistState.first).second;
		string state = interWghtDistState.second;
		// if either protein is bound return weight so that interaction can be further assessed by bindingProbability() in BindingProtein.cpp
		if(  weight != 0 && (prot1->isBound() || prot2->isBound()) && abs(pos1-pos2) <= dist ) { // is the 2nd protein close enough to interact
			cout <<"PBS: line 424 if either one bound state: "<< prot1->uniqueID()<<" could interact "<<prot2->uniqueID()<<endl;//TTT
			//XXX need to include state
			wghtNrg = make_pair( make_pair(weight, state), make_pair(prot2->energyAt(pos2), prot2->thresholdBindingEnergy() ) );
			// set potential interaction, to be erased in BindingProtein bindingProbability() if fails requirements
			// then has to be checked again if testPro() does not ok binding for tetherToBound or remainBound
			prot1->setInteractingProt(prot2);
			prot2->setInteractingProt(prot1);
			return wghtNrg;
		}
	}
	return wghtNrg;
}

BindingProtein* ProteinBindingSimulator::isOccupiedPosition(unsigned int pos) const
{
    return myOccupied[pos];
}

inline void ProteinBindingSimulator::setOccupiedPosition(unsigned int pos, BindingProtein* prot)
{
    myOccupied[pos] = prot;
}

inline void ProteinBindingSimulator::setUnoccupiedPosition(unsigned int pos)
{
    myOccupied[pos] = NULL;
}

BindingProtein* ProteinBindingSimulator::isOccupiedRegion( unsigned int start, unsigned int end) const
{
    BindingProtein* occupied = NULL;
	// in simulate() make sure start and end are within sequence bounds
	for (register unsigned int pos = start; pos <= end; pos++) {
    	if (myOccupied[pos] != NULL) {
	    	occupied = myOccupied[pos];
	    	break;
		}
    }
    return occupied;
}

void ProteinBindingSimulator::setOccupiedRegion( unsigned int start, unsigned int end, BindingProtein* prot)
{
    for (register unsigned int pos = start; pos <= end; pos++)
    	myOccupied[pos] = prot;
}

void ProteinBindingSimulator::setUnoccupiedRegion( unsigned int start, unsigned int end)
{
    for (register unsigned int pos = start; pos <= end; pos++)
    	myOccupied[pos] = NULL;
}


// XXX SHould this even be here... perhaps it should be in BindingProtein.cpp
int ProteinBindingSimulator::translocate( BindingProtein* prot, unsigned int seqLen)
{
	// return new position from hop or jump
	if( prot->freeOrigin() != -1 && prot->maxHopDistance() > 0 ) { //jump if it diffused off the DNA earlier
		if( testProb( prot->hopnotjumpProb() ) ) { //TRUE=common event ie hop (roughly 70% of decisions). 
			// HOP (myHopProb p=0.70 default)
			int hopDir = rand() % 2;
			if (hopDir == 0)   hopDir = -1;
			/* XXX inversely relate hop distance to proteinLength? since presumably diffusion
					is related to size, and the footprint might tell us a protein is big or not? */
			int newPos = prot->freeOrigin() + (rand() % (prot->maxHopDistance() ))*hopDir;
			if( newPos < 0 || newPos > int(seqLen - prot->length()) ) {
				return -1; // hopped off the DNA
			}
			else {
				return newPos; // hopped to a location
			}
		}else{
			// jump. p=30%: allow simulator to randomly choose initial position
			// not "blanking out" local DNA, so it could jump local
		}
	}else{
		// jump if it had already left the DNA in prior step
		//  set some random initial position to test DNA.
	}
	return rand() % (seqLen - prot->length());		
}

// XXX update full region or just start position (start position would also need record a direction)
// XXX Full region would mimic ChIP footprint 
void ProteinBindingSimulator::incTetheredCounts( const string &uniqueID, unsigned int start, unsigned int end)
{
    for (register unsigned int pos = start; pos <= end; pos++)
    	myTetheredCounts[uniqueID][pos]++;
}

// XXX update full region or just start position (start position would also need record a direction)
// XXX Full region would mimic ChIP footprint
void ProteinBindingSimulator::incBoundCounts( const string &uniqueID, unsigned int start, unsigned int end)
{
    for (register unsigned int pos = start; pos <= end; pos++)
    	myBoundCounts[uniqueID][pos]++;
}

