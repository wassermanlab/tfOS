// launch simulator:  ./sim_initiate.cpp
// 2014.07.17 hopefully have allowed input and output file names to be provided by user

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/ProteinBindingSimulator.h"
#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/ProteinInteractions.h"
#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/SequenceAccessibility.h"
#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/BindingProtein.h"

using namespace std;


// sequence length should be no more than 10,000
/* previously we had 60,000 simulate objects (each with several proteins), and the user could say
how many jobs to run i.e. 10 jobs, means 10 sets of 60,000 */
	 
//static const char PATH_INPUT_PARAMETERS[] = "/raid6/rebecca/project_SIMULATION/occupancySim_C++/INPUT_parameters.txt"; 
//static const char PATH_OUTPUT_FILE[] =		"/raid6/rebecca/project_SIMULATION/occupancySim_C++/simulation/simOut_count_file.txt";

static const unsigned int	NUM_ITER = 60000; // how often simulate is called i.e. protein steps
static const unsigned int	DFLT_PROT_COPYNUMBER = 2; 

ProteinBindingSimulator* read_input(const char*);  

int main(int argc, char* argv[])
{
	if( argc != 3 && argc != 4 ){
        cout << "Missing arguments:" << endl;
        cout << "./sim_initiate  arg1  arg2" << endl;
        cout << "arg1  /absolute path/INPUT_parameters.txt" << endl;
        cout << "arg2  /absolute path/name_simulation_output_file.txt" << endl; 
        cout << "arg3  (optional) number of iterations (default 60,000)" << endl; 
        exit(1);
    }
    const char* infilename = argv[1];  // parameter file to run simulation
    const char* outfilename = argv[2];
	//std::string outfilename = argv[2]; # can be modified apparently
	unsigned int numIters = NUM_ITER;
	if( argc == 4){
		char* p;
		unsigned int i = strtol(argv[3], &p, 10);
		if( *p == '\0') { // got to the end of the argument (which means they were numbers as strtol checks in some way
			numIters = i;
		}
	}

    //Just for info, print the time simulator starts and ends
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
	// check if outfile was created, if so then keep going
	// the read_input function will check if infile can open
	ofstream myfileOut;  //ofstream myfileOut (PATH_OUTPUT_FILE);
	myfileOut.open(outfilename);
 	if( myfileOut.is_open()) {
    	cout << "Start sim_initiate: " << asctime (timeinfo);
		ProteinBindingSimulator *sim = read_input( infilename); // add seq + proteins etc from perl file
        // cout << "coop dist: " << sim->cooperativeDistance() << endl;
        // cout << "sequence: " << sim->sequence() << endl;
        sim->simulate( numIters ); //this is one call i.e. one simulation where proteins interact for 60,000 steps (iterations) 
        
		// print sim. results to file
        //XXX sending to both simOut_count_file.txt and to standard out for now until figure out how to put separate jobs into
        // seperate files.
        if (myfileOut) {
            for (unsigned int k = 0; k < sim->bindingProteins().size() ; k++) {
                //myfileOut << "#proteinID\tposition\tCounts\tstate" << endl;
                const BindingProtein *prot1 = sim->bindingProteins()[k]; //the first protein in list
                for (unsigned int i = 0; i < sim->proteinTetheredCounts(prot1->uniqueID()).size(); i++) {
                    myfileOut <<prot1->uniqueID()<<"\t"<< i <<"\t"<<sim->proteinTetheredCounts(prot1->uniqueID())[i]<<"\ttethered" <<endl;
                }
                //myfileOut << endl << "proteinID\tposition\tCounts\tstate" << endl;
                for (unsigned int i = 0; i < sim->proteinBoundCounts(prot1->uniqueID()).size(); i++) {
                    myfileOut <<prot1->uniqueID()<<"\t"<< i <<"\t"<<sim->proteinBoundCounts(prot1->uniqueID())[i]<<"\tbound" <<endl;
                }
            }
        }
		myfileOut.close();
    }else{
		cout << "Error: Unable to create outfile - " << outfilename << endl;
	    return 0;
	}

    time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	cout << endl << "End sim_initiate: " << asctime (timeinfo) << endl;
    
	return 0;
}

ProteinBindingSimulator* read_input(const char* filename )
{
    const int MAX_BUF = 10000+100; //big enough for seqs. Added 100 to buffer the buffer
    const int MAX_SEQ = 10009; // allowing sequence of length 10,000 (allow a few nucleotides miscount)
    const int MAX_PROT_NAME = 30;
    char buf[MAX_BUF];
	int protCopyNumber = DFLT_PROT_COPYNUMBER;
    ProteinBindingSimulator *sim = NULL;
    ProteinInteractions* interac = NULL;
    SequenceAccessibility* access = NULL;
    BindingProtein* prot = NULL;
	BindingProtein* copyprot = NULL;
    vector<bp_energy_t> *fep = NULL;
    vector<bp_energy_t> *rep = NULL;
	
	//ifstream myinfile(PATH_INPUT_PARAMETERS);
	ifstream myinfile(filename);
    if (myinfile.is_open()) {
        int skipgetline = 0;
		while (!myinfile.eof()) {
            if(skipgetline == 0){ // at one or file lines I set skipgetline=1 so that buf retains the line header
				myinfile.getline(buf, MAX_BUF);
            }
			skipgetline=0; //reset back to 0 everytime as should only have value 1 for one loop
			if (strncmp(buf, "seq=", 4) == 0) {  // occurs only once for the simulation
                char seqstr[MAX_SEQ];
                sscanf(buf, "seq=%s", seqstr);
                sim = new ProteinBindingSimulator(seqstr);
            } else if (strncmp(buf, "maxallcoopdist=",15) == 0) { // occurs once; may need to do pairs of proteins sometime 
                unsigned int coopDist;
                sscanf(buf, "maxallcoopdist=%d", &coopDist);
                sim->maxAllCoopDistances( coopDist );
            } else if (strncmp(buf, "falloff_distribution=",19) == 0) { // occurs once for simulation (applies to all proteins) 
                unsigned int fallOffDistr;
                sscanf(buf, "falloff_distribution=%d", &fallOffDistr);
                sim->fallOffDistribution( fallOffDistr );
			} else if (strncmp(buf, "epigendensitythreshold=",24) == 0) { // occurs once for simulation (applies to all proteins) 
                // NOT USING YET. At moment code relies on epigentic scores being >0 (in either active or repressive context)
				unsigned int histThr;                                      
                sscanf(buf, "epigendensitythreshold=%d", &histThr);
                //sim->epigenThreshold( histThr );  // no function yet
            } else if (strncmp(buf, "coop_interactions", 18) == 0) {
				myinfile.getline(buf, MAX_BUF);
				if( strncmp(buf, "end_interactions", 17) == 0 ){
					continue;
				}
				interac = new ProteinInteractions();
				while( strncmp(buf, "end_interactions", 17) != 0){
					// prior computeBindingEnergies.pl code ensured 4 elements to a line and that ordering was consistent, so don't need to check
					string line(buf);
					string prot1, prot2, strweight, strdist, state;
					istringstream ss(line);  //line becomes an input stream; Then each direction operator >> directs the next element in the stream
										// can then get each word in line; as character based streams are deemed to be separated by whitespace \n,\t,' '
					ss >> prot1 >> prot2 >> strweight >> strdist >> state;
					cout << "protein interaction added: "<<prot1<<" "<<prot2  //TTT
						 << "      weight: "<<strweight<<" dist: "<<strdist<<" state: "<<state <<endl; //TTT
					interac->setInteractionVal(prot1, prot2, atof(strweight.c_str()), atoi(strdist.c_str()), state);
					myinfile.getline(buf, MAX_BUF);
				}
				sim->addInteractionsObj(interac);
			  	skipgetline=1; // retain "end_interactions" in buf so that outer while loop
            } else if (strncmp(buf, "protein=", 8) == 0) {  
                char protName[MAX_PROT_NAME];
                sscanf(buf, "protein=%s", protName);
                // the recently populated protein object is added, before the next protein is generated 
                if (prot) {  
                    cout << "add " << protCopyNumber << " "<< prot->name() << " protein(s)\n";
                    // add shallow copies of protein to simulator (shallow ok, all fields are values not pointers)
                    if( protCopyNumber > 1) {
                        for (int num = 1; num < protCopyNumber; num++) {
                            std::ostringstream in;
                            in << num;// << std::ends;
                            std::string protnum = in.str();
                            copyprot = new BindingProtein(*prot); // copies of  original
                            copyprot->uniqueID(prot->name() + "-" + protnum); //e.g. SRF-1, SRF-2,
                            sim->addBindingProtein(copyprot);
                        }
                    }
                    prot->uniqueID(prot->name() + "-0"); // e.g. SRF-0
                    sim->addBindingProtein(prot);
                }
                // create the next protein, which needs to be populated
                prot = new BindingProtein();
                prot->name(protName);
                protCopyNumber = DFLT_PROT_COPYNUMBER;
            } else if (strncmp(buf, "copy_number=", 12) == 0) {
                //protCopyNumber defines loop iterations (above) adding protein objects
                sscanf(buf, "copy_number=%d", &protCopyNumber);
    //			cout << prot->name() << " copynum "<< protCopyNumber<<endl;
            } else if (strncmp(buf, "length=", 7) == 0) {
                int protLen;
                sscanf(buf, "length=%d", &protLen);
                prot->length(protLen);
            } else if (strncmp(buf, "pioneer=", 8) == 0) {
				bool pioneer;
                sscanf(buf, "pioneer=%d", &pioneer);
            	prot->isPioneer( pioneer);
			} else if (strncmp(buf, "sliding_distance=", 17) == 0) {
                unsigned int protMaxSlide;
                sscanf(buf, "sliding_distance=%d", &protMaxSlide);
            	prot->maxSlide( protMaxSlide);
            } else if (strncmp(buf, "maxhopdistance=",15) == 0) { // occurs once for simulation (applies to all proteins) 
                unsigned int maxHopDist;
                sscanf(buf, "maxhopdistance=%d", &maxHopDist);
                prot->maxHopDistance( maxHopDist );
            } else if (strncmp(buf, "hopnotjumpChance=",17) == 0) { // occurs once for simulation (applies to all proteins) 
                unsigned int hopChance;
                sscanf(buf, "hopnotjumpChance=%d", &hopChance);
                prot->hopnotjumpProb( hopChance );
            } else if (strncmp(buf, "max_energymatrix=", 17) == 0) {
                double protMaxEnergy;
                sscanf(buf, "max_energymatrix=%lf", &protMaxEnergy);
                prot->maxEnergy((bp_energy_t) protMaxEnergy);
            } else if (strncmp(buf, "threshold_energy=", 17) == 0) {
                double protThresholdEnergy;
                sscanf(buf, "threshold_energy=%lf", &protThresholdEnergy);
                prot->thresholdBindingEnergy((bp_energy_t) protThresholdEnergy);
            } else if (strncmp(buf, "palindromic=", 12) == 0) { // specifies if protein has orientation
                char c;
                sscanf(buf, "palindromic=%c", &c);
                bool protPalindromic; // BindingProtein.cpp default is false
                if (c == 'y') { 	  // thus only need to worry if user set "true"
                    protPalindromic = true;
                    prot->isPalindromic(protPalindromic);
                }	
            } else if (strncmp(buf, "forward_energies:", 17) == 0) {
				// add energy one by one
				myinfile.getline(buf, MAX_BUF);    
            	while( strncmp(buf, "reverse_energies:", 17) != 0 ){
					prot->addForwardEnergy((bp_energy_t) strtod(buf, NULL));
					myinfile.getline(buf, MAX_BUF);    
				}
				skipgetline=1;
			} else if (strncmp(buf, "reverse_energies:", 17) == 0) {
				char* bufptr;
				myinfile.getline(buf, MAX_BUF);    
				bufptr = buf; // need char* to use strlen() to make sure not empty line
				while(strlen(bufptr)>0 && strncmp(buf, "PROTEIN_SPECIFIC:", 17) !=0 && !myinfile.eof() ) {
					prot->addReverseEnergy((bp_energy_t) strtod(buf, NULL));	
					myinfile.getline(buf, MAX_BUF);    
					bufptr = buf;
				}
				skipgetline=1;
		  	}else if( strncmp(buf, "active_epigeneticScores:", 25) == 0 ){
		  		// put histone marks into vector <double> then send
				// like energies, each hist mod position is on a line
				char* bufptr;
				myinfile.getline(buf, MAX_BUF);    
				bufptr = buf; // need char* to use strlen() to make sure not empty line
				if(strlen(bufptr)==0 || strncmp(buf, "END_PARAMETER_FILE:", 19) ==0 && myinfile.eof() ) {
					continue;
				}
				access = new SequenceAccessibility();
				while(strlen(bufptr)>0 && strncmp(buf, "END_PARAMETER_FILE:", 19) !=0 && !myinfile.eof() ) {
					access->addActiveEpigen((double) strtod(buf, NULL));	
					myinfile.getline(buf, MAX_BUF);    
					bufptr = buf;
				}
                sim->addAccessibilityObj(access);
				skipgetline=1;
				//was reading single line of histValues into a vector...but vector must be pass-by-ref so useless
				/*string hist;
				vector <double> histVal;
				getline(myinfile, hist);
				istringstream ss( hist );
				string hist2;
				while(ss>>hist2){
					double histNum;
					histNum = atof(hist2.c_str());
					// could send each value to sim to store like sequence (this is how energies are sent to BindingProtein)
					// or can store as a vector<double> and send to sim
					histVal.push_back(histNum);
				}
				access->setActiveHist( &histVal);*/
		  	}else if( strncmp(buf, "repress_epigeneticScores:", 26) == 0 ){
		  		// put histone marks into vector <double> then send
				// like energies, each hist mod position is on a line
				char* bufptr;
				myinfile.getline(buf, MAX_BUF);    
				bufptr = buf; // need char* to use strlen() to make sure not empty line
				if(strlen(bufptr)==0 || strncmp(buf, "END_PARAMETER_FILE:", 19) ==0 && myinfile.eof() ) {
					continue;
				}
				access = new SequenceAccessibility();
				while(strlen(bufptr)>0 && strncmp(buf, "END_PARAMETER_FILE:", 19) !=0 && !myinfile.eof() ) {
					access->addRepressiveEpigen((double) strtod(buf, NULL));	
					myinfile.getline(buf, MAX_BUF);    
					bufptr = buf;
				}
                sim->addAccessibilityObj(access);
				skipgetline=1;
		  	}else{
				//cout << "no match: "<<buf << endl;
				// do nothing	
			}
            if( myinfile.eof() ) {
                cout << "Last protein:\nAdd "<< protCopyNumber << " " << prot->name() << " protein(s)" <<endl; 
                if( protCopyNumber > 1) {
                    for (int num = 1; num < protCopyNumber; num++) {
                        std::ostringstream in;
                        in << num; // << std::ends;
                        std::string protnum = in.str();
                        copyprot = new BindingProtein(*prot);
                        copyprot->uniqueID(prot->name() + "-" + protnum); //e.g. the copies, SRF-1, SRF-2
                        sim->addBindingProtein(copyprot);
                    }
                }
                prot->uniqueID(prot->name() + "-0"); // e.g. the original, SRF-0
                sim->addBindingProtein(prot);
            	cout << "End of parameter file," << endl << "Proteins are populated." << endl;
            }
        } // end while
    } else{
        cout << "Error: Unable to open infile - " << filename << endl;
        return 0;
    }
	myinfile.close();
    return sim;
}

