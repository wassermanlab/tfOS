//sim_initiate.cpp

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>

#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/ProteinBindingSimulator.h"
#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/BindingProtein.h"

// previously we had 60,000 simulate objects (each with several proteins), and the user could say
// how many jobs to run i.e. 10 jobs, means 10 sets of 60,000
	 
static const char PATH_INPUT_PARAMETERS[] = "/raid6/rebecca/project_SIMULATION/occupancySim_C++/INPUT_parameters.txt"; 
static const char PATH_OUTPUT_FILE[] =		"/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/count_file.txt";

static const unsigned int	NUM_ITER = 60000; // how often simulate is called i.e. protein steps
static const unsigned int	DFLT_PROT_COPYNUMBER = 1; 
//static const char			DFLT_PALINDROMIC = 'n'; // not needed. Set in BindingProtein constructor

ProteinBindingSimulator* read_input();  

using namespace std;

int main()
{
	//Just for info, print the time simulator starts and ends
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	ProteinBindingSimulator *sim = read_input(); // add seq + proteins etc from perl file

//    cout << "coop dist: " << sim->cooperativeDistance() << endl;
//    cout << "sequence: " << sim->sequence() << endl;

    sim->simulate(NUM_ITER); //this is one call i.e. one simulation where proteins interact for X steps (iterations) 
	
	//XXX sending to both count_file.txt and to standard out for now until figure out how to put separate jobs into
	// seperate files.
	ofstream myoutfile (PATH_OUTPUT_FILE);
	if (myoutfile) {
		for (unsigned int k = 0; k < sim->bindingProteins().size() ; k++) {
			myoutfile << "Tethered Counts" << endl;
			const BindingProtein *prot1 = sim->bindingProteins()[k]; //the first protein in list
    		for (unsigned int i = 0; i < sim->proteinTetheredCounts(prot1->uniqueID()).size(); i++) {
				myoutfile <<prot1->uniqueID()<<" "<< i <<":\t"<<sim->proteinTetheredCounts(prot1->uniqueID())[i]<<endl;
			}
			myoutfile << endl << "Bound Counts" << endl;
			for (unsigned int i = 0; i < sim->proteinBoundCounts(prot1->uniqueID()).size(); i++) {
    			myoutfile <<prot1->uniqueID()<<" "<< i <<":\t"<<sim->proteinBoundCounts(prot1->uniqueID())[i]<<endl;
			}
		}
    }
	myoutfile.close();

    cout << "Start time and date: " << asctime (timeinfo);
    time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	cout << "End time and date: " << asctime (timeinfo) << endl;
    
	return 0;
}

ProteinBindingSimulator* read_input()
{
    const int MAX_BUF = 10000+100; //added 100 because seq sometimes DNA length > 9995 nt
    const int MAX_PROT_NAME = 30;
    const int MAX_SEQ = 10000; // using sequence of length 10,000
	
    char buf[MAX_BUF];
	int protCopyNumber = DFLT_PROT_COPYNUMBER;
    ProteinBindingSimulator *sim = NULL;
    BindingProtein* prot = NULL;
	BindingProtein* copyprot = NULL;
	
    vector<bp_energy_t> *fep = NULL;
    vector<bp_energy_t> *rep = NULL;
    
	ifstream infile(PATH_INPUT_PARAMETERS);
	
    int forwRev = 0; // just a flag to decide whether at forward or reverse energies 
    while (!infile.eof()) {
    	infile.getline(buf, MAX_BUF);
		if (strncmp(buf, "seq=", 4) == 0) {  // occurs only once for the simulation
			char seqstr[MAX_SEQ];
			sscanf(buf, "seq=%s", seqstr);
			sim = new ProteinBindingSimulator(seqstr);
		} else if (strncmp(buf, "fallOff_distribution=",19) == 0) { // occurs once for simulation (applies to all proteins) 
			unsigned int fallOffDistr;
			sscanf(buf, "fallOff_distribution=%d", &fallOffDistr);
			sim->fallOffDistribution( fallOffDistr );
		} else if (strncmp(buf, "cooperative_distance=",21) == 0) { // occurs once; may need to do pairs of proteins sometime 
			unsigned int coopDist;
			sscanf(buf, "cooperative_distance=%d", &coopDist);
			sim->cooperativeDistance( coopDist );
		} else if (strncmp(buf, "protein=", 8) == 0) {  
			char protName[MAX_PROT_NAME];
			forwRev = 0;
			sscanf(buf, "protein=%s", protName);
			// the recently populated protein object is added, before the next protein is generated 
			if (prot) {  
				cout << "add " << protCopyNumber << " "<< prot->name() << " protein(s)\n";
				// add shallow copies of protein to simulator (shallow ok, all fields are values not pointers)
				if( protCopyNumber > 1) {
					for (int num = 1; num < protCopyNumber; num++) {
						std::ostringstream in;
						in << num << std::ends;
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
		} else if (strncmp(buf, "max_energy=", 11) == 0) {
	    	double protMaxEnergy;
	    	sscanf(buf, "max_energy=%lf", &protMaxEnergy);
	    	prot->maxEnergy((bp_energy_t) protMaxEnergy);
		} else if (strncmp(buf, "threshold_energy=", 11) == 0) {
			double protThresholdEnergy;
			sscanf(buf, "threshold_energy=%lf", &protThresholdEnergy);
			prot->thresholdBindingEnergy((bp_energy_t) protThresholdEnergy);
		} else if (strncmp(buf, "palindromic=", 12) == 0) {
	    	char c;
	    	sscanf(buf, "palindromic=%c", &c);
	    	bool protPalindromic; // BindingProtein.cpp default is false
			if (c == 'y') { 	  // thus only need to worry if user set "true"
	    		protPalindromic = true;
	    	 	prot->isPalindromic(protPalindromic);
			}	
		} else if (strncmp(buf, "forward_energies:", 17) == 0) {
	    	forwRev = 1;
		} else if (strncmp(buf, "reverse_energies:", 17) == 0) {
	    	forwRev = -1;
		} else if (forwRev == 1) {
	    	prot->addForwardEnergy((bp_energy_t) strtod(buf, NULL));
		} else if (forwRev == -1) {
	    	prot->addReverseEnergy((bp_energy_t) strtod(buf, NULL));
		}

    	if(infile.eof()) {
			cout << "End of File. Add "<< protCopyNumber << " " << prot->name() << " protein(s)" <<endl; 
			if( protCopyNumber > 1) {
				for (int num = 1; num < protCopyNumber; num++) {
					std::ostringstream in;
					in << num << std::ends;
					std::string protnum = in.str();
					copyprot = new BindingProtein(*prot);
					copyprot->uniqueID(prot->name() + "-" + protnum); //e.g. SRF-1, SRF-2
					sim->addBindingProtein(copyprot);
				}
			}
			prot->uniqueID(prot->name() + "-0"); // e.g. SRF-0
			sim->addBindingProtein(prot);
		}
	} // end while
	
    return sim;
}

