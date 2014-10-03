//
// ProteinInteractions.h
//
// Author:
//	Rebecca Worsley Hunt
//
// Copywrite:
//	Wasserman Lab
//	Centre For Molecular Medicine and Therapeutics
//	Child & Family Research Institute
//	University of British Columbia
//
#ifndef _ProteinInteractions_h_
#define _ProteinInteractions_h_

#include <string>
#include <utility>
#include <map>

using namespace std;

typedef std::pair< std::pair<float, int>, string> MyPairs; //to take interaction info: state, weight, distance
typedef	std::map<string, MyPairs> MyInnerMap;
typedef std::map<string, MyInnerMap > MyMap; 

class ProteinInteractions {
public:
    ProteinInteractions();
    ProteinInteractions(const string prot1, 
					const string Prot2, 
					float weight, 
    				int distance,
					string state);
	~ProteinInteractions();

	void setInteractionVal(const string prot1, const string interactingProtein, float weight, int distance, string state);
    MyPairs getInteractionValues(const string prot1, const string interactingProt); 
	MyPairs noInteractionVal() const {return myNoInteraction;}
private:
	MyMap interacMap;
	MyMap::iterator iterMymap;
	MyInnerMap innermap;
	MyInnerMap::iterator iterInnermap;
	MyPairs myNoInteraction;
};

#endif	// _ProteinInteractions_h_
