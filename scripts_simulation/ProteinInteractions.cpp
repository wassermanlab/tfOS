//
// ProteinInteractions.cpp
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
#include <cstdlib>
#include <iostream>
#include <utility>
#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/ProteinInteractions.h"

using namespace std;

// default constructor
ProteinInteractions::ProteinInteractions()
{
	myNoInteraction = make_pair(make_pair(0,0), ""); //empty weight, distance, state
}

// Constructor
ProteinInteractions::ProteinInteractions(string prot1, string prot2, float weight, int distance, string state) 
{
	setInteractionVal(prot1, prot2, weight, distance, state);	
	
	myNoInteraction = make_pair(make_pair(0,0),""); // weight=0 is the key
}

// Destructor
ProteinInteractions::~ProteinInteractions()
{
}

// set protein interactions
void ProteinInteractions::setInteractionVal( string prot1, string prot2, float weight, int distance, string state)
{
	// insert in map both prot1->prot2 and prot2->prot1 so that iterative lookups don't have to test both
	// forward
	innermap.insert( map<string, MyPairs>::value_type( prot2, make_pair(make_pair(weight, distance),state) ) );
	interacMap.insert( MyMap::value_type( prot1, innermap ) );
	innermap.clear();
	// reverse
	innermap.insert( map<string, MyPairs>::value_type( prot1, make_pair(make_pair(weight, distance),state) ) );
	interacMap.insert( MyMap::value_type( prot2, innermap ) );
	innermap.clear();
	return;
}

// get interaction
MyPairs ProteinInteractions::getInteractionValues(const string prot1, const string interactingProt)  
{
	iterMymap = interacMap.find(prot1);
	if( iterMymap != interacMap.end()){
		// first protein found
		iterInnermap = (iterMymap->second).find(interactingProt);
		if( iterInnermap != (iterMymap->second).end()){
			// interacting protein found
			return iterInnermap->second; // return the weight, distance, state as MyPairs
		}else{
			return myNoInteraction; //not found. So return a useless set of values
		}
	}else{
		return myNoInteraction; //not found
	}
	return myNoInteraction; //not found
}



