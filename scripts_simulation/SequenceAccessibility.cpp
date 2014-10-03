//
// SequenceAccessibility.cpp
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
#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/SequenceAccessibility.h"

using namespace std;

// default constructor
SequenceAccessibility::SequenceAccessibility()
{
	myNoCoverage = -1;  //  this empty value was used to fill in the vector elements corresponding to sequence positions
					    //  not represented in the epigenModification file
						//  Also used when sequence accessibility data was not provided at all
}

// Constructor
SequenceAccessibility::SequenceAccessibility(const string modtype, double epigenval) 
{
	myNoCoverage = -1;
	if( modtype == "active"){
		addActiveEpigen( epigenval );
	}else{
		addRepressiveEpigen( epigenval );
	}
}

// Destructor
SequenceAccessibility::~SequenceAccessibility()
{
}

// set active epigenetic or open chromatin (at this stage only use one epigenone mod at a time) 
void SequenceAccessibility::addActiveEpigen( double epigenval )
{
 	myActive.push_back( epigenval );
}

// set repressive epigenetic or MNase?
void SequenceAccessibility::addRepressiveEpigen( double epigenval )
{
 	myRepressed.push_back( epigenval );
}

double SequenceAccessibility::getAccessibility( const int pos, const unsigned int protLen)  
{
	//look up epigenetic value at this location from pos to pos+protLen 
	// Do I return active or repress?
	// do I have a function to decide on accessibility for the day when use more than one?
	unsigned int plen=0;
	float sum = 0;
	if( pos > 0 && !myActive.empty() ){
		while( plen < protLen){ // data 'peaks' have values, between peaks is -1 
			sum += myActive[pos+plen]; // if on the "edge" of an accessibile region -1 will pull down the value but not abolish it
			plen++;
		}
	}
	if( pos > 0 && myActive.empty() && !myRepressed.empty() ){
		while( plen < protLen){ // data 'peaks' have values, between peaks is -1 
			sum += myRepressed[pos+plen]; // if on the "edge" of an accessibile region -1 will pull down the value but not abolish it
			plen++;
		}
		// reverse the repressed sign: if region has repressed scores then tell simulator -ve value to represent repressed,
		//   if was -ve sum, because had no repressive data, then report +ve value to represent "not repressed"
		sum *= -1;  
	}
	return sum/protLen;
}



