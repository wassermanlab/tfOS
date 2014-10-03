//
// Sequnce.cpp
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
#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/Sequnce.h"

using namespace std;

// default constructor
Sequnce::Sequnce()
{
	myLength = 0;
	myStrand = 1; // 1 is +ve strand
	myPosition = 0;
}

// Destructor
Sequnce::~Sequnce()
{
}

// set sequence occupied positions
BindingProtein* Sequence::isOccupiedPosition(unsigned int pos) const
{
    return myOccupied[pos];
}

inline void Sequence::setOccupiedPosition(unsigned int pos, BindingProtein* prot)
{
    myOccupied[pos] = prot;
}

inline void Sequence::setUnoccupiedPosition(unsigned int pos)
{
    myOccupied[pos] = NULL;
}

BindingProtein* Sequence::isOccupiedRegion( unsigned int start, unsigned int end) const
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

void Sequence::setOccupiedRegion( unsigned int start, unsigned int end, BindingProtein* prot)
{
    for (register unsigned int pos = start; pos <= end; pos++)
    	myOccupied[pos] = prot;
}

void Sequence::setUnoccupiedRegion( unsigned int start, unsigned int end)
{
    for (register unsigned int pos = start; pos <= end; pos++)
    	myOccupied[pos] = NULL;
}



