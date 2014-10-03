//
// SequenceAccessibility.h
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
#ifndef _SequenceAccessibility_h_
#define _SequenceAccessibility_h_

#include <string>
#include <vector>

using namespace std;

class SequenceAccessibility {
public:
    SequenceAccessibility();
	SequenceAccessibility(const string modtype, double epigenval);
	~SequenceAccessibility();

	void addActiveEpigen( double epigenval );
	void addRepressiveEpigen( double epigenval );
	double noCoverageVal() const {return myNoCoverage;}
	double getAccessibility( const int pos, const unsigned int protLen);

private:
	vector<double> myActive;
	vector<double> myRepressed;
	int myNoCoverage; 
};

#endif	// _SequenceAccessibility_h_
