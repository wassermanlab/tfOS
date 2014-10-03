//
// Sequence.h
//
// Author:
//     Rebecca Hunt 
//
// Copywrite:
//      Wasserman Lab
//      Centre For Molecular Medicine and Therapeutics
//      Child & Family Research Institute
//      University of British Columbia
//
#ifndef _Sequence_h_
#define _Sequence_h_

#include <string>
#include <map>
#include <vector>


using namespace std;

class Sequence {
public:
    Sequence(const string &seq);
    ~Sequence();
    void 			sequence(const string &seq);
    const string& 	sequence() const {return mySeq;} // if a call for seq comes in, return mySeq
    //void 			histoneMod(const vector<double> &hist);
    //bool addAccessibilityObj(SequenceAccessibility *access); // add data about active sequence regions

private:
	//static const unsigned int  		DFLT_SEQ_ACCESSIBILITY; // anything >0 is considered open to binding 
    string							mySeq;
	//SequenceAccessibility*			mySeqAccess;
    vector<BindingProtein*>				myOccupied;
	BindingProtein* isOccupiedPosition(unsigned int pos) const;
    BindingProtein* isOccupiedRegion(unsigned int start, unsigned int end) const;
    void 			setOccupiedPosition(unsigned int pos, BindingProtein* prot);
    void 			setUnoccupiedPosition(unsigned int pos);
    void 			setOccupiedRegion(unsigned int start, unsigned int end, BindingProtein* prot);
    void 			setUnoccupiedRegion(unsigned int start, unsigned int end);
	//double isSeqAccess(BindingProtein *prot, SequenceAccessibility* seqAccess) const; 
};

#endif	// _Sequence_h_
