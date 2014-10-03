//
// ProteinBindingSimulator.h
//
// Author:
//      David Arenillas (dave@cmmt.ubc.ca)
//
// Copywrite:
//      Wasserman Lab
//      Centre For Molecular Medicine and Therapeutics
//      Child & Family Research Institute
//      University of British Columbia
//
#ifndef _ProteinBindingSimulator_h_
#define _ProteinBindingSimulator_h_

#include <string>
#include <map>
#include <vector>

#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/BindingProtein.h"
#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/ProteinInteractions.h"
#include "/raid6/rebecca/project_SIMULATION/occupancySim_C++/scripts_simulation/SequenceAccessibility.h"

using namespace std;

class ProteinBindingSimulator {
public:
    ProteinBindingSimulator(const string &seq);
    ~ProteinBindingSimulator();
    void 			sequence(const string &seq);
    const string& 	sequence() const {return mySeq;} // if a call for seq comes in, return mySeq
    void 			histoneMod(const vector<double> &hist);
    void			fallOffDistribution(unsigned int fallOffDistr);
	unsigned int	fallOffDistribution() const {return myFallOffDistr;} 
	void            maxAllCoopDistances(unsigned int coopDist);
	unsigned int    maxAllCoopDistances() const {return myMaxCoopDists;} // give back the distance, if asked
	const vector<unsigned int>& proteinTetheredCounts(const string &name)
								{return myTetheredCounts[name];}
    const vector<unsigned int>& proteinBoundCounts(const string &name)
								{return myBoundCounts[name];}
	bool addInteractionsObj(ProteinInteractions *interac); // add interactions obj so sim can use them
    bool addAccessibilityObj(SequenceAccessibility *access); // add data about active sequence regions
	bool addBindingProtein( BindingProtein *prot); // there is private vector of BindingProteins, so sim can use them 
    const vector< BindingProtein* >& bindingProteins() const{return myProteins;}
    bool simulate(unsigned int iters); // the code that runs the interations

private:
	static const unsigned int		DFLT_FALLOFF_DISTR;
	static const unsigned int  		DFLT_SEQ_ACCESSIBILITY; // anything >0 is considered open to binding 
	static const pair<pair<float, string>, pair<double,double> > DFLT_NOINTERACTION; // <<weight,bindingstate>, <nrg, nrgThr>> 
	int								myIters;
    string							mySeq;
	unsigned int					myFallOffDistr; //default value
	unsigned int                    myMaxCoopDists; // default value
	ProteinInteractions*			myInteractions;
	SequenceAccessibility*			mySeqAccess;
	vector<BindingProtein*>				myProteins;
    vector<BindingProtein*>				myOccupied;
    map<string, vector<unsigned int> >	myTetheredCounts;
    map<string, vector<unsigned int> >	myBoundCounts;
    bool 			testProb(float prob) const;
	BindingProtein* isOccupiedPosition(unsigned int pos) const;
    BindingProtein* isOccupiedRegion(unsigned int start, unsigned int end) const;
    void 			setOccupiedPosition(unsigned int pos, BindingProtein* prot);
    void 			setUnoccupiedPosition(unsigned int pos);
    void 			setOccupiedRegion(unsigned int start, unsigned int end, BindingProtein* prot);
    void 			setUnoccupiedRegion(unsigned int start, unsigned int end);
    int				translocate( BindingProtein* prot, unsigned int seqLen);
    void incTetheredCounts(
						const string &protName, unsigned int start, unsigned int end);
    void incBoundCounts(
						const string &protName, unsigned int start, unsigned int end);
	double isSeqAccess(BindingProtein *prot, SequenceAccessibility* seqAccess) const; 
	pair<pair<float,string>, pair<double, double> > isInteraction(BindingProtein *prot1, ProteinInteractions *inter ) const;
};

#endif	// _ProteinBindingSimulator_h_
