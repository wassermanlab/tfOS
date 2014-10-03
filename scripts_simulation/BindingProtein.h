//
// BindingProtein.h
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
#ifndef _BindingProtein_h_
#define _BindingProtein_h_

#include <string>
#include <vector>
#include <utility>
using namespace std;

typedef enum {BP_FREE, BP_TETHERED, BP_BOUND} bp_state_t;
typedef double bp_energy_t;

using namespace std;

class BindingProtein {
public:
    BindingProtein();
    BindingProtein(const string &name, 
				const string &uniqueID, 
				unsigned int length, 
    			unsigned int distance,
				unsigned int maxHop,
				double hopProb,
				double energy,
		    	double thresholdBinding,
				const vector<bp_energy_t> &fep,
		    	const vector<bp_energy_t> &rep);
	~BindingProtein();
	BindingProtein(const BindingProtein &prot);
	BindingProtein& operator=(const BindingProtein &prot);

	void name(const string &name);
	const string& name() const {return myName;}
	void uniqueID(const string &uniqueID);
	const string& uniqueID() const {return myUniqueID;}
    void isPioneer( bool pioneer);
    bool isPioneer() const{ return myPioneerStatus;}
	void length(unsigned int length);
    unsigned int length() const {return myLength;}
    void position(int pos);
    int position() const {return myPosition;}
    void freeOrigin(int pos);
    int freeOrigin() const {return myFreeOrigin;}
    void direction(int dir);
    int direction() const {return myDirection;}
    void maxSlide(unsigned int distance);
	int maxSlide() const {return myMaxSlide;}
	void			maxHopDistance( unsigned int maxHopDist);
	unsigned int	maxHopDistance() const {return myMaxHop;}
	void   hopnotjumpProb( float hopProb);
    double hopnotjumpProb() const; 
	void isPalindromic(bool palindromic);
    bool isPalindromic() const {return myPalindromic;}
    void strandOrientation (int strandOrientation);
    int strandOrientation() const {return myStrandOrientation;}
    void state(bp_state_t newState);
    bp_state_t state() const {return myState;}
    bp_energy_t energyAt(int pos) const;
	void maxEnergy (bp_energy_t energy);
	void thresholdBindingEnergy (bp_energy_t thresholdBinding);
    bp_energy_t maxEnergy() const {return myMaxEnergy;}
	bp_energy_t thresholdBindingEnergy() const {return myThresholdBindingEnergy;}
    void forwardEnergyProfile(const vector<bp_energy_t> &fep);
    const vector<bp_energy_t>& forwardEnergyProfile() const {return myFEP;}
    void reverseEnergyProfile(const vector<bp_energy_t> &rep);
    const vector<bp_energy_t>& reverseEnergyProfile() const {return myREP;}
    void addForwardEnergy(bp_energy_t val);
    void addReverseEnergy(bp_energy_t val);
    unsigned int slidCount() const {return mySlidCount;}
    bool isFree() const {return (myState == BP_FREE);}
    bool isTethered() const {return (myState == BP_TETHERED);}
    bool isBound() const {return (myState == BP_BOUND);}
    bool maximumSlid() const {return (mySlidCount >= myMaxSlide);}
    void tetherAt(int pos, int dir, int strandOrientation);
    void bindAt(int pos);
    void bindAt(int pos, int strandOrientation);
    void fallOff();
    void slideRight();
    void slideLeft();
    void slideCurrent();
    double freeToTetheredProb(int pos) const;
    double freeToBoundProb(int pos) const;
    double tetheredToFreeProb(unsigned int fallOffDistr) const;
    double boundToFreeProb() const;
    double boundToTetheredProb() const;
    double tetheredToBoundProb(double, pair<pair<float,string>, pair<double, double> > ); // not const because calls bindingProbability(....)
    double remainBoundProb(  double, pair<pair<float,string>, pair<double, double> > ); // not const because calls bindingProbability(....)
	void setInteractingProt( BindingProtein* protein);
	BindingProtein* isInteractingProt() const {return myInteractingProt;}
private:
    static const unsigned int	BP_MAX_SLIDE;
	static const unsigned int  	BP_MAX_HOP;
	static const double		  	BP_HOP_PROB;
    static const bp_energy_t	BP_R;	// gas constant J/K*mol
    static const bp_energy_t	BP_k; 	// R/avagadro = k boltzmann constant J/K*molecule
    static const bp_energy_t	BP_T;	// temperature Kelvin
	bool			myPalindromic;
	bool			myPioneerStatus;
    int				myPosition;
    int				myFreeOrigin;
    int				myDirection;
	int				myStrandOrientation;
	unsigned int	myLength;
    unsigned int	mySlidCount;
	unsigned int    myMaxSlide;
	unsigned int	myMaxHop;
	double			myHopProb;
	bp_state_t			myState;
    bp_energy_t			myMaxEnergy;
	bp_energy_t			myThresholdBindingEnergy;
    string				myName;
	string				myUniqueID;
	vector<bp_energy_t>		myFEP;
    vector<bp_energy_t>		myREP;
    double bindingProbability() const;
	double bindingProbability(double seqaccess, pair<pair<float,string>, pair<double, double> > weightNrg); //changes this->myInteractionProt
	BindingProtein* myInteractingProt;
};

#endif	// _BindingProtein_h_
