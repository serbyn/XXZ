#ifndef obs_band_h
#define obs_band_h
#include <hsp.h>
#include <ham.h>
#include <dmt.h>
#include <obs.h>
#include <iostream>
#include <sstream>
using namespace std;
// class includes all observables that depend only on a SINGLE HAMILTONIAN and do not include any time evolution
// quantities included are:
// entanglement entropy in the band
// statistics of spectrum
// ...

class obs_band: public obs{
    public:
    obs_band(hsp *thehsp1, ham* theham1, Ran* theran1, int nruns1, FILE* filep1, dmt* thedmt);
    ~obs_band();
    void measuredata(); // measures data for each instance of Hamiltonian
    void writedata(); // average data and writes data at the end of the run
    
    private:
    dmt* thedmt; // density matrix needed for some routines
    int band; // number of portions of the band
    int band2; // number of portions in the energy band
    int dt; // number of states in a portion of the band
    int ez; // zero energy
    int hw; // half-width
    void findez();// searching for zero energy

    // stroring the data I
    double **entropypos; // ent entropy
    double **mmav; // magnetization
    double **lspos; //level statistics
    int **lshist; // levels histogram
    double **iprepos; //energy IPR
    int **hystdmte;
    int **hystspace; // levels spacing of entanglement
    double **snapent;
    double **hystSz; // histogram of Sz
    double **sentl;
    
    //generating the data I
    void ententpos(); // entanglement entropy at several points
    void ententposE(); // entanglement entropy at several points and Energy fluctuations
    void szavav(); // average magnetization
    void szav(double*, int); // helper function calculating average for a given state
    void levelstat();// level spacing statistics
    double levelstatband(int si, int ei);// level spacing statistics helper
    double levelstatband(int si, int ei,int *hist);// level spacing statistics helper with histogram
    
    // stroring the data for histogram II-1
    int nbins; // number of bins
    int nobs; // number of operators
    double *A2; // matrix after the action of a operator
    double *A2c; // matrix with coefficients of decomposition, the main object
    double **hyst1; // array with hystogram data
    double **hyst2; // array with hystogram data
    double **mdE; // array with hystogram data
    int ** cE;
    double **mdE2; // array with hystogram data
    int ** cE2;
    double **ae0; // array for distribution of matrix elements with energy
    double **ae1; // array for distribution of matrix elements with energy
    double **aeoE; // array for distribution of matrix elements with energy
    double **ae0L; // array for distribution of matrix elements with energy
    double **aeoEL; // array for distribution of matrix elements with energy
    double **aew; // array for distribution of matrix elements with energy
    double **snapv1; // one snapshot vector per run to be saved
    double **mdd; // array with distribution of m/d in the middle
    double **mdc;
    double **mdce;
    // generating the data II-1
      //spin action

    void szA2(); // constructing the matrix after action of sz
    void szszA2(); // constructing the matrix after action of szsz
    void spsmA2(); // constructing the matrix after action of S+S- + S-S+
    void ssA2(); // constructing the matrix after action of 1/2+SS
    void sp__smA2(); // constructing the matrix after action of S+S- + S-S+ at the opposite ends
    void sp__smMA2();// constructing the matrix after action of S+S- - S-S+ at the opposite ends
    void szmA2(int ns); // constructing the matrix after action of sz on the generic spin
    void spsmmA2(int ns); // constructing the matrix after action of S+S- + S-S+ on generic spin ns and ns+1
    void szA2vec(int ns); // constructing the matrix after action of sz using optimized code
    void spsmmA2vec(int ns); // constructing the matrix after action of S+S- + S-S+ on generic spin ns and ns+1
    void sp__smA2vec(); // constructing the matrix after action of S+S- + S-S+ at the opposite ends
    void currA2();// action of current operator
     //measuring functions
    void buildA2c(); // fills in A2c -- matrix of coefficients
    void fillhist2sort(double *h1,double *h2,double *md, int *ce, double *md2, int *ce2,double *mddistr, double *mdc, double *mdde); // fill histogram of matrix elements and adjacent energy levels
    void fillhist2(double *h1,double *h2,double *md, int *ce, double *md2, int *ce2,double *mddistr, double *mdc, double *mdde); // fill histogram of matrix elements and adjacent energy levels
    void fillhist2fast(double *h2,double *md, int *ce, double *md2, int *ce2,double *mddistr); // fill histogram of matrix elements and adjacent energy levels
    void fillae(double *ae,double *ae1,double *ae2,double *,double *,double *); // fills ae
    
    // measuring statistics of matrix elements
    void measureAstat2(); // measures data set II-1 and II-2 simultaneously
};
#endif
