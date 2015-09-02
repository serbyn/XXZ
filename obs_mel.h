#ifndef obs_mel_h
#define obs_mel_h
#include <hsp.h>
#include <ham.h>
#include <dmt.h>
#include <obs.h>
#include <iostream>
#include <sstream>
using namespace std;
// class includes matrix elements observables; does not include any time evolution
// quantities included are:

class obs_mel: public obs{
    public:
    obs_mel(hsp *thehsp1, ham* theham1, Ran* theran1, int nruns1, FILE* filep1);
    ~obs_mel();
    void measuredata(); // measures data for each instance of Hamiltonian
    void writedata(); // average data and writes data at the end of the run
    void szA2vec(int ns); // constructing the matrix after action of sz using optimized code
    void buildA2c(); // fills in A2c -- matrix of coefficients
    double *A2c; // matrix with coefficients of decomposition, the main object
    private:
    int band2; // number of portions in the energy band
    int dt; // number of states in a portion of the band
    int ez; // zero energy
    int hw; // half-width
    void findez();// searching for zero energy
    
    // stroring the data for histograms 
    int nbins; // number of bins
    int nobs; // number of operators
    int iter;// storing current operator
    double *A2; // matrix after the action of a operator
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
    double **mdce2;
    double **beta,**betaL,**m,**mL,**moES,**moESL;
    double *iprM, **iprMstat;
    // different operators
    void szA2(); // constructing the matrix after action of sz
    void szszA2(); // constructing the matrix after action of szsz
    void spsmA2(); // constructing the matrix after action of S+S- + S-S+
    void ssA2(); // constructing the matrix after action of 1/2+SS
    void sp__smA2(); // constructing the matrix after action of S+S- + S-S+ at the opposite ends
    void sp__smMA2();// constructing the matrix after action of S+S- - S-S+ at the opposite ends
    void szmA2(int ns); // constructing the matrix after action of sz on the generic spin
    void spsmmA2(int ns); // constructing the matrix after action of S+S- + S-S+ on generic spin ns and ns+1
    void spsmmA2vec(int ns); // constructing the matrix after action of S+S- + S-S+ on generic spin ns and ns+1
    void sp__smA2vec(); // constructing the matrix after action of S+S- + S-S+ at the opposite ends
    void currA2();// action of current operator
    void curr1A2(int spinpos); // current operator on given 2 sites
    //measuring functions
    void fillhist2sort(double *h1,double *h2,double *md, int *ce, double *md2, int *ce2,double *mddistr, double *mdc, double *mdde, double *mdde2); // fill histogram of matrix elements and adjacent energy levels
    void fillhist2(double *h1,double *h2,double *md, int *ce, double *md2, int *ce2,double *mddistr, double *mdc, double *mdde); // fill histogram of matrix elements and adjacent energy levels
    void fillhist2fast(double *h2,double *md, int *ce, double *md2, int *ce2,double *mddistr); // fill histogram of matrix elements and adjacent energy levels
    void fillae(double *ae,double *ae1,double *ae2,double *,double *,double *); // fills ae
    void fillaN(double *beta,double *betaL,double *m,double *mL,double *moeS,double *moeSL); // fills the same as ae but with index
    void fillIPR(double *iprM,double *iprMstat);//filling IPR's
    // measuring statistics of matrix elements
    void measureAstat2(); // measures data set II-1 and II-2 simultaneously
};
#endif
