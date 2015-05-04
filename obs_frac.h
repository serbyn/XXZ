#ifndef diagl_obs_frac_h
#define diagl_obs_frac_h

#include <hsp.h>
#include <ham.h>
#include <ran.h>
#include <obs.h>
#include <iostream>
#include <sstream>
using namespace std;
// class includes all observables that depend only on a SINGLE HAMILTONIAN and do not include any time evolution
// quantities included are:
// entanglement entropy in the band
// statistics of spectrum
// ...

class obs_frac: public obs{
public:
    obs_frac(hsp *thehsp1, ham* theham1,Ran* theran1, int nruns1, FILE* filep1);
    ~obs_frac();
    void measuredata(); // measures data for each instance of Hamiltonian
    void writedata(); // average data and writes data at the end of the run
    
private:
    int nbins; // number of bins
    int dt; // number of states in a portion of the band
    int ez; // zero energy
    void findez();// searching for zero energy
    
    // stroring the data
    double **hist_wf; // wf coefficients histogram
    double **wf_overlap; // wf overlaps with energy
    double **wf_overlap_typ; // typical wf overlaps with energy
    
    // generating the data
    void measure_hist_wf(double* );
    void measure_wf_overlap(double*,double*);
};

#endif
