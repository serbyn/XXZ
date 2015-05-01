#ifndef obs_h
#define obs_h
#include <hsp.h>
#include <ham.h>
#include <iostream>
#include <sstream>

// class to be included in ANY OBSERVABLE
// contains:
// references to the Hilbert space, Hamiltonian, Random generator, number of runs, and pointer to the file with data;
// functions to collect, average and write the data;
using namespace std;

class obs{
    
    public:
    obs(hsp*,ham*,Ran*, int nruns, FILE* filep);
    virtual ~obs();
    virtual void measuredata() = 0; // measures data for each instance of Hamiltonian
    //virtual void averagedata() = 0; // averages data, preparing for dumping all data
    virtual void writedata() = 0; // writes data at the end of the run
    void writebasicdata(); // writes basic data regarding the hamiltonian
    FILE* filep;
    
    protected:
    // variables
    hsp* thehsp;
    ham* theham;
    Ran* theran;
    int nspin;
    int nspin2;
    integer hdim0; // Naive Hilbert space dimension
    integer hdim, hdim2, hdimh; // True Hilert space dimension in sector Sz=0
    int nruns; // number of instances for averaging
    int crun; // current run
    // functions for writing data
    std::ostringstream sline1, sline2;
    void fput(string comment, double* array, int arrl);
    void fput(string comment, int n);
    void fput(string comment, double n);
    void stringput(string, double);
    void stringput(string, int);
    void stringwrite();
    // functions for debugging
    void printvec(double*,int);
    void printmat(double*,int,int);
    void printvec(cmplx*, int);
    
};
#endif
