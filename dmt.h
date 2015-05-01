#ifndef dmt_h
#define dmt_h
using namespace std;
#include <stdlib.h>
#include <iostream>
#include <globalfunctions.h>
#include <hsp.h>
#include <ham.h>
#include <matrix.h>
class dmt : public matrix
{
public:
    dmt(hsp*);
    ~dmt();
    void constructrho(int,ham*);
    void constructrho(double*);
    double entent();
    // versions that construct DMT for part of the system
    void constructrho(int,ham*,int);
    void constructrho(double*,int);
    double entent(int);
    double entropy;
    // calculating entanglement entropy at several points
    //void ententpos();

private:
    hsp* thehsp;
    int nspin;
    int nspin2;
    // Naive Hilbert space dimension
    int hdim0;
    // True Hilert space dimension in sector Sz=0
    int hdim, hdim2, hdimh;
};
#endif
