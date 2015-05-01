#ifndef dmtc_h
#define dmtc_h
using namespace std;
#include <stdlib.h>
#include <iostream>
#include <globalfunctions.h>
#include <hsp.h>
#include <ham.h>
#include <matrix.h>
class dmtc : public matrix
{
public:
    dmtc(hsp*);
    ~dmtc();
    void constructrho(int,ham*);
    void constructrho(doublereal*);
    void constructrho(cmplx*);
    void constructrho(doublereal*,doublereal*);
    double entent();
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
