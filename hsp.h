#ifndef hsp_h
#define hsp_h
using namespace std;
#include <stdlib.h>
#include <globalfunctions.h>
#include <iostream>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#include <vecLib/clapack.h>
typedef __CLPK_doublereal doublereal;
typedef __CLPK_integer integer;
typedef __CLPK_doublecomplex cmplx;
#define mkl false
#else
#include <mkl.h>
typedef MKL_INT integer;
#endif

class hsp{
public:
    //________________________________Functions_________________________________
    hsp( int systemSize );
    ~hsp();
    int array_to_int( int *);
    void int_to_array( int *, int);
    void printhspace();
    integer gethdim();
    integer gethdimh();
    //________________________________Variables_________________________________
    integer nspin;
    integer nspin2;
    // Naive Hilbert space dimension
    integer hdim0;
    // True Hilert space dimension in sector Sz=0
    integer hdim, hdim2, hdimh;
    // hspace and dictionary
    int *hspace;
    int *dict;
    int *mult; //array of multiplicities for half-chain
    int **multgen; //array of multiplicities for first spins from 1.. nspin2
private:
    // helper array
    int *a;
    void makehspace();
};
#endif

/*
 void set_conf(int **);
 int** get_conf();
 int **conf;
 */
