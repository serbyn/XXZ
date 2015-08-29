#ifndef matrix_h
#define matrix_h
#include <iostream>

using namespace std;
// dealing with LAPACK in a universal way
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
//#include <vecLib/clapack.h>
typedef __CLPK_doublereal doublereal;
typedef __CLPK_integer integer;
typedef __CLPK_doublecomplex cmplx;
#else
#include <mkl.h>
typedef double doublereal;
typedef MKL_Complex16 cmplx;
typedef double doublereal;
typedef MKL_INT integer;
// or typedef long int integer; on some systems
#endif
class matrix{
    public:
    matrix(int size,bool flag);
    ~matrix();
    void printA();
    void printM();
    void printW();
    void backup(); // copies matrix A to M, to save initial Hamiltonian!
    void diagonalize();
    void svd();
    void svd(int); //svdying submatrix -- for entanglement!!
    void svdc();
    void getv(int,double*);
    void makec();
    integer getinfo();
    void resetinfo();
    //variables
    integer size;
    integer size2;
    doublereal *A;
    doublereal *M;
    doublereal *W;
    
    cmplx* Ac;

    private:
    // variables for diagonalization
    bool isc;
    doublereal *Work;
    integer *Iwork;
    integer n, lda, nr, ldar, lwork, liwork, info;
    char jobz, uplo, jobzr;
    // separate variables for SVD routine
    doublereal *Srho;
    doublereal *Urho;
    doublereal *VTrho;
    integer LDUrho, LDVTrho, nrho, lworkrho;
    integer *Iworkrho;
    // separate variables for complex SVD
    cmplx* Urhoc;
    cmplx* VTrhoc;
    cmplx* Workrhoc;
    doublereal* RWorkc;
    

};
#endif
