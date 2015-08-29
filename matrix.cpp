#include <iostream>
#include <cstdlib>
#include <memory.h>
#include <globalfunctions.h>
#include <matrix.h>
//#include <complex>
//#define mkl true
//typedef complex<double> dcmplx;
// dealing with LAPACK in a universal way
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
//#include <vecLib/clapack.h>
typedef __CLPK_doublereal doublereal;
typedef __CLPK_integer integer;
typedef __CLPK_doublecomplex cmplx;
#include <math.h>
#define mkl false
#else
typedef MKL_INT integer;
#include <mkl.h>
#include <stdlib.h>
#include "/usr/include/math.h"
#define mkl true
#endif
matrix::matrix(int size1, bool flag)
{
    if (flag == false){
        isc = false;
        size = size1;
        size2 = (integer)size*(integer)size;
        A = new doublereal[size2];
        M = new doublereal[1];
        W = new doublereal[size];
        jobz = 'V';
        uplo = 'U';
        n = size;
        lda = size;
        //workspace query for dsyevd_ for a given size
        lwork = -1;
        Work = new doublereal[1];
        Iwork = new integer[1];
        dsyevd_(&jobz, &uplo, &n, A, &lda, W, Work, &lwork, Iwork, &liwork, &info);
        //cout <<"%lwork= "<< Work[0]<<endl;
        lwork = (integer)(Work[0]);
        liwork =  Iwork[0];
        //cout <<"%lwork= "<< Work[0]<<" iwork= "<< Iwork[0]<<endl;
        delete [] Work;
        delete [] Iwork;
        Work = new doublereal[lwork];
        Iwork = new integer[liwork];
        //workspace query for dgesdd_ for a given size
        jobzr = 'N';
        nr = size;
        ldar = size;
        Iworkrho = new integer[8*size];
        doublereal * Workrho;
        Workrho = new doublereal[1];
        LDUrho = size;
        LDVTrho = size;
        nrho = size;
        lworkrho = -1;
        dgesdd_(&jobzr, &nrho, &nrho, A, &nrho, W, Urho, &LDUrho, VTrho, &LDVTrho, Workrho, &lworkrho, Iworkrho, &info);
        lworkrho = (integer)(Workrho[0]);
        //cout <<"lworkrho= "<< Workrho[0]<<endl;
        delete [] Workrho;
        //Workrho = new doublereal[lworkrho];
    }else{
        isc = true;
        size = size1;
        size2 = (integer)size*(integer)size;
        Ac = new cmplx[size2];
        W = new doublereal[size];
        //workspace query for cgesdd_ for a given size
        jobzr = 'N';
        nr = size;
        ldar = size;
        Workrhoc = new cmplx[1];
        RWorkc = new doublereal[7*size];
        Iworkrho = new integer[8*size];
        LDUrho = size;
        LDVTrho = size;
        nrho = size;
        lworkrho = -1;
        zgesdd_(&jobzr, &nrho, &nrho, Ac, &nrho, W, Urhoc, &LDUrho, VTrhoc, &LDVTrho, Workrhoc, &lworkrho, RWorkc, Iworkrho, &info);
#if mkl
        lworkrho =  (integer)(Workrhoc[0].real);
#else
        lworkrho =  (integer)(Workrhoc[0].r);
#endif
        //cout <<"lworkrho= "<< lworkrho<<endl;
        delete [] Workrhoc;
        Workrhoc = new cmplx[lworkrho];
    }
}

matrix::~matrix() {
    if (isc==false){
        //cout << "Destructing matrix... bye!\n";
        delete[] A;
        delete[] M;
        delete[] Iworkrho;
        delete[] W;
        delete[] Work;
        delete[] Iwork;
    }else{
        //cout << "Destructing complex matrix... bye!\n";
        delete[] Ac;
        delete[] W;
        delete[] Workrhoc;
        delete[] RWorkc;
        delete[] Iworkrho;
    }
}



// Print matrix A
void matrix::printA(){
    if (isc == false){
        cout << "A=" << endl;
        for(int row=0; row<size; row++) {
            for(int col=0; col<size; col++) {
                cout << A[row+col*size] << " ";
            }
            cout << endl;
        }
        cout << endl;
    } else{
        cout << "Ac=" << endl;
        for(int row=0; row<size; row++) {
            for(int col=0; col<size; col++) {
#if mkl
          cout << Ac[row+col*size].real << "+I*"<<Ac[row+col*size].imag<<"; ";
#else
          cout << Ac[row+col*size].r << "+I*"<<Ac[row+col*size].i<<"; ";
#endif
            }
            cout << endl;
        }
        cout << endl;
    }
}

// Print matrix M
void matrix::printM(){
    if (isc == false){
        cout << "M=" << endl;
        for(int row=0; row<size; row++) {
            for(int col=0; col<size; col++) {
                cout << M[row+col*size] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}


// Print eigenvalues W
void matrix::printW(){
    cout << "W=" << endl;
    for(int row=0; row<size; row++) {
        cout << W[row] << " ";
    }
    cout << endl;
}



// returning given eigenvector (after diagonalization)
void matrix::getv(int vecnum, double* vec){
    for (int i=0; i<size; i++) {
        vec[i] = (double) A[i+vecnum*size];
    }
}


//backup matrix A to M
void matrix::backup(){
    delete[] M;
    M = new double[size2];
    //copy matrix A to M
    memcpy(M,A,size2*sizeof(doublereal));
}

// Diagonalization of A
void matrix::diagonalize(){
     dsyevd_(&jobz, &uplo, &n, A, &lda, W, Work, &lwork, Iwork, &liwork, &info);
    /*// output x
     cout << "Few eigenvalues(not more than 5)=" << endl;
     for(int i=0; i<min(num_rows,50); i++) {
     cout << W[i] << "; ";
     }
     */
    //cout << "GS energy="<<W[0]<< endl;
    //cout <<"Work[0]=" <<Work[0] << endl;
    if (info<0){
        cout << "matrix::diagonalize()::dsyevd_::error with arguments, exiting.."<< endl;
        exit(1);
    }else if (info>0){
        cout << "matrix::diagonalize()::dsyevd_::failed to converge, info="<<info<< endl;
    }
}


// SVD of A
void matrix::svd(){
    //copy matrix A to M
    //memcpy(M,A,size2*sizeof(doublereal));
    doublereal* Workrho;
    Workrho = new doublereal[lworkrho];
    dgesdd_(&jobzr, &nrho, &nrho, A, &nrho, W, Urho, &LDUrho, VTrho, &LDVTrho, Workrho, &lworkrho, Iworkrho, &info);
    delete[] Workrho;
    if (info<0){
        cout << "matrix::svd()::dgesdd_::error with arguments, exiting.."<< endl;
        exit(1);
    }else if (info>0){
        cout << "matrix::svd()::dgesdd_::failed to converge, info="<<info<< endl;
    }
}

// SVD of part of A for ns spins
void matrix::svd(int ns){
    integer hdiml=(integer)pow(2,ns);
    //copy matrix A to M
    //memcpy(M,A,size2*sizeof(doublereal));
    doublereal* Workrho;
    Workrho = new doublereal[lworkrho];
    dgesdd_(&jobzr, &hdiml, &hdiml, A, &hdiml, W, Urho, &LDUrho, VTrho, &LDVTrho, Workrho, &lworkrho, Iworkrho, &info);
    delete[] Workrho;
    if (info<0){
        cout << "matrix::svd()::dgesdd_::error with arguments, exiting.."<< endl;
        exit(1);
    }else if (info>0){
        cout << "matrix::svd()::dgesdd_::failed to converge, info="<<info<< endl;
    }
}

// complex SVD of A
void matrix::svdc(){
    if (isc == false){
        cout<<"error";
        //exit(1);
    }else{
    zgesdd_(&jobzr, &nrho, &nrho, Ac, &nrho, W, Urhoc, &LDUrho, VTrhoc, &LDVTrho, Workrhoc, &lworkrho, RWorkc, Iworkrho, &info);
    }
    if (info<0){
        cout << "matrix::svdc()::zgesdd_::error with arguments, exiting.."<< endl;
        exit(1);
    }else if (info>0){
        cout << "matrix::svdc()::zgesdd_::failed to converge, info="<<info<< endl;
    }
}

integer matrix::getinfo(){
    return info;
}

void matrix::resetinfo(){
    info = 0;
}






