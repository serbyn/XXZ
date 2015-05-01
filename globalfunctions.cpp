#include "globalfunctions.h"
#ifdef __APPLE__
#include <math.h>
#else
#include "/usr/include/math.h"
#endif



// Factorial of n
double fact(int n) {
    int i;
    double f = 1;
    for(i=1; i<=n; ++i)
        f *= i;
    return f;
}

// Binomial coefficients
int C(int n, int k) {
    return ( (int)(fact(n)/(fact(k)*fact(n-k))) );
}

// absolut value

double absv(double x){
    if (x<0) return -1.*x;
    return x;
}


// averages one-dimensional array of double with dispersion in a naive way
// helper function: used by getav
void getav0(double * arr, int size, double *av, double *disp){
    int i;
    double av0, disp0;
    av0 = 0.;
    disp0 = 0.;
    for (i=0; i<size; i++) {
        av0 += arr[i];
    }
    av0 = av0/((double)size);
    for (i=0; i<size; i++) {
        disp0 += (arr[i]-av0)*(arr[i]-av0);
    }
    disp0 = disp0/((double)size);
    disp0 = sqrt(disp0);
    
    av[0] = av0;
    disp[0] = disp0;
}

// averages one-dimensional array of double with dispersion
void getav(double * arr, int size, double *av, double *disp){
    int N = 5;
    int M, M0;
    if (size>N){
        M = size/N;
        M0 = M+(size%N);
    } else {
        M=1;
        M0=1;
        N=1;
    }
    int i,k;
    double * av0;
    av0 = new double[N];
    for (k=0; k<N; k++){
        av0[k] = 0.;
    }
    for (i=0; i<size; i++) {
        av0[(i/M)%N] += arr[i];
    }
    av0[0] = av0[0]/((double)M0);
    for (k=1; k<N; k++){
        av0[k] = av0[k]/((double)M);
    }
    getav0(av0, N, av, disp);
    delete[] av0;
}

// averages one-dimensional array of int with dispersion
void getav(int * arr, int size, double *av, double *disp){
    int N = 5;
    int M, M0;
    if (size>N){
        M = size/N;
        M0 = M+(size%N);
    } else {
        M=1;
        M0=1;
        N=1;
    }
    int i,k;
    double * av0;
    av0 = new double[N];
    for (k=0; k<N; k++){
        av0[k] = 0.;
    }
    for (i=0; i<size; i++) {
        av0[(i/M)%N] += arr[i];
    }
    av0[0] = av0[0]/((double)M0);
    for (k=1; k<N; k++){
        av0[k] = av0[k]/((double)M);
    }
    getav0(av0, N, av, disp);
    delete[] av0;
}

// averages one-dimensional array of double WITHOUT dispersion
void getav(double * arr, int size, double *av){
    int i;
    double av0;
    av0 = 0.;
    for (i=0; i<size; i++) {
        av0 += arr[i];
    }
    av0 = av0/((double)size);
    av[0] = av0;
}


// averages 2dim array of doubles without dispersion
void getav(double ** arr, int size1, int size2, double *av){
    int i,j;
    double av0;
    for (j=0; j<size2; j++) {
        av0 = 0.;
        for (i=0; i<size1; i++) {
            av0 += arr[i][j];
        }
        av0 = av0/((double)size1);
        av[j] = av0;
    }
}

// averages 2dim array of doubles over first dimension, results are written to arrays av and disp
void getav(double ** arr, int size1, int size2, double *av, double *disp){
    int N = 5;
    int M, M0;
    if (size1>N){
        M = size1/N;
        M0 = M+(size1%N);
    } else {
        M=1;
        M0=1;
        N=1;
    }
    int i,j,k;
    double * av0;
    av0 = new double[N];
    for (j=0; j<size2; j++) {
        for (k=0; k<N; k++){
            av0[k] = 0.;
        }
        for (i=0; i<size1; i++) {
            av0[(i/M)%N] += arr[i][j];
        }
        av0[0] = av0[0]/((double)M0);
        for (k=1; k<N; k++){
            av0[k] = av0[k]/((double)M);
        }
        getav0(av0, N, av+j, disp+j);
    }
    delete[] av0;
}

// averages 2dim array of int over first dimension, results are written to arrays av and disp
void getav(int ** arr, int size1, int size2, double *av, double *disp){
    int N = 5;
    int M, M0;
    if (size1>N){
        M = size1/N;
        M0 = M+(size1%N);
    } else {
        M=1;
        M0=1;
        N=1;
    }    int i,j,k;
    double * av0;
    av0 = new double[N];
    for (j=0; j<size2; j++) {
        for (k=0; k<N; k++){
            av0[k] = 0.;
        }
        for (i=0; i<size1; i++) {
            av0[(i/M)%N] += arr[i][j];
        }
        av0[0] = av0[0]/((double)M0);
        for (k=1; k<N; k++){
            av0[k] = av0[k]/((double)M);
        }
        getav0(av0, N, av+j, disp+j);
    }
    delete[] av0;
}



//Allocates and returns a square 2d array of integers
int** createint(int n) {
  int** ret = new (int(*[n]));
  for(int i=0;i<n;i++) ret[i] = new int[n];
  return ret;
}
int** createint(int n, int m) {
  int** ret = new (int(*[n]));
  for(int i=0;i<n;i++) ret[i] = new int[m];
  return ret;
}

//Allocates and returns a square 2d array of real numbers
double** createdouble(int n) {
  double** ret = new (double(*[n]));
  for(int i=0;i<n;i++) ret[i] = new double[n];
  return ret;
}

double** createdouble(int n, int m) {
    double** ret = new (double(*[n]));
    for(int i=0;i<n;i++) ret[i] = new double[m];
    return ret;
}

//Deallocates the 2d array ptr
void destroy( int** ptr, int n) {
  for(int i=0;i<n;i++) delete[] ptr[i];
  delete[] ptr;
}
void destroy( double** ptr, int n) {
  for(int i=0;i<n;i++) delete[] ptr[i];
  delete[] ptr;
}

void printmat(double* M, int size){
        cout << "M=" << endl;
        for(int row=0; row<size; row++) {
            for(int col=0; col<size; col++) {
                cout << M[row+col*size] << " ";
            }
            cout << endl;
        }
        cout << endl;
}


