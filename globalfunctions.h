#ifndef globalfunctions_h
#define globalfunctions_h

#include <iostream>

using namespace std;

// factorial and binomial
int    C(int, int);
double fact(int);

// absv

double absv(double);

// functions for average and dispersion of a given quantity

// one dim array
void getav(double * arr, int size, double *av, double *disp);
void getav0(double * arr, int size, double *av, double *disp);
void getav(int * arr, int size, double *av, double *disp);

void getav(double * arr, int size, double *av); // without dispersion
// two dim array, average goes over first index!
void getav(double ** arr, int size1, int size2, double *av, double *disp);
void getav(int ** arr, int size1, int size2, double *av, double *disp);

void getav(double ** arr, int size1, int size2, double *av); // without dispersion


// functions for conversion of state into number


//2D array creation functions
int **createint(int);
int **createint(int,int);

double **createdouble(int);
double **createdouble(int n, int m);

void destroy(int**, int);
void destroy(double**, int);

/*//matrix print functions
void write_m(double**, int);
void write_m(complex<double>**, int);
void write_m(int**, int);
*/

// Print matrix M
void printmat(double* M, int);

#endif

