#ifndef ham_h
#define ham_h
using namespace std;
// local libraries
#include <hsp.h>
#include <matrix.h>
#include <globalfunctions.h>
// global libraries
#include <stdlib.h>
#include <iostream>
class ham: public matrix // inherits matrix class!
{
public:
    //________________________________Functions_________________________________
    ham(hsp*);
    ~ham();
    void setparam(double *, double *, double *);
    void setrandhz(double lowlim, double uplim); 
    void setrandJz(double lowlim, double uplim); 
    void setrandJp(double lowlim, double uplim);
    double gethz();
    double getJz();
    integer createH2();// creates Hamiltonian for the second half of spins
    double *Ml, *Ul;
    void createH();
    //void diagonalizeH(); // version of diagonalization which handles exceptions
    //double szav(int pos, int numb);
    //void szav(double *sz, int numb);
    //________________________________Variables_________________________________

/*    
    void set_conf(int **); 
    int** get_conf();
	int **conf;
*/
private:
    hsp* thehsp;
    int nspin;
    int nspin2;
    integer hdim0;
    integer hdim, hdim2, hdimh;
     // parameters of Hamiltonian
    double * hz;    
    double * Jz;
    double * Jp;
    // constructing Hamiltonian
    double Hdiag(double * hz, double * Jz,int * ain);
    void Hoffdiag(double *Jp, int *ain, int *aout, int *numout, double *coefout);
    // constructing Hamiltonian for the half-chain
    // constructing Hamiltonian
    double Hdiag2(double * hz, double * Jz,int * ain);
    void Hoffdiag2(double *Jp, int *ain, int *aout, int *numout, double *coefout);
    int *ain;
    int *aout;
    int *numout;
    double *coefout;
    double hzmax;
};
#endif
