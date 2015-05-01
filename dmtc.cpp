#include <iostream>
#include <stdio.h>
#ifdef __APPLE__
#include <math.h>
#else
#include "/usr/include/math.h"
#endif
#include <memory.h>
#include <dmtc.h>
#include <globalfunctions.h>
//constructing a density matrix
dmtc::dmtc(hsp *thehsp1) : matrix(thehsp1->hdimh,true)
{
    thehsp = thehsp1;
    //cout << "Constructing complex density matrix...\n";
    nspin = thehsp->nspin;
    nspin2 = nspin/2;
    hdim0 = thehsp->hdim0;
    hdim = thehsp->hdim;
    hdim2 = thehsp->hdim2;
    hdimh = thehsp->hdimh;
    
}

// construct density matrix rho, at the input have Hamiltionian
void dmtc::constructrho(int numb, ham* H){
    int i,j,k,ist;
    int ind1, ind2;
    
    /*// print corresponding eigenvector
    cout << "vector["<<numb<<"]= " << endl;
    for (i=0;i<hdim;i++){
        cout<< A[i+numb*hdim]<<"; ";
    }
    cout <<"\n";
    */
    // zero out matrix rho
    for(int i=0; i<hdimh*hdimh; i++){
#if mkl
        Ac[i].real=0.0;
        Ac[i].imag=0.0;
#else
        Ac[i].r=0.0;
        Ac[i].i=0.0;
#endif
    }
    // trace out first half of spins
    ist = 0;
    for (i=0; i<hdimh; i++) {
        for (j=0; j<thehsp->mult[i]; j++) {
            for (k=0; k<thehsp->mult[i]; k++) {
                ind1 = (thehsp->hspace[ist+j])%hdimh;
                ind2 = (thehsp->hspace[ist+k])%hdimh;
#if mkl
          Ac[ind1+ind2*hdimh].real = Ac[ind1+ind2*hdimh].real+ H->A[j+ist+numb*hdim]*H->A[k+ist+numb*hdim];
#else
                Ac[ind1+ind2*hdimh].r = Ac[ind1+ind2*hdimh].r+ H->A[j+ist+numb*hdim]*H->A[k+ist+numb*hdim];
#endif
            }
        }
        ist = ist + thehsp->mult[i];
    }
}

// construct rho, at the input have REAL vector
void dmtc::constructrho(doublereal* vect){
    int i,j,k,ist;
    int ind1, ind2;
    // zero out matrix rho
    for(int i=0; i<hdimh*hdimh; i++){
#if mkl
        Ac[i].real=0.0;
        Ac[i].imag=0.0;
#else
        Ac[i].r=0.0;
        Ac[i].i=0.0;
#endif
    }
    // trace out first half of spins
    ist = 0;
    for (i=0; i<hdimh; i++) {
        for (j=0; j<thehsp->mult[i]; j++) {
            for (k=0; k<thehsp->mult[i]; k++) {
                ind1 = (thehsp->hspace[ist+j])%hdimh;
                ind2 = (thehsp->hspace[ist+k])%hdimh;
#if mkl
          Ac[ind1+ind2*hdimh].real = Ac[ind1+ind2*hdimh].real+ vect[j+ist]*vect[k+ist];
#else
          Ac[ind1+ind2*hdimh].r = Ac[ind1+ind2*hdimh].r+ vect[j+ist]*vect[k+ist];
#endif
            }
        }
        ist = ist + thehsp->mult[i];
    }
}

// construct rho, at the input have COMPLEX vector
void dmtc::constructrho(cmplx* vect){
    int i,j,k,ist;
    int ind1, ind2;
    // zero out matrix rho
    for(int i=0; i<hdimh*hdimh; i++){
#if mkl
        Ac[i].real=0.0;
        Ac[i].imag=0.0;
#else
        Ac[i].r=0.0;
        Ac[i].i=0.0;
#endif
    }
    // trace out first half of spins
    ist = 0;
    for (i=0; i<hdimh; i++) {
        for (j=0; j<thehsp->mult[i]; j++) {
            for (k=0; k<thehsp->mult[i]; k++) {
                ind1 = (thehsp->hspace[ist+j])%hdimh;
                ind2 = (thehsp->hspace[ist+k])%hdimh;
#if mkl
                Ac[ind1+ind2*hdimh].real = Ac[ind1+ind2*hdimh].real + vect[j+ist].real*vect[k+ist].real + vect[j+ist].imag*vect[k+ist].imag;
                Ac[ind1+ind2*hdimh].imag = Ac[ind1+ind2*hdimh].imag + vect[j+ist].real*vect[k+ist].imag - vect[j+ist].imag*vect[k+ist].real;
                
#else
                Ac[ind1+ind2*hdimh].r = Ac[ind1+ind2*hdimh].r + vect[j+ist].r*vect[k+ist].r + vect[j+ist].i*vect[k+ist].i;
                Ac[ind1+ind2*hdimh].i = Ac[ind1+ind2*hdimh].i + vect[j+ist].r*vect[k+ist].i - vect[j+ist].i*vect[k+ist].r;
                
#endif
            }
        }
        ist = ist + thehsp->mult[i];
    }
}

// construct rho, at the input have COMPLEX vector
void dmtc::constructrho(doublereal* vectr,doublereal* vecti){
    int i,j,k,ist;
    int ind1, ind2;
    // zero out matrix rho
    for(int i=0; i<hdimh*hdimh; i++){
#if mkl
        Ac[i].real=0.0;
        Ac[i].imag=0.0;
#else
        Ac[i].r=0.0;
        Ac[i].i=0.0;
#endif 
    }
    // trace out first half of spins
    ist = 0;
    for (i=0; i<hdimh; i++) {
        for (j=0; j<thehsp->mult[i]; j++) {
            for (k=0; k<thehsp->mult[i]; k++) {
                ind1 = (thehsp->hspace[ist+j])%hdimh;
                ind2 = (thehsp->hspace[ist+k])%hdimh;
#if mkl
                Ac[ind1+ind2*hdimh].real = Ac[ind1+ind2*hdimh].real + vectr[j+ist]*vectr[k+ist] + vecti[j+ist]*vecti[k+ist];
                Ac[ind1+ind2*hdimh].imag = Ac[ind1+ind2*hdimh].imag + vectr[j+ist]*vecti[k+ist] - vecti[j+ist]*vectr[k+ist];
                
#else
                Ac[ind1+ind2*hdimh].r = Ac[ind1+ind2*hdimh].r + vectr[j+ist]*vectr[k+ist] + vecti[j+ist]*vecti[k+ist];
                Ac[ind1+ind2*hdimh].i = Ac[ind1+ind2*hdimh].i + vectr[j+ist]*vecti[k+ist] - vecti[j+ist]*vectr[k+ist];
#endif
            }
        }
        ist = ist + thehsp->mult[i];
    }
}


// calculate entanglemenet entropy after SVD'ing
double dmtc::entent(){
    double tr,ent,dt;
    ent = 0.;
    tr = 0.;
    if (W[0]>1.){
        ent=0.;
    }else{
        for (int i=0; i<hdimh; i++) {
            dt = W[i];
            tr+= dt;
            ent+= -dt*log((dt==0)?1:dt);
        }
    }
    entropy = ent;
    //cout<<"\n Trace of rho= "<<tr<<endl;
    //cout<<"\n Entropy= "<<ent<<endl;
    if ((ent<-0.01)||(ent!=ent)){
        cout<<"error!!! entropy is negative="<<ent;
        exit(1);
    }
    return ent;
}


dmtc::~dmtc() {
    //cout << "Destructing complex density matrix... bye!\n";
}




