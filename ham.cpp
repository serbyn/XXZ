#include <iostream>
#include <time.h>
#include <stdio.h>
#ifdef __APPLE__
#include <math.h>
#else
#include "/usr/include/math.h"
#endif
#include "hsp.h"
#include "ham.h"
#include <fstream>
#include <sstream>
#include <memory.h>
#include <globalfunctions.h>

//constructing a spin-half chain of a given size
ham::ham(hsp *thehsp1) : matrix(thehsp1->hdim,false)
{
    thehsp = thehsp1;
    cout << "%Constructing Hamiltonian for a spin chain...\n";
    nspin = thehsp->nspin;
    nspin2 = nspin/2;
    hdim0 = thehsp->hdim0;
    hdim = thehsp->hdim;
    hdim2 = thehsp->hdim2;
    hdimh = thehsp->hdimh;
    
// Creating parameters of Hamiltonian
    hz = new double[nspin];
    Jz = new double[nspin];
    Jp = new double[nspin];
    
// Variables for creation of H
    ain = new int[nspin];
    aout = new int[nspin];
    numout = new int[nspin+1];
    coefout = new double[nspin+1];
    
// storing M and U
    Ml = new double[hdimh*hdimh];
    Ul = new double[hdimh];
}


ham::~ham() {
    cout << "%Destructing Hamiltonian... bye!\n";
    delete[] hz;
    delete[] Jz;
    delete[] Jp;
    delete[] ain;
    delete[] aout;
    delete[] numout;
    delete[] coefout;
    delete[] Ml;
    delete[] Ul;
}


// Set parameters
void ham::setparam(double * h0, double * Jz0, double * Jp0){
    memcpy(hz,h0,sizeof(double)*nspin);
    memcpy(Jz,Jz0,sizeof(double)*(nspin));
    memcpy(Jp,Jp0,sizeof(double)*(nspin));
}

// Get parameters
double ham::gethz(){
    return hz[0];
}
double ham::getJz(){
    return Jz[0];
}
// Diagonal part of H, returns energy
double ham::Hdiag(double * hz, double * Jz,int * ain){
    double e=hz[nspin-1]*(ain[nspin-1]-0.5)+Jz[nspin-1]*(ain[nspin-1]-0.5)*(ain[0]-0.5);
    for(int i = 0;i<nspin-1;i++) {
        e+=hz[i]*(ain[i]-0.5)+Jz[i]*(ain[i]-0.5)*(ain[i+1]-0.5);
    }
    return  e;
}

// Off-diagonal part of H, returns vector of states and coefficients
void ham::Hoffdiag(double *Jp, int *ain, int *aout, int *numout, double *coefout){
    int ind=0;
    for(int i = 0;i<nspin-1;i++) {
        if (ain[i]+ain[i+1] == 1){
            memcpy(aout,ain,sizeof(int)*nspin);
            aout[i+1] = ain[i];
            aout[i]=ain[i+1];
            numout[ind]=thehsp->array_to_int(aout);
            coefout[ind]=0.5*Jp[i];
            ind= ind+1;
        }
    }
    if (ain[nspin-1]+ain[0] == 1){
        memcpy(aout,ain,sizeof(int)*nspin);
        aout[0] = ain[nspin-1];
        aout[nspin-1]=ain[0];
        numout[ind]=thehsp->array_to_int(aout);
        coefout[ind]=0.5*Jp[nspin-1];
        ind= ind+1;
    }
    numout[ind] = -1;
}

// Creates Hamiltonian
void ham::createH(){
    //zero out matrix A
    for(integer i=0; i<hdim2; i++){
        A[i]=0.0;
    }
    //fill A
    integer col, row, i, psin, ind;
    for(col=0; col<hdim; col++) {
        psin = thehsp->hspace[col];
        thehsp->int_to_array(ain,psin);
        // filling diagonal element
        ind = col+col*hdim;
        if (ind<0) {cout<< "ham::createH():: integer overflow exiting"; exit(-1);}
        A[col+col*hdim]=Hdiag(hz,Jz,ain);
        // finding off-diagonal elements
        Hoffdiag(Jp,ain,aout,numout,coefout);
        // filling matrix entities
        i = 0;
        while (numout[i]>-1) {
            row = (integer)thehsp->dict[numout[i]];
            if (row<col+1) {
                ind = row+col*hdim;
                if (ind<0) {cout<< "ham::createH():: integer overflow exiting"; exit(-1);}
                A[ind]=coefout[i];
            }
            i = i+1;
        }
    }

}


// Diagonal part of H2, returns energy, open BC
double ham::Hdiag2(double * hz, double * Jz,int * ain){
    double e=hz[nspin2-1]*(ain[nspin2-1]-0.5);
    for(int i = 0;i<nspin2-1;i++) {
        e+=hz[i]*(ain[i]-0.5)+Jz[i]*(ain[i]-0.5)*(ain[i+1]-0.5);
    }
    return  e;
}

// Off-diagonal part of H2, returns vector of states and coefficients
void ham::Hoffdiag2(double *Jp, int *ain, int *aout, int *numout, double *coefout){
    int ind=0;
    for(int i = 0;i<nspin2-1;i++) {
        if (ain[i]+ain[i+1] == 1){
            memcpy(aout,ain,sizeof(int)*nspin);
            aout[i+1] = ain[i];
            aout[i]=ain[i+1];
            numout[ind]=thehsp->array_to_int(aout);
            coefout[ind]=0.5*Jp[i];
            ind= ind+1;
        }
    }
    numout[ind] = -1;
}



// creates Hamiltonian for the firs half of the spins...
// diagonalizes it and copies it into M and U matrices respectively...
integer ham::createH2(){
    matrix mtmp(hdimh,false);
    // zero out matrix A
    for(int i=0; i<hdimh*hdimh; i++){
        mtmp.A[i]=0.0;
    }
    // fill A
    int col, row, i, psin;
    for(col=0; col<hdimh; col++) {
        psin = col;
        thehsp->int_to_array(ain,psin);
        // filling diagonal element
        mtmp.A[col+col*hdimh]=Hdiag2(hz,Jz,ain);
        // finding off-diagonal elements
        Hoffdiag2(Jp,ain,aout,numout,coefout);
        // filling matrix entities
        i = 0;
        while (numout[i]>-1) {
            row = numout[i];
            if (row<col+1) {
                mtmp.A[row+col*hdimh]=coefout[i];
            }
            i = i+1;
        }
    }
    // diagonalize
    mtmp.diagonalize();
    // copy
    memcpy(Ml,mtmp.A,hdimh*hdimh*sizeof(double));
    memcpy(Ul,mtmp.W,hdimh*sizeof(double));
    return mtmp.getinfo();
}

/*void ham::diagonalizeH(){
    diagonalize();
    integer info1 = getinfo();
    if (info1>0) {
        cout<<"positive info, trying to rediagonalize!!!"<<endl;
        createH();
        diagonalize();
        info1 = getinfo();
        cout<<"after rediagonalization info is"<<getinfo()<<endl;
    }
}
*/



