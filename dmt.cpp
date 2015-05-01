#include <iostream>
#include <stdio.h>
#ifdef __APPLE__
#include <math.h>
#else
#include "/usr/include/math.h"
#endif
#include <memory.h>
#include <dmt.h>
#include <globalfunctions.h>
//constructing a density matrix
dmt::dmt(hsp *thehsp1) : matrix(thehsp1->hdimh,false)
{
    thehsp = thehsp1;
    //cout << "Constructing Density matrix for a spin chain...\n";
    nspin = thehsp->nspin;
    nspin2 = nspin/2;
    hdim0 = thehsp->hdim0;
    hdim = thehsp->hdim;
    hdim2 = thehsp->hdim2;
    hdimh = thehsp->hdimh;
    
}

// construct density matrix rho, at the input have Hamiltionian
void dmt::constructrho(int numb, ham* H){
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
        A[i]=0.0;
    }
    // trace out second half of spins (remember that the binary representation is reversed)
    ist = 0;
    for (i=0; i<hdimh; i++) { // loop over blocks with fixed Sz for a half of spins
        for (j=0; j<thehsp->mult[i]; j++) { // loop within each block 
            for (k=0; k<thehsp->mult[i]; k++) {
                ind1 = (thehsp->hspace[ist+j])%hdimh;
                ind2 = (thehsp->hspace[ist+k])%hdimh;
                A[ind1+ind2*hdimh] = A[ind1+ind2*hdimh] + H->A[j+ist+numb*hdim]*H->A[k+ist+numb*hdim];
            }
        }
        ist = ist + thehsp->mult[i];
    }
}

// construct rho, at the input have REAL vector
void dmt::constructrho(double* vect){
    int i,j,k,ist;
    int ind1, ind2;
    // zero out matrix rho
    for(int i=0; i<hdimh*hdimh; i++){
        A[i]=0.;
    }
    // trace out second half of spins
    ist = 0;
    for (i=0; i<hdimh; i++) {
        for (j=0; j<thehsp->mult[i]; j++) {
            for (k=0; k<thehsp->mult[i]; k++) { // sum goes over all states of the leftmost spins (smallest numbers)
                ind1 = (thehsp->hspace[ist+j])%hdimh;
                ind2 = (thehsp->hspace[ist+k])%hdimh;
                A[ind1+ind2*hdimh] = A[ind1+ind2*hdimh]+ vect[j+ist]*vect[k+ist];
            }
        }
        ist = ist + thehsp->mult[i];
    }
}

// calculate entanglemenet entropy
double dmt::entent(){
    double tr,ent,dt;
    double thr0 = 1./((double)hdimh)/((double)hdimh);
    double thr1 = 1-1./((double)hdimh);
    int numz1 = 0;
    int numz2 = 0;
    ent = 0.;
    tr = 0.;
    if (W[0]>1.){
        ent=0.;
        numz1 = hdimh-1;
        numz2 = 1;
    }else{
        for (int i=0; i<hdimh; i++) {
            dt = W[i];
            tr+= dt;
            ent+= -dt*log((dt==0)?1:dt);
            if (dt<thr0) {
                numz1+=1;
            } else if (dt>thr1){
                numz2+=1;
            }
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

// versions that work for arbitrary number of spins
void dmt::constructrho(int numb, ham* H, int ns){
    int i,j,k,ist,hdimhl;
    int ind1, ind2;
    if ((ns<1)||(ns>nspin-1)){
        cout<<"dmt::improper number of spins!!!"<<endl;
        exit(-1);
    }
    hdimhl = (int)pow(2,ns);
    // zero out matrix rho
    for(int i=0; i<hdimh*hdimh; i++){
        A[i]=0.0;
    }
    // trace out second nspin-ns of spins (remember that the binary representation is reversed)
    ist = 0;
    i=-1;
    while (thehsp->multgen[ns-1][i+1]>0) {// loop over blocks with fixed Sz for a ns spins
        i++;
        for (j=0; j<thehsp->multgen[ns-1][i]; j++) { // loop within each block
            for (k=0; k<thehsp->multgen[ns-1][i]; k++) {
                ind1 = (thehsp->hspace[ist+j])%hdimhl;
                ind2 = (thehsp->hspace[ist+k])%hdimhl;
                A[ind1+ind2*hdimhl] = A[ind1+ind2*hdimhl] + H->A[j+ist+numb*hdim]*H->A[k+ist+numb*hdim];
            }
        }
        ist = ist + thehsp->multgen[ns-1][i];
    }
}

// versions that work for arbitrary number of spins
void dmt::constructrho(double *vect, int ns){
    int i,j,k,ist,hdimhl;
    int ind1, ind2;
    if ((ns<1)||(ns>nspin-1)){
        cout<<"dmt::improper number of spins!!!"<<endl;
        exit(-1);
    }
    hdimhl = (int)pow(2,ns);
     // zero out matrix rho
    for(int i=0; i<hdimh*hdimh; i++){
        A[i]=0.0;
    }
    // trace out second nspin-ns of spins (remember that the binary representation is reversed)
    ist = 0;
    i=-1;
    while (thehsp->multgen[ns-1][i+1]>0) {// loop over blocks with fixed state for large part of numbers
        i++;
        //cout<<"i="<<i<<"; mult="<<thehsp->multgen[ns-1][i]<<endl;
        for (j=0; j<thehsp->multgen[ns-1][i]; j++) { // loop within each block
            for (k=0; k<thehsp->multgen[ns-1][i]; k++) {
                ind1 = (thehsp->hspace[ist+j])%hdimhl;
                ind2 = (thehsp->hspace[ist+k])%hdimhl;
                A[ind1+ind2*hdimhl] = A[ind1+ind2*hdimhl] + vect[j+ist]*vect[k+ist];
                //cout<<ind1<<" "<<ind2<<endl;
            }
        }
        ist = ist + thehsp->multgen[ns-1][i];
    }
}

// calculate entanglemenet entropy
double dmt::entent(int ns){
    if ((ns<1)||(ns>nspin-1)){
        cout<<"dmt::improper number of spins!!!"<<endl;
        exit(-1);
    }
    int hdiml = (int)pow(2,ns);
    double tr,ent,dt;
    ent = 0.;
    tr = 0.;
    if (W[0]>1.){
        ent=0.;
    }else{
        for (int i=0; i<hdiml; i++) {
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

dmt::~dmt() {
    //cout << "Destructing density matrix... bye!\n";
}




