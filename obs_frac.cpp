//
//  obs_frac.cpp
//  diagl
//
//  Created by Maksym Serbyn on 5/2/15.
//  Copyright (c) 2015 Maksym Serbyn. All rights reserved.
//
#include "obs_frac.h"
#include <obs.h>
#include <hsp.h>
#include <ham.h>
#include <globalfunctions.h>

obs_frac::obs_frac(hsp *thehsp1, ham* theham1, Ran* theran1, int nruns1, FILE* filep1): obs(thehsp1,theham1,theran1,nruns1,filep1)
{
    dt = 2*( (int) (hdim/nspin/2));
    if (dt<2) dt=2;
    nbins = 1000;
    ez=hdim/2;
    cout<<"% in the middle we have "<< dt <<" states to study"<<endl;
    // observables
    hist_wf = createdouble(nruns,nbins);
    wf_overlap = createdouble(nruns, nbins);
    wf_overlap_typ = createdouble(nruns, nbins);
}


obs_frac::~obs_frac() {
    // QUANTITIES I
    destroy(hist_wf,nruns);
    destroy(wf_overlap,nruns);
    destroy(wf_overlap_typ,nruns);
}

void obs_frac::findez(){
    for (int i=0; i<hdim; i++) {
        if(theham->W[i]>0){
            ez=i;
            break;
        }
    }
    //cout<<"%found zero energy at="<<ez<<endl;
    if ((ez-dt/2<0)||(ez+dt/2>=hdim)){
        cout<< "obs_band::error -- too large dt, choosing ez to be hdim/2"<<endl;
        ez = hdim/2;
    }
    // searching for half-width of the band
    int indup=ez+(int)(hdim/4);
    if (indup>hdim) indup=hdim-1;
    int inddo=ez-(int)(hdim/4);
    if (indup<0) inddo=0;
}


void obs_frac::measure_hist_wf(double* hwf){
    int col,row,ind;
    integer cel;
    for (col=0; col<nbins; col++) {
        hwf[col]=0;
    }
    for (col=ez-dt/2; col<ez+dt/2; col++) {
        for (row=ez-dt/2; row<ez+dt/2; row++) {
            cel = row+((integer)hdim)*col;
            if (absv(theham->A[cel])>1e-15){
                ind = (int)floor(-log(absv(theham->A[cel]))/30*nbins);
            } else {
                ind = nbins-1;
            }
            if (ind>=nbins) ind = nbins-1;
            if (ind<0) ind = 0;
            hwf[ind]+=1.;
        }
    }
};

void obs_frac::measure_wf_overlap(double* hwf,double* hwf_log){
    int col1, col2, row, ind;
    integer cel1,cel2;
    int *cntr = new int[nbins];
    for (row=0; row<nbins; row++) {
        hwf[row] = 0.;
        hwf_log[row] = 0.;
        cntr[row] = 0;
    }    
    double ps12,ov,ov_log;
    // column specifies eigenvector number
    for (col1=ez-dt/2; col1<ez+dt/2; col1++) {
        for (col2=ez-dt/2; col2<ez+dt/2; col2++) {
            ind = (int)((theham->W[col1]-theham->W[col2]+2.5)/5.*(nbins));
            if ((ind>=0)&(ind<nbins)){
                ov =0;
                //doing loop over all coefficients
                for (row=0; row<hdim; row++) {
                    cel1 = row+hdim*col1;
                    cel2 = row+hdim*col2;
                    ps12 = pow(theham->A[cel1],2)*pow(theham->A[cel2],2);
                    ov+=ps12;
                    ov_log+=log(ps12);
                }
                ov/=(double)hdim;
                ov_log/=(double)hdim;
                hwf[ind]+=ov;
                hwf_log[ind]+=ov_log;
                cntr[ind]++;
            }
        }
    }
    for (row=0; row<nbins; row++) {
        if (cntr[row]>0){
            hwf[row] /= (double)cntr[row];
            hwf_log[row] /= (double)cntr[row];
        }
    }
    delete[] cntr;
}



void obs_frac::measuredata(){
    findez();
    measure_hist_wf(hist_wf[crun]);
    measure_wf_overlap(wf_overlap[crun],wf_overlap_typ[crun]);
    crun += 1;
}


// writing data to file=filep, which must be initialized!!!
void obs_frac::writedata() {
    double *darr3, *darr4;
    darr3 = new double[nbins];
    darr4 = new double[nbins];
    getav(hist_wf, crun, nbins, darr3, darr4);
    fput("hist for wf coefficients", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(hist_wf, crun, nbins, darr3, darr4);
    fput("mean overlap with energy", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(hist_wf, crun, nbins, darr3, darr4);
    fput("typical overlap with energy", darr3, nbins);
    fput("disp", darr4, nbins);
    
    delete[] darr3;
    delete[] darr4;
    
}


