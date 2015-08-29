#include <iostream>
#include <time.h>
#include <stdio.h>
#include <ran.h>
#include <obs.h>
#include <obs_band.h>
#include <hsp.h>
#include <ham.h>
#include <dmt.h>
#include <dmtc.h>
#include <fstream>
#include <sstream>
#include <memory.h>
#include <globalfunctions.h>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
//#include <vecLib/cblas.h>
typedef __CLPK_doublereal doublereal;
typedef __CLPK_integer integer;
typedef __CLPK_doublecomplex cmplx;
#include <math.h>
#else
#include <mkl.h>
#include "/usr/include/math.h"
#endif
// Preamble to use sorting in the code:
#include <algorithm>
#include <vector>
#include <map>
#define STAT true // if true, calculating all statistics, if false -- only between NN

obs_band::obs_band(hsp *thehsp1, ham* theham1, Ran* theran1, int nruns1, FILE* filep1, dmt* thedmt1): obs(thehsp1,theham1,theran1,nruns1,filep1)
{
    thedmt = thedmt1;
    band = 1;
    band2 = 101;
    dt = 2*( (int) (hdim/nspin/2));
    if (dt<2) dt=2;
    nbins = 1000;
    nobs = 4;
    ez=hdim/2;
    hw = 1;
    cout<<"% in the middle we have "<< dt <<" states to study"<<endl;
    // QUANTITIES I
    entropypos = createdouble(nruns,band);
    hystdmte = createint(nruns*band, nbins); // stores entanglement histogram
    hystspace = createint(nruns*band, nbins); // stores entanglement levels spacing histogram
    hystSz = createdouble(nruns*band2,17);
    snapent = createdouble(nobs,hdim);
    lspos = createdouble(nruns,2);
    lshist = createint(nruns,nbins);
    mmav = createdouble(nruns,band);
    sentl = createdouble(nruns,nspin2);
    // QUANTITIES II-1
    A2 = new double[hdim2];
    if (STAT){
        A2c = new double[hdim2];
    }else{
        A2c = new double[hdim];
    }
    hyst1 = createdouble(nobs*nruns,nbins); // histogram of matrix elements
    hyst2 = createdouble(nobs*nruns,nbins); // histogram of matrix elements between nn
    mdE = createdouble(nobs*nruns,band2); // average mdE as a function of relative energy
    mdE2 = createdouble(nobs*nruns,band2); // average mdE as a function of absolute energy
    mdd = createdouble(nobs*nruns,nbins);
    mdc = createdouble(nobs*nruns,nbins);
    mdce = createdouble(nobs*nruns,band2);
    cE = createint(nruns,band2); // average mdE as a function of relative energy
    cE2 = createint(nruns,band2); // average mdE as a function of absolute energ
    ae0 = createdouble(nruns*nobs,nbins);
    ae1 = createdouble(nruns*nobs,nbins);
    aeoE = createdouble(nruns*nobs,nbins);
    ae0L = createdouble(nruns*nobs,nbins);
    aeoEL = createdouble(nruns*nobs,nbins);
    aew = createdouble(nruns*nobs,nbins);
    snapv1 = createdouble(nobs,hdim);
}


obs_band::~obs_band() {
    // QUANTITIES I
    destroy(entropypos,nruns);
    destroy(mmav,nruns);
    destroy(lspos,nruns);
    destroy(lshist,nruns);
    destroy(hystdmte,nruns*band);
    destroy(hystspace,nruns*band);
    destroy(hystSz,nruns*band2);
    destroy(snapent,nobs);
    destroy(sentl,nruns);

    // QUANTITIES II-1
    delete[] A2;
    delete[] A2c;
    destroy(hyst1,nobs*nruns);
    destroy(hyst2,nobs*nruns);
    destroy(mdE,nobs*nruns);
    destroy(cE,nruns);
    destroy(mdE2,nobs*nruns);
    destroy(cE2,nruns);
    destroy(ae0,nobs*nruns);
    destroy(ae1,nobs*nruns);
    destroy(aeoE,nobs*nruns);
    destroy(ae0L,nobs*nruns);
    destroy(aeoEL,nobs*nruns);
    destroy(aew,nobs*nruns);
    destroy(snapv1,nobs);
    destroy(mdd,nobs*nruns);
    destroy(mdc,nobs*nruns);
    destroy(mdce,nobs*nruns);
}
// QUANTITIES II-1
// creating matrix A2 = 2sz A, where sz acts on the largest index or rightmost spin
void obs_band::szA2(){
    integer col,row;
    int t = (int)pow(2,nspin-1);
    memcpy(A2, theham->A, sizeof(double)*hdim2); //start with the matrix A
    for (row=0; row<hdim; row++) {
        if ((thehsp->hspace[row]/t)%2==0){
            for (col=0; col<hdim; col++){
                A2[row+hdim*col] = -A2[row+hdim*col];
            }
        }
    }
}
// sz in the generic spin ns, where ns is counted from the right of the chain/from the largest index
// ns runs from 1 to nspin
void obs_band::szA2vec(int ns){
    integer row;
    int t = (int)pow(2,nspin-ns);
    memcpy(A2, theham->A, sizeof(double)*hdim2); //start with the matrix A
    for (row=0; row<hdim; row++) {
        if ((thehsp->hspace[row]/t)%2==0){
            cblas_dscal(hdim, -1., A2+row, hdim);
        }
    }
}

void obs_band::szmA2(int ns){
    integer col,row;
    int t = (int)pow(2,nspin-ns);
    memcpy(A2, theham->A, sizeof(double)*hdim2); //start with the matrix A
    for (row=0; row<hdim; row++) {
        if ((thehsp->hspace[row]/t)%2==0){
            for (col=0; col<hdim; col++){
                A2[row+hdim*col] = -A2[row+hdim*col];
            }
        }
    }
}

// creating matrix A2 = 2sz 2sz A
void obs_band::szszA2(){
    integer col,row,tmp;
    int t = (int)pow(2,nspin-2);
    memcpy(A2, theham->A, sizeof(double)*hdim2); //start with the matrix A
    for (row=0; row<hdim; row++) {
        tmp = thehsp->hspace[row]/t;
        if (((tmp%2==0)&((tmp/2)%2==1))||((tmp%2==1)&((tmp/2)%2==0))){
            for (col=0; col<hdim; col++){
                A2[row+hdim*col] = -A2[row+hdim*col];
            }
        }
    }
}

// creating matrix A2 = (S+S-+S-S+) A
void obs_band::spsmA2(){
    integer col,row,tmp,s1,s2,b,row2;
    int t = (int)pow(2,nspin-2);
    //memcpy(A2, theham->A, sizeof(double)*hdim2); //start with the matrix A
    for(col=0; col<hdim2; col++) {
        A2[col]=0.;
    }
    for (row=0; row<hdim; row++) {
        tmp = thehsp->hspace[row]/t;
        b = thehsp->hspace[row]%t;
        s1 = tmp%2;
        s2 = (tmp/2)%2;
        if (s1!=s2){
            row2 = thehsp->dict[b+t*s2+2*t*s1];
            //if ((row2>hdim-1)||(row2<0)) cout<<"error with row2!!"<<endl;
            //cout<<"tmp="<<tmp<<", row = "<<2*(2*(b/2)+s1)+s2<<endl;
            for (col=0; col<hdim; col++){
                A2[row+hdim*col] =  theham->A[row2+hdim*col];
            }
        }
    }
    /*
    cout << "A2=" << endl;
    int size = hdim;
    for(row=0; row<size; row++) {
        for(col=0; col<size; col++) {
            cout << A2[row+col*size] << " ";
        }
        cout << endl;
    }
    cout << endl;
*/
}

// creating matrix A2 = (S+S-+S-S+) A in the spins ns and ns+1
// ns should run from 1 to nspin (the measurment across the cut is implemented below)
void obs_band::spsmmA2(int ns){
    integer col,row,tmp,s1,s2,row2;
    int t = (int)pow(2,ns-1);
    if (ns==nspin) sp__smA2();
        else{
            //memcpy(A2, theham->A, sizeof(double)*hdim2); //start with the matrix A
            for(col=0; col<hdim2; col++) {
                A2[col]=0.;
            }
            for (row=0; row<hdim; row++) {
                tmp = thehsp->hspace[row]/t;
                s1 = tmp%2;
                s2 = (tmp/2)%2;
                if (s1!=s2){
                    // cout<<"i1="<<thehsp->hspace[row]<<"; i2="<<thehsp->hspace[row] -s1*t+s2*t -s2*2*t+s1*2*t<<endl;
                    row2 = thehsp->dict[thehsp->hspace[row] -s1*t+s2*t -s2*2*t+s1*2*t];
                    //if ((row2>hdim-1)||(row2<0)) cout<<"error with row2!!"<<endl;
                    //cout<<"tmp="<<tmp<<", row = "<<2*(2*(b/2)+s1)+s2<<endl;
                    for (col=0; col<hdim; col++){
                        A2[row+hdim*col] =  theham->A[row2+hdim*col];
                    }
                }
            }
        }
}

// creating matrix A2 = (S+S-+S-S+) A where S+ and S- are at the opposite ends
void obs_band::sp__smA2(){
    integer col,row,tmp,s1,s2,b,row2;
    int t = (int)pow(2,nspin-1);
    //memcpy(A2, theham->A, sizeof(double)*hdim2); //start with the matrix A
    for(col=0; col<hdim2; col++) {
        A2[col]=0.;
    }
    for (row=0; row<hdim; row++) {
        tmp = thehsp->hspace[row]/t;
        b = 2*((int)(thehsp->hspace[row]%t)/2);
        s1 = tmp%2;// spin at one end
        s2 = thehsp->hspace[row]%2; // spin at the other end
        if (s1!=s2){
            row2 = thehsp->dict[b+s1+t*s2];
            //if ((row2>hdim-1)||(row2<0)) cout<<"error with row2!!"<<endl;
            //cout<<"tmp="<<tmp<<", row = "<<2*(2*(b/2)+s1)+s2<<endl;
            for (col=0; col<hdim; col++){
                A2[row+hdim*col] =  theham->A[row2+hdim*col];
            }
        }
    }
}


// creating matrix A2 = (S+S- - S-S+) A where S+ and S- are at the opposite ends
void obs_band::sp__smMA2(){
    integer col,row,tmp,s1,s2,b,row2;
    int t = (int)pow(2,nspin-1);
    //memcpy(A2, theham->A, sizeof(double)*hdim2); //start with the matrix A
    for(col=0; col<hdim2; col++) {
        A2[col]=0.;
    }
    for (row=0; row<hdim; row++) {
        tmp = thehsp->hspace[row]/t;
        b = 2*((int)(thehsp->hspace[row]%t)/2);
        s1 = tmp%2;// spin at one end
        s2 = thehsp->hspace[row]%2; // spin at the other end
        if ((s1!=s2)){
            row2 = thehsp->dict[b+s1+t*s2];
            if (s1==0){
                //if ((row2>hdim-1)||(row2<0)) cout<<"error with row2!!"<<endl;
                //cout<<"tmp="<<tmp<<", row = "<<2*(2*(b/2)+s1)+s2<<endl;
                for (col=0; col<hdim; col++){
                    A2[row+hdim*col] =  theham->A[row2+hdim*col];
                }
            } else {
                for (col=0; col<hdim; col++){
                    A2[row+hdim*col] =  -theham->A[row2+hdim*col];
                }
            }
                
        }
    }
}

// creating matrix A2 = j A where j is current operator
void obs_band::currA2(){
    integer col,row,tmp,s1,s2,b,row2;
    int t;
    //memcpy(A2, theham->A, sizeof(double)*hdim2); //start with the matrix A
    for(col=0; col<hdim2; col++) {
        A2[col]=0.;
    }
    // loop over spin position before end
    for (int spinpos = 0; spinpos<nspin-1; spinpos++) {
        t = (int)pow(2,spinpos);// spinpos ranges from 0 to nspin-1, hence t goes from 1 to 2^(nspin-1)
        for (row=0; row<hdim; row++) {
            tmp = thehsp->hspace[row]/t;
            s1 = tmp%2;
            s2 = (tmp/2)%2;
            b = 2*((int)(thehsp->hspace[row]%t)/2);
            if ((s1!=s2)){
                row2 = thehsp->dict[thehsp->hspace[row] -s1*t+s2*t -s2*2*t+s1*2*t];
                if (s1==0){
                    cblas_daxpy(hdim, 1., theham->A+row2, hdim, A2+row, hdim);
                } else {
                    cblas_daxpy(hdim, -1., theham->A+row2, hdim, A2+row, hdim);
                }
            }
        }
    }
    // repeat across the cut
    t = (int)pow(2,nspin-1);// spinpos ranges from 0 to nspin-1, hence t goes from 1 to 2^(nspin-1)
    for (row=0; row<hdim; row++) {
        tmp = thehsp->hspace[row]/t;
        s1 = tmp%2;
        s2 = thehsp->hspace[row]%2; // spin at the other end
        b = 2*((int)(thehsp->hspace[row]%t)/2);
        if ((s1!=s2)){
            row2 = thehsp->dict[b+s1+t*s2];
            if (s1==0){
                cblas_daxpy(hdim, 1., theham->A+row2, hdim, A2+row, hdim);
            } else {
                cblas_daxpy(hdim, -1., theham->A+row2, hdim, A2+row, hdim);
            }
        }
    }
}


// creating matrix A2 = (S+S-+S-S+) A where S+ and S- are at the opposite ends
void obs_band::sp__smA2vec(){
    integer col,row,tmp,s1,s2,b,row2;
    int t = (int)pow(2,nspin-1);
    for(col=0; col<hdim2; col++) {
        A2[col]=0.;
    }
    for (row=0; row<hdim; row++) {
        tmp = thehsp->hspace[row]/t;
        b = 2*((int)(thehsp->hspace[row]%t)/2);
        s1 = tmp%2;// spin at one end
        s2 = thehsp->hspace[row]%2; // spin at the other end
        if (s1!=s2){
            row2 = thehsp->dict[b+s1+t*s2];
            cblas_dcopy(hdim, theham->A+row2, hdim, A2+row, hdim);
        }
    }
}

// creating matrix A2 = (S+S-+S-S+) A in the spins ns and ns+1
// ns should run from 1 to nspin (the measurment across the cut is implemented below)
void obs_band::spsmmA2vec(int ns){
    integer col,row,tmp,s1,s2,row2;
    int t = (int)pow(2,ns-1);
    if (ns==nspin) sp__smA2vec();
    else{
        for(col=0; col<hdim2; col++) {
            A2[col]=0.;
        }
        for (row=0; row<hdim; row++) {
            tmp = thehsp->hspace[row]/t;
            s1 = tmp%2;
            s2 = (tmp/2)%2;
            if (s1!=s2){
                row2 = thehsp->dict[thehsp->hspace[row] -s1*t+s2*t -s2*2*t+s1*2*t];
                cblas_dcopy(hdim, theham->A+row2, hdim, A2+row, hdim);
            }
        }
    }
}

// creating matrix A2 = (S.S) A
void obs_band::ssA2(){
    integer col,row,tmp,s1,s2,b,row2;
    int t = (int)pow(2,nspin-2);
    //memcpy(A2, theham->A, sizeof(double)*hdim2); //start with the matrix A
    for(col=0; col<hdim2; col++) {
        A2[col]=0.;
    }
    for (row=0; row<hdim; row++) {
        tmp = thehsp->hspace[row]/t;
        b = thehsp->hspace[row]%t;
        s1 = tmp%2;
        s2 = (tmp/2)%2;
        if (s1==s2){
            for (col=0; col<hdim; col++){
                A2[row+hdim*col] = theham->A[row+hdim*col];
            }
        }else{
            row2 = thehsp->dict[b+t*s2+2*t*s1];
            //if ((row2>hdim-1)||(row2<0)) cout<<"error with row2!!"<<endl;
            //cout<<"tmp="<<tmp<<", row = "<<2*(2*(b/2)+s1)+s2<<endl;
            for (col=0; col<hdim; col++){
                A2[row+hdim*col] =  theham->A[row2+hdim*col];
            }
        }
    }
    /*
     cout << "A2=" << endl;
     int size = hdim;
     for(row=0; row<size; row++) {
     for(col=0; col<size; col++) {
     cout << A2[row+col*size] << " ";
     }
     cout << endl;
     }
     cout << endl;
     */
}


// builds A2c
void obs_band::buildA2c(){
    int m,n,k,lda,ldb,ldc;
    double alpha,beta;
    if (STAT){
        // first, do A2=A^T.A2
        m = hdim;
        n = hdim;
        k = hdim;
        alpha = 1.;
        beta = 0.;
        lda = m;
        ldb = k;
        ldc = m;
        cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,alpha,theham->A,lda,A2,ldb,beta,A2c,ldc);//using dgemm to multiply: Matrix A^T: hdim*hdim times A2 and write to A2c
        /*for (int j=0; j<hdim-1; j++) {
            cout<<"v1="<<cblas_ddot(hdim, theham->A+(j+1)*hdim, 1,A2+j*hdim, 1)<<"; v2="<<A2c[j+1+hdim*j]<<endl;
        }*/
    } else{ // multiplying just vectors to get nearest neighbor matrix elements
        for (int j=0; j<hdim-1; j++) {
            A2c[j]=cblas_ddot(hdim, theham->A+(j+1)*hdim, 1,A2+j*hdim, 1);
        }
    }
    /*cout << "A2c=" << endl;
    int size = hdim;
    double t;
    for(int row=0; row<size; row++) {
        for(int col=0; col<size; col++) {
            t=A2c[row+col*size];
            if (absv(t)>0.0001) cout << t << " ";
            else cout<<"0 ";
        }
        cout << endl;
    }
    cout << endl;
     */
}

template <class RAIter, class Compare>
void argsort(RAIter iterBegin, RAIter iterEnd, Compare comp,
             std::vector<size_t>& indexes) {
    
    std::vector< std::pair<size_t,RAIter> > pv ;
    pv.reserve(iterEnd - iterBegin) ;
    
    RAIter iter ;
    size_t k ;
    for (iter = iterBegin, k = 0 ; iter != iterEnd ; iter++, k++) {
        pv.push_back( std::pair<int,RAIter>(k,iter) ) ;
    }
    
    std::sort(pv.begin(), pv.end(),
              [&comp](const std::pair<size_t,RAIter>& a, const std::pair<size_t,RAIter>& b) -> bool
              { return comp(*a.second, *b.second) ; }) ;
    
    indexes.resize(pv.size()) ;
    std::transform(pv.begin(), pv.end(), indexes.begin(),
                   [](const std::pair<size_t,RAIter>& a) -> size_t { return a.first ; }) ;
}
// procedure that works with usual spectrum: computes:
// distribution of Mnn 
// average of Mnn*deltaE with deltaE
// average of Mnn^2/deltaE with deltaE
// curvature




// procedure that works with reshuffled spectrum: computes
// (1) nn mel distribution: ==== h1 *
// (2) average m/delta:     ==== md *
// (3) mel over delta E     ==== h2 *
// (4) average m/delta over delta E = md2 *
// (3) m/delta distribution in the middle of the band ==== mddistr
// (4) fraction of resonances with N in the middle of the band ==== resf
// (5) average number of resonances with energy density ==== resav

// fills in the histogram of matrix elements and adjacent energy levels
void obs_band::fillhist2sort(double *h1, double *h2, double *md, int *counts2, double *md2, int *counts3,double *mddistr,double *resf,double *resav){
    integer row, col, ind,ind2, inde1,inde2, cel;
    double inc  = 1./(hdim);
    double Z=0.25; // factor over which we multiply matrix element
    double DE = 0.25; // energy range over which we average
    int *counts4 =  new int[band2];
    int *countsR =  new int[nbins];
    // reset h1
    for (col=0; col<nbins; col++){
        h1[col]=0.;
        h2[col]=0.;
        mddistr[col]=0.;
        resf[col] = 0.;
        countsR[col]=0;
    }
    for (col=0; col<band2; col++){
        counts2[col]=0;
        counts3[col]=0;
        counts4[col]=0;
        md[col]=0;
        md2[col]=0;
        resav[col]=0.;
    }
    // getting reshuffled energy spectrum
    // we define new array, coinciding with the old spectrum+correction
    double *W =  new double[hdim];
    for (row= 0; row<hdim; row++) {
        W[row] = theham->W[row]+Z*A2c[row+((integer)hdim)*row];
             //cout<<W[row]<<",";
    }
    // sorting it: now we have the correct order in indices, array W is left intact!!! To go in the correct order after reshuffling, we now have to use indexes[] as an index!
    std::vector<size_t> indexes;
    argsort(W, W + hdim, std::less<double>(), indexes);
    
    double dloclog =0.; //local typical level spacing
    double dloc =0.; //local level spacing
    double mtyp=0;
    int ldos=0;
    int nres;
    integer irow,icol;
    for (row= 0; row<hdim; row++) { // loop over individual states within the band
        irow = (integer)indexes[row];
        mtyp = 0.;
        dloclog=0.;
        dloc=0.;
        nres = 0;
        ldos = 0;
        // calculating ldos
        for (col=0; col<hdim; col++) { // loop over vectors for these states; we exclude symmetric copy
            icol = (integer)indexes[col];
            if((absv(W[irow]-W[icol])<DE)&((irow!=icol))) ldos++;
        }
        // index for the energy histogram
        inde1 = (int) (((theham->W[irow])/((double)nspin)/2.+0.5)*(band2-1));
        if (inde1>=band2-1) inde1=band2-2;
        if (inde1<0) inde1=0;
        for (col=0; col<hdim; col++) { // loop over vectors for these states; we exclude symmetric copy
            icol = (integer)indexes[col];
            // doing something only if we are within energy range
            if((absv(W[irow]-W[icol])<DE)&((irow!=icol))){
                cel = irow+((integer)hdim)*icol;
                if (absv(A2c[cel])>1e-15){
                    ind = (int)floor(-log(absv(A2c[cel]))/30*nbins);
                } else {
                    ind = nbins-1;
                }
                if (ind>=nbins) ind = nbins-1;
                if (ind<0) ind=0;
                // filling distribution over NN
                if (irow==icol+1)  h1[ind] += inc;
                // filling distribution over energy range:
                h2[ind] += 1./((double) ldos);
                // typical matrix element
                mtyp +=log(absv(A2c[cel]));
                // ENERGY resolved probes
                // md: over NN
                if (irow==icol+1){
                    md[inde1]+=log(absv(A2c[cel]/(W[irow]-W[icol])));
                    counts2[inde1]++;
                    // filling histogram if we are in the middle
                    if ((irow>=ez-dt/2)&(irow<ez+dt/2)){
                        inde2=(log(absv(A2c[cel]/(W[irow]-W[icol])))+30)/60*nbins;
                        if (inde2>=nbins) inde2=nbins-1;
                        if (inde2<0) inde2=0;
                        mddistr[inde2]+=1.;
                    }
                }
                // md2: over energy range: calculating level spacing AND matrix element
                if (icol>0){
                    dloclog += log(W[icol]-W[icol-1]);
                    dloc += W[icol]-W[icol-1];
                }
                mtyp+=log(absv(A2c[cel]));
            }
            // calculating resonances (energy resolved)
            if (absv(A2c[cel])>absv(W[icol]-W[irow])){
                nres++;
            }
            // probability of resonance for the middle -- over all levels
            if((row>=ez-dt/2)&(row<ez+dt/2)){
                ind2  = abs(irow-icol);
                if (ind2>nbins-1) ind2=nbins-1;
                countsR[ind2]+=1;
                if (absv(A2c[cel])>absv(W[icol]-W[irow])){
                    resf[ind2]+=1;
                }
                //if (ind==1){
                //    cout<<absv(A2c[cel])<<","<<absv(W[(integer)indexes[col]]-W[(integer)indexes[row]])<<";";
                //}
            }
        }
        // filling energy resolved md2
        if (ldos>0){
            dloc /= (double)ldos;
            dloclog/=(double)ldos;
            mtyp/=(double)ldos;
            md2[inde1]+=mtyp-dloclog;
            counts3[inde1]+=1;
        }
        // filling number of resonances
        counts4[inde1]+=1;
        resav[inde1] += nres;
    }
    // END OF LOOP OVER BAND!!!
    // now NORMALIZING what is nessecary 
    // dividing by the number of measurments
    for (col=0; col<band2-1; col++){
        md[band2-1]+=md[col];
        if (counts2[col]>0) md[col]/=(double)counts2[col];
        md2[band2-1]+=md2[col];
        if (counts3[col]>0) md2[col]/=(double)counts3[col];
        resav[band2-1]+=resav[col];
        if (counts4[col]>0) resav[col]/=(double)counts4[col];
        //cout << md[col]<<"; "<< counts2[col]<<"; ";
    }
    //cout<< endl;
    md[band2-1]/= (double) hdim-1;
    md2[band2-1]/= (double) hdim-1;
    resav[band2-1]/= (double) hdim-1;
    
    for (col=0; col<nbins; col++){
        if (countsR[col]>0) resf[col]/=(double)countsR[col];
        //cout<<resf[col]<<",";
    }
    //cout<<endl;
    delete [] counts4;
    delete [] countsR;
/*    for (col=0; col<nbins; col++){
    //cout<<mdc[col]<<"/"<<counts4[col]<<endl;
//        if (counts4[col]>0) mdc[col]=((double) mdc[col])/(double)counts4[col];
        //cout<<mdc[col]<<endl;
    }
    // printing mdde for a given disorder
  cout<<"mde("<<crun+1<<",:)=["<<mdde[0];
    for (col=1; col<nbins; col++){
        cout<<","<<mdde[col];
    }
    cout<<"];"<<endl;
    cout<<"mddistr("<<crun+1<<",:)=["<<mddistr[0];
    for (col=1; col<nbins; col++){
        cout<<","<<mddistr[col];
    }
    cout<<"];"<<endl;
    cout<<"mdc("<<crun+1<<",:)=["<<mdc[0];
    for (col=1; col<2; col++){
        cout<<","<<mdc[col];
    }
    cout<<"];"<<endl;*/
}


// fills in the histogram of matrix elements and adjacent energy levels
void obs_band::fillhist2fast(double *h2, double *md, int *counts2, double *md2, int *counts3,double *mddistr){
    integer row, col, ind, inde1,inde2;
    double inc  = 1./(hdim-1);
    // reset h1
    for (col=0; col<nbins; col++){
        h2[col]=0.;
        mddistr[col]=0.;
    }
    for (col=0; col<band2; col++){
        counts2[col]=0;
        counts3[col]=0;
        md[col]=0;
        md2[col]=0;
    }
    for (col=0; col<hdim-1; col++) { // loop over vectors for these states
        row=col+1;
        if (absv(A2c[col])>1e-15){
            ind = (int)floor(-log(absv(A2c[col]))/30*nbins);
        } else {
            ind = nbins-1;
        }
        //filling matrix elements for nearest neighbors
        if (ind<nbins){
            h2[ind] += inc;
        } else h2[nbins-1] +=inc;
        // filling md scaled in the band (over 1/2 width)
        inde1 =(int) (((theham->W[col])/(3*hw)+0.5)*(band2-1));
        if (inde1>=band2-1) inde1=band2-2;
        if (inde1<0) inde1=0;
        md[inde1]-=log(absv(A2c[col])/(theham->W[row]-theham->W[col]));
        counts2[inde1]+=1;
        // filling md2 vs energy density
        inde1 = (int) (((theham->W[col])/((double)nspin)/3.+0.5)*(band2-1));
        if (inde1>=band2-1) inde1=band2-2;
        if (inde1<0) inde1=0;
        md2[inde1]-=log(absv(A2c[col])/(theham->W[row]-theham->W[col]));
        counts3[inde1]+=1;
        // filling histogram if we are in the middle
        if ((row>=ez-dt/2)&(row<ez+dt/2)){
            inde2=(log(absv(A2c[col])/(theham->W[row]-theham->W[col]))+30)/60*nbins;
            if (inde2>=nbins) inde2=nbins-1;
            if (inde2<0) inde2=0;
            mddistr[inde2]+=1.;
        }
        
    }
    // dividing by the number of measurments
    for (col=0; col<band2-1; col++){
        md[band2-1]+=md[col];
        if (counts2[col]>0) md[col]/=(double)counts2[col];
        md2[band2-1]+=md2[col];
        if (counts3[col]>0) md2[col]/=(double)counts3[col];
        //cout << md[col]<<"; "<< counts2[col]<<"; ";
    }
    //cout<< endl;
    md[band2-1]/= (double) hdim-1;
    md2[band2-1]/= (double) hdim-1;
}



//fills in distribution of matrix elements with energy difference in the middle of the band: power: linear, square spectral function
// elements are filled for the middle of the band
void obs_band::fillae(double *ae0,double *ae1,double *aeoE,double *ae0L,double *aeoEL,double *aew){
    int col, row, ind;
    integer cel;
    //double ae1n;
    int *cntr = new int[nbins];
    for (row=0; row<nbins; row++) {
        ae0[row] = 0.;
        ae1[row] = 0.;
        aeoE[row] = 0.;
        ae0L[row] = 0.;
        aeoEL[row] = 0.;
        aew[row] = 0.;
        cntr[row] = 0;
    }
    //cout<<"matrix elements with energy"<<endl;
    //cout<<"mat=[";
    //ae1n = 0;
    double mel,mel2;
    double melL;
    for (col=ez-dt/2; col<ez+dt/2; col++) {
        for (row=0; row<hdim; row++) {
            cel = row+hdim*col;
            if (row!=col){
                ind = (int)((theham->W[row]-theham->W[col]+2.5)/5.*(nbins));
                if (ind<0) ind = 0;
                if (ind>=nbins) ind = nbins-1;
                mel  = absv(A2c[cel]);
                mel2 = pow(mel,2);
                melL = log(mel);
                ae0[ind] += mel;
                ae1[ind] += mel2;
                aeoE[ind] += mel2/absv(theham->W[row]-theham->W[col]);
                ae0L[ind] += melL;
                aeoEL[ind] += 2*melL - log(absv(theham->W[row]-theham->W[col]));
                aew[ind] += mel2;
                cntr[ind]+=1;
                //ae1n += pow(A2c[row+hdim*col],2);
            }
            //cout<<theham->W[row]-theham->W[col]<<","<<pow(A2c[row+hdim*col],2)<<";";
        }
    }
    for (row=0; row<nbins; row++) {
        if (cntr[row]>0){
            ae0[row] /= (double)cntr[row];
            ae1[row] /= (double)cntr[row];
            aeoE[row] /= (double)cntr[row];
            ae0L[row] /= (double)cntr[row];
            aeoEL[row] /= (double)cntr[row];
        }
    }
    delete[] cntr;
    //cout<<"];";
}



// measuring distribution of the matrix elements for all operators
void obs_band::measureAstat2(){
    //loop over spins
    int cind=0;
    int cspin,cspin2;
    // sz in the beginning
    cspin=1;
    cspin2=1;
    cind = crun;
    szA2vec(cspin);    //szmA2(cspin);
    buildA2c();
    if (STAT){
        fillae(ae0[cind],ae1[cind],aeoE[cind],ae0L[cind],aeoEL[cind],aew[cind]);
        fillhist2sort(hyst1[cind],hyst2[cind],mdE[cind],cE[crun],mdE2[cind],cE2[crun],mdd[cind],mdc[cind],mdce[cind]);
        //if ((crun==0)&&(cspin==1)){
        //    memcpy(snapv1[0], A2c+(int)ez*hdim, sizeof(doublereal)*hdim);
        //}
    }else{
        fillhist2fast(hyst2[cind],mdE[cind],cE[crun],mdE2[cind],cE2[crun],mdd[cind]);
    }
    
    // spsm in the beginning
    cind = nruns+crun;
    spsmmA2vec(cspin2);
    buildA2c();
    if (STAT){
        fillae(ae0[cind],ae1[cind],aeoE[cind],ae0L[cind],aeoEL[cind],aew[cind]);
        fillhist2sort(hyst1[cind],hyst2[cind],mdE[cind],cE[crun],mdE2[cind],cE2[crun],mdd[cind],mdc[cind],mdce[cind]);
        //if ((crun==0)&&(cspin==1)){
        //    memcpy(snapv1[1], A2c+(int)ez*hdim, sizeof(doublereal)*hdim);
        //}
    }else{
        fillhist2fast(hyst2[cind],mdE[cind],cE[crun],mdE2[cind],cE2[crun],mdd[cind]);
    }
    /*
    // current
    cspin=nspin2;
    cspin2=nspin2;
    cind = 2*nruns+crun;
    currA2();    //szmA2(cspin);
    buildA2c();
    if (STAT){
        //fillae(ae1[cind]);
        fillhist2(hyst1[cind],hyst2[cind],mdE[cind],cE[crun],mdE2[cind],cE2[crun],mdd[cind],mdc[cind],mdce[cind]);
        //if ((crun==0)&&(cspin==1)){
        //    memcpy(snapv1[0], A2c+(int)ez*hdim, sizeof(doublereal)*hdim);
        //}
    }else{
        fillhist2fast(hyst2[cind],mdE[cind],cE[crun],mdE2[cind],cE2[crun],mdd[cind]);
    }
    // spsm-spsm at opposite ends
    cind = 3*nruns+crun;
    sp__smMA2();
    buildA2c();
    if (STAT){
        //fillae(ae1[cind]);
        fillhist2(hyst1[cind],hyst2[cind],mdE[cind],cE[crun],mdE2[cind],cE2[crun],mdd[cind],mdc[cind],mdce[cind]);
        //if ((crun==0)&&(cspin==1)){
        //    memcpy(snapv1[1], A2c+(int)ez*hdim, sizeof(doublereal)*hdim);
        //}
    }else{
        fillhist2fast(hyst2[cind],mdE[cind],cE[crun],mdE2[cind],cE2[crun],mdd[cind]);
    }
     */
}


// QUANTITIES I

// calculate entanglement entropy for several positions as well as energy fluctuations for a half-system!
void obs_band::ententposE(){
    double *ent;
    ent = new double[dt];
    double inc = 1./((double)dt);
    int hind;
        // zero out hyst
    for (int numv=0;numv<nbins; numv++) {
        hystdmte[crun][numv] = 0;
        hystspace[crun][numv] = 0;
    }
    for (int numv=0;numv<nspin2; numv++) {
        sentl[crun][numv] = 0;
    }
    int cspin;
    //cout<<"Sent={"<<endl;
    for (int i=ez-dt/2; i<ez+dt/2; i++) { // loop over eigenstates
        thedmt->constructrho(i,theham);
        // histogram of the eigenvalues of the density matrix
        thedmt->svd();
          //cout<< "W("<<i+1<<",:)=[";
        for (int j=0; j<hdimh; j++) {
            if (absv(thedmt->W[j])>1e-17)
            {
                hind = (int)((-log(absv(thedmt->W[j])))/35.*(nbins));
                if (hind<0) hind = 0;
                if (hind>=nbins) hind = nbins-1;
            } else hind = nbins-1;
            hystspace[crun][hind]+=1;
            //cout<< thedmt->W[j]<<", ";
        }
        //cout<<"];"<<endl;
        // histogram of the entanglement entropy
        ent[i-ez+dt/2] = thedmt->entent();
        //cout<<i-ez+dt/2<<"; "<< ent[i-ez+dt/2]<<endl;
        hind  = (int)floor(ent[i-ez+dt/2]/5.*(nbins));
        if (hind>nbins-1) hind = nbins-1;
        hystdmte[crun][hind] += 1;
        if (crun<2) snapent[crun][i-ez+dt/2]=ent[i-ez+dt/2];
        //looping to fill sentl - entanglement with L
        for (cspin=0; cspin<nspin2; cspin++) {
            thedmt->constructrho(i, theham, cspin+1);
            thedmt->svd(cspin+1);
            sentl[crun][cspin]+=thedmt->entent(cspin+1)*inc;
        }
    }
    //averaging and storing
    getav(ent, dt, entropypos[crun]);
    //cout<<"dt="<<dt<<"; "<< entropypos[crun][0]<<endl;
    delete[] ent;
}

//helper for calculating level statistics and histogram
double obs_band::levelstatband(int si, int ei,int * hst){
    double r=0.;
    double d1, d2;
    double bin = sqrt((double)nspin)*50/((double)hdim);
    int ind;
    for (int i=0; i<nbins; i++) {
        hst[i] = 0;
    }
    d1 = theham->W[si+1]-theham->W[si];
    for (int i=si+1; i<ei; i++) {
        //filling the histogram
        ind = (int)(d1/bin*(nbins));
        if (ind>=nbins) ind = nbins-1;
        hst[ind]+=1;
        d2 = theham->W[i+1]-theham->W[i];
        r += (d2>d1)?d1/d2:d2/d1;
        d1 = d2;
    }
    r = r/(double(ei-si-1));
    return r;
}

//helper for calculating just level statistics
double obs_band::levelstatband(int si, int ei){
    double r=0.;
    double d1, d2;
    double bin = sqrt((double)nspin)*50/((double)hdim);
    int ind;
    d1 = theham->W[si+1]-theham->W[si];
    for (int i=si+1; i<ei; i++) {
        //filling the histogram
        ind = (int)(d1/bin*(nbins));
        if (ind>=nbins) ind = nbins-1;
        d2 = theham->W[i+1]-theham->W[i];
        r += (d2>d1)?d1/d2:d2/d1;
        d1 = d2;
    }
    r = r/(double(ei-si-1));
    return r;
}

//calculate level statistics 
void obs_band::levelstat(){
    int i=0;
    lspos[crun][i] = levelstatband(ez-dt/2-1,ez+dt/2,lshist[crun+i*nruns]);
    i=1;
    int third=(int)(hdim/6);
    lspos[crun][i] = levelstatband(ez-third-1,ez+third);
}

void obs_band::findez(){
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

    hw = theham->W[indup]-theham->W[inddo];
    //cout<<"% hw="<<hw<<endl;
}

void obs_band::measuredata(){
/* // testing Sent
    int s1,s2;
    bool flag;
    double *vec = new double[hdim];
    for (int i=0; i<hdim; i++) {
        flag=true;
        for (int j=0; j<nspin2; j++) {
            s1=(int)(thehsp->hspace[i]/(int)pow(2,2*j))%2;
            s2=(int)(thehsp->hspace[i]/(int)pow(2,2*j+1))%2;
            if (s1+s2!=1) flag=false;
        }
        if (flag) vec[i]=pow(0.5,0.25*(double)nspin);
        else vec[i]=0;
    }
    double norm=0;
    for (int i=0; i<hdim; i++) {
        norm+=pow(vec[i],2);
    }
    cout<<"test vector prepared, norm="<<norm<<endl;
    for (int i=1; i<nspin2+1; i++) {
    cout<<"L="<<i<<endl;
    thedmt->constructrho(vec, i);
//    thedmt->printA();
    thedmt->svd(i);
//    thedmt->printW();
    cout<<"sent"<<thedmt->entent(i)<<endl;
    }
 */
    // measure quantities -I
    findez();
    //ententposE();
    //levelstat();
    measureAstat2();
    // increase number of current runs by 1
    crun += 1;
}



// writing data to file=filep, which must be initialized!!!
void obs_band::writedata() {
    double *darr1, *darr2;
    darr1 = new double[nspin];
    darr2 = new double[nspin];
    // averaging first the data and saving it into array to be dumped for all quantities
    
    getav(entropypos, crun, band, darr1, darr2);
    fput("Ent entropy in the band: av", darr1, band);
    fput("Ent entropy in the band: disp", darr2, band);
    
    getav(lspos, crun, 2, darr1, darr2);
    fput("levelspacings: av", darr1, 2);
    fput("levelspacings: disp", darr2, 2);
    
    getav(sentl, crun, nspin2, darr1, darr2);
    fput("Ent entropy with size: av", darr1, nspin2);
    fput("Ent entropy with size: disp", darr2, nspin2);
    
    double *darr3, *darr4;
    darr3 = new double[nbins];
    darr4 = new double[nbins];
    
    for (int i=0; i<band; i++) {
        getav(lshist, crun+nruns*i, nbins, darr3, darr4);
        fput("hist for levelstat", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    for (int i=0; i<band; i++) {
        getav(hystdmte, crun+nruns*i, nbins, darr3, darr4);
        fput("hist for entanglement entropy", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    for (int i=0; i<band; i++) {
        getav(hystspace, crun+nruns*i, nbins, darr3, darr4);
        fput("hist for dmt levels", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    fput("Sent snapshot1", snapent[0],hdim);
    fput("Sent snapshot2", snapent[1],hdim);
    //  fput("Sent snapshot3", snapent[2],hdim);
    
    getav(hyst1, crun, nbins, darr3, darr4);
    fput("hist-all for sz at edge", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(hyst1+nruns, crun, nbins, darr3, darr4);
    fput("hist-all for spsm at edge", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(hyst1+2*nruns, crun, nbins, darr3, darr4);
    fput("hist-all for spsm in midlle", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(hyst1+3*nruns, crun, nbins, darr3, darr4);
    fput("hist-all for spsm at opposite ends", darr3, nbins);
    fput("disp", darr4, nbins);
    
    
    getav(hyst2, crun, nbins, darr3, darr4);
    fput("hist-nn for sz", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(hyst2+nruns, crun, nbins, darr3, darr4);
    fput("hist-nn for spsm", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(hyst2+2*nruns, crun, nbins, darr3, darr4);
    fput("hist-nn for spsm", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(hyst2+3*nruns, crun, nbins, darr3, darr4);
    fput("hist-nn for spsm", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(cE, crun, band2, darr3, darr4);
    fput("dos scaled", darr3, band2);
    fput("disp", darr4, band2);
    
    getav(cE2, crun, band2, darr3, darr4);
    fput("dos-vs-e", darr3, band2);
    fput("disp", darr4, band2);
    
    getav(mdE, crun, band2, darr3, darr4);
    fput("mdE = M/sqrt d", darr3, band2);
    fput("disp", darr4, band2);
    
    getav(mdE+nruns, crun, band2, darr3, darr4);
    fput("mdE for spsm", darr3, band2);
    fput("disp", darr4, band2);
    
    getav(mdE+2*nruns, crun, band2, darr3, darr4);
    fput("mdE for spsm", darr3, band2);
    fput("disp", darr4, band2);
    
    getav(mdE+2*nruns, crun, band2, darr3, darr4);
    fput("mdE for spsm", darr3, band2);
    fput("disp", darr4, band2);
    
    getav(mdE2, crun, band2, darr3, darr4);
    fput("mdE2 for Sz vs e", darr3, band2);
    fput("disp", darr4, band2);
    
    getav(mdE2+nruns, crun, band2, darr3, darr4);
    fput("mdE2 for spsm", darr3, band2);
    fput("disp", darr4, band2);
    
    getav(mdE2+2*nruns, crun, band2, darr3, darr4);
    fput("mdE2 for spsm", darr3, band2);
    fput("disp", darr4, band2);
    
    getav(mdE2+3*nruns, crun, band2, darr3, darr4);
    fput("mdE2 for spsm", darr3, band2);
    fput("disp", darr4, band2);
    
    
    getav(mdd, crun, nbins, darr3, darr4);
    fput("V/delta-distr in the middle of the band", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(mdd+nruns, crun, nbins, darr3, darr4);
    fput("V/delta-distr in the middle of the band Sz", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(mdd+2*nruns, crun, nbins, darr3, darr4);
    fput("V/delta-distr in the middle of the band Sz", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(mdd+3*nruns, crun, nbins, darr3, darr4);
    fput("V/delta-distr in the middle of the band Sz", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(mdc, crun, nbins, darr3, darr4);
    fput("prob of resonances", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(mdc+nruns, crun, nbins, darr3, darr4);
    fput("prob of resonances", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(mdc+2*nruns, crun, nbins, darr3, darr4);
    fput("prob of resonances ", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(mdc+3*nruns, crun, nbins, darr3, darr4);
    fput("prob of resonances", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(mdce, crun, band2, darr3, darr4);
    fput("Nres in the middle of the band", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(mdce+nruns, crun, band2, darr3, darr4);
    fput("Nres in the middle of the band", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(mdce+2*nruns, crun, band2, darr3, darr4);
    fput("Nres in the middle of the band", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(mdce+3*nruns, crun, band2, darr3, darr4);
    fput("Nres in the middle of the band", darr3, nbins);
    fput("disp", darr4, nbins);
    
    /*
     fput("ae-Mhist for 0th obs snapshot", ae1[0],nbins);
     fput("ae-Mhist for 1th obs snapshot", ae1[nruns*2],nbins);
     //fput("ae-Mhist for 2th obs", ae1[2*nruns],nbins);
    */
    getav(ae0, crun, nbins, darr3, darr4);
    fput("M average for Sz", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(ae0+nruns, crun, nbins, darr3, darr4);
    fput("M^ average for SpSm", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(ae1, crun, nbins, darr3, darr4);
     fput("M^2 average for Sz", darr3, nbins);
     fput("disp", darr4, nbins);
    
    getav(ae1+nruns, crun, nbins, darr3, darr4);
     fput("M^2 average for SpSm", darr3, nbins);
     fput("disp", darr4, nbins);

    getav(aeoE, crun, nbins, darr3, darr4);
    fput("M^2/E average for Sz", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(aeoE+nruns, crun, nbins, darr3, darr4);
    fput("M^2/E average for SpSm", darr3, nbins);
    fput("disp", darr4, nbins);

    getav(ae0L, crun, nbins, darr3, darr4);
    fput("Log(M) average for Sz", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(ae0L+nruns, crun, nbins, darr3, darr4);
    fput("log(M) average for SpSm", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(aeoEL, crun, nbins, darr3, darr4);
    fput("log(M^2/E) average for Sz", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(aeoEL+nruns, crun, nbins, darr3, darr4);
    fput("log(M^2/E) average for SpSm", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(aew, crun, nbins, darr3, darr4);
    fput("M^4 average for Sz", darr3, nbins);
    fput("disp", darr4, nbins);
    
    getav(aew+nruns, crun, nbins, darr3, darr4);
    fput("M^4 average for SpSm", darr3, nbins);
    fput("disp", darr4, nbins);

    
     /*
    //getav(ae1+2*nruns, crun, nbins, darr3, darr4);
     //fput("ae av S.S", darr3, nbins);
     //fput("disp", darr4, nbins);
     
     fput("snapshot Sz", snapv1[0],hdim);
     fput("snapshot Spsm", snapv1[1],hdim);
     //fput("snapshot SzSz", snapv1[2],hdim);
     
     */
    
    
    
    delete[] darr1;
    delete[] darr2;
    delete[] darr3;
    delete[] darr4;
    
    
    
    
}

