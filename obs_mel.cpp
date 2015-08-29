#include <iostream>
#include <time.h>
#include <stdio.h>
#include <ran.h>
#include <obs.h>
#include <obs_mel.h>
#include <hsp.h>
#include <ham.h>
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

obs_mel::obs_mel(hsp *thehsp1, ham* theham1, Ran* theran1, int nruns1, FILE* filep1): obs(thehsp1,theham1,theran1,nruns1,filep1)
{
    band2 = 101;
    dt = 2*( (int) (hdim/nspin/2));
    if (dt<2) dt=2;
    nbins = 1000;
    nobs = 3;
    ez=hdim/2;
    hw = 1;
    cout<<"% in the middle we have "<< dt <<" states to study"<<endl;
    // QUANTITIES II
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
    mdce2 = createdouble(nobs*nruns,band2);
    cE = createint(nruns,band2); // average mdE as a function of relative energy
    cE2 = createint(nruns,band2); // average mdE as a function of absolute energ
    ae0 = createdouble(nruns*nobs,nbins);
    ae1 = createdouble(nruns*nobs,nbins);
    aeoE = createdouble(nruns*nobs,nbins);
    ae0L = createdouble(nruns*nobs,nbins);
    aeoEL = createdouble(nruns*nobs,nbins);
    aew = createdouble(nruns*nobs,nbins);
    snapv1 = createdouble(nobs,hdim);
    beta = createdouble(nobs*nruns,nbins);
    betaL = createdouble(nobs*nruns,nbins);
    m = createdouble(nobs*nruns,nbins);
    mL = createdouble(nobs*nruns,nbins);
    moES = createdouble(nobs*nruns,nbins);
    moESL = createdouble(nobs*nruns,nbins);
    iprM = new double[nobs*nruns];
    iprMstat = createdouble(nobs*nruns,nbins);
}


obs_mel::~obs_mel() {
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
    destroy(mdce2,nobs*nruns);
    destroy(beta,nobs*nruns);
    destroy(betaL,nobs*nruns);
    destroy(m,nobs*nruns);
    destroy(mL,nobs*nruns);
    destroy(moES,nobs*nruns);
    destroy(moESL,nobs*nruns);
    delete[] iprM;
    destroy(iprMstat,nobs*nruns);
}
// QUANTITIES II-1
// creating matrix A2 = 2sz A, where sz acts on the largest index or rightmost spin
void obs_mel::szA2(){
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
void obs_mel::szA2vec(int ns){
    integer row;
    int t = (int)pow(2,nspin-ns);
    memcpy(A2, theham->A, sizeof(double)*hdim2); //start with the matrix A
    for (row=0; row<hdim; row++) {
        if ((thehsp->hspace[row]/t)%2==0){
            cblas_dscal(hdim, -1., A2+row, hdim);
        }
    }
}

void obs_mel::szmA2(int ns){
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
void obs_mel::szszA2(){
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
void obs_mel::spsmA2(){
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
void obs_mel::spsmmA2(int ns){
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
void obs_mel::sp__smA2(){
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


// creating matrix A2 = (S+S- - S-S+) A where S+ and S- are at site 1 and nspin2+1
void obs_mel::sp__smMA2(){
    integer col,row,tmp,s1,s2,b,row2;
    int t = (int)pow(2,nspin2);
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

// creating matrix A2 = j A where j is current operator on firts 2 sites
void obs_mel::curr1A2(int spinpos){
    integer col,row,tmp,s1,s2,b,row2;
    int t;
    if ((spinpos<0)||(spinpos>nspin-1)){
        cout<<"obs_mel::curr1A2(int spinpos):: error: illegal value of spinpos="<<spinpos<<endl;
        exit(-1);
    }
    //memcpy(A2, theham->A, sizeof(double)*hdim2); //start with the matrix A
    for(col=0; col<hdim2; col++) {
        A2[col]=0.;
    }
    // loop over spin position before end
    if(spinpos<nspin-1) {
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
    } else{ // assuming that spin operator goes across the cut
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
}

// creating matrix A2 = j A where j is current operator
void obs_mel::currA2(){
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
void obs_mel::sp__smA2vec(){
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
void obs_mel::spsmmA2vec(int ns){
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
void obs_mel::ssA2(){
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
void obs_mel::buildA2c(){
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
/*
template<typename T> class CompareIndicesByAnotherVectorValues {
    std::vector<T>* _values;
    public: CompareIndicesByAnotherVectorValues(std::vector<T>* values): _values(values) {}
    public: bool operator() (const int& a, const int& b) const
    {
        return (*_values)[a] > (*_values)[b];
    }
};

template <class RAIter, class Compare>
void argsort(RAIter iterBegin, RAIter iterEnd, Compare comp,
             std::vector<size_t>& indexes) {
    
    std::vector< std::pair<size_t,RAIter> > pv ;
    pv.reserve(iterEnd - iterBegin) ;
    
    RAIter iter ;
    int k ;
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
*/


typedef std::pair<int,double> mypair;
bool comparator ( const mypair& l, const mypair& r)
{ return l.second < r.second; };

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
void obs_mel::fillhist2sort(double *h1, double *h2, double *md, int *counts2, double *md2, int *counts3,double *mddistr,double *resf,double *
    resav,double *resav2){
    integer row, col, ind,ind2,inde,inde1, cel;
    double inc  = 1./(hdim);
    double Z=0.25; // factor over which we multiply matrix element
    int *counts4 =  new int[band2];
    int *counts5 =  new int[band2];
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
        counts5[col]=0;
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
    pair<int,double> *mypair = new pair<int,double>[hdim];
    for (row=0; row<hdim; row++) {
        mypair[row].first = row;
        mypair[row].second = W[row];
        //cout<< mypair[ii].first<<","<< mypair[ii].second<<endl;
    }
    std::sort(mypair, mypair+hdim, comparator);
    int* indexes=new int[hdim];
    for (row=0; row<hdim; row++) {
        //cout<< mypair[row].first<<","<< mypair[row].second<<endl;
        indexes[row] =mypair[row].first;
    }
   // argsort(W, W + hdim, std::less<double>(), indexes);
    double bw = theham->W[hdim-1]-theham->W[0];
    double bw2 = 0.5*(theham->W[hdim-1]-theham->W[0]);
    double e0 =theham->W[0];
    double ec = theham->W[ez];
    double moverd;
    int nres;
    integer irow,icol;
    for (row= 0; row<hdim-1; row++) { // loop over individual states within the band
        irow = (integer)indexes[row];
        // index for the energy histogram like Allet
        inde = (int) (((theham->W[irow]-e0)/(bw))*(band2-1));
        if (inde>=band2-1) inde=band2-2;
        if (inde<0) inde=0;
        // index for the energy histogram w.r.t. the band center
        inde1 = (int) (((theham->W[irow]-ec)/bw2+0.5)*(band2-1));
        if (inde1>=band2-1) inde1=band2-2;
        if (inde1<0) inde1=0;
        // measuring nn probes:
        col = row+1;
        icol = (integer)indexes[col];
        cel = irow+((integer)hdim)*icol;
        if (absv(A2c[cel])>1e-15){
            ind = (int)floor(-log(absv(A2c[cel]))/30*nbins);
        } else {
            ind = nbins-1;
        }
        if (ind>=nbins) ind = nbins-1;
        if (ind<0) ind=0;
        // distribution of nn spacings
        h1[ind] += inc;
        moverd = log(absv(A2c[cel]/(W[irow]-W[icol])));
        md[inde]+=moverd;
        counts2[inde]++;
        md2[inde1]+=moverd;
        counts3[inde1]++;
        // distribution if we are in the middle
        if ((irow>=ez-dt/2)&(irow<ez+dt/2)){
            ind2=(moverd+30)/60*nbins;
            if (ind2>=nbins) ind2=nbins-1;
            if (ind2<0) ind2=0;
            mddistr[ind2]+=1.;
        }
        // taking care of resonances for given level
        nres=0;
        for (col=0; col<hdim; col++) { // loop over vectors for these states; we exclude symmetric copy
            if (col!=row){
                icol = (integer)indexes[col];
                cel = irow+((integer)hdim)*icol;
                // calculating resonances (energy resolved)
                if (absv(A2c[cel])>absv(W[icol]-W[irow])){
                    nres++;
                }
                // probability of resonance for the middle -- over all levels
                if((row>=ez-dt/2)&(row<ez+dt/2)){
                    ind2  = abs(row-col);
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
        }
        // filling number of resonances
        counts4[inde1]++;
        resav[inde1] += nres;
        counts5[inde1]++;
        resav2[inde1] += nres;
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
        resav2[band2-1]+=resav2[col];
        if (counts5[col]>0) resav2[col]/=(double)counts5[col];
        
        //cout << md[col]<<"; "<< counts2[col]<<"; ";
    }
    //cout<< endl;
    md[band2-1]/= (double) hdim-1;
    md2[band2-1]/= (double) hdim-1;
    resav[band2-1]/= (double) hdim-1;
    resav2[band2-1]/= (double) hdim-1;
    
    for (col=0; col<nbins; col++){
        if (countsR[col]>0) resf[col]/=(double)countsR[col];
        //cout<<resf[col]<<",";
    }
    //cout<<";resf="<<endl;
    delete [] counts4;
    delete [] counts5;
    delete [] countsR;
   /* cout<<"mddistr("<<crun+1<<",:)=["<<mddistr[0];
    for (col=1; col<10; col++){
        cout<<","<<resf[col];
    }
    cout<<"];"<<endl;*/
    delete[] mypair;
    delete[] indexes;
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
void obs_mel::fillhist2fast(double *h2, double *md, int *counts2, double *md2, int *counts3,double *mddistr){
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

void obs_mel::fillIPR(double* iprM, double *iprMstat){
    int col, row, ind;
    integer cel;
    iprM[0]=0.;
    for (row=0; row<nbins; row++) {
        iprMstat[row]=0;
    };
    matrix iprmat(hdim,false);
    memcpy(iprmat.A, A2c, hdim2*sizeof(double));
    for (row=0; row<hdim; row++) {
        iprmat.A[row+hdim*row] += theham->W[row];
    }
    iprmat.diagonalize();
    double cef;
    int ni=0;
    double iprl;
    for (col=ez-dt/2; col<ez+dt/2; col++) {
        iprl=0;
        for (row=0; row<hdim; row++) {
            cel = row+hdim*col;
            cef= absv(iprmat.A[cel]);
            if (cef>1e-15){
                ind = (int)floor(-log(cef)/30*nbins);
            } else {
                ind = nbins-1;
            }
            if (ind<0) ind =0;
            iprMstat[ind]+=1;
            iprl+=pow(cef,4);
        }
        //cout<<iprl<<endl;
        ni++;
        iprM[0]+=1./iprl;
    }
    iprM[0]/=(double)ni;
}

//fills in distribution of matrix elements with energy difference in the middle of the band: power: linear, square spectral function
// elements are filled for the middle of the band
// aeoE now is just diagonal matrix element with energy difference!!!
void obs_mel::fillae(double *ae0,double *ae1,double *aeoE,double *ae0L,double *aeoEL,double *aew){
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
    double mel,mel2,mel0,mel00;
    double melL;
    for (col=ez-dt/2; col<ez+dt/2; col++) {
        mel0 = A2c[col+hdim*col];
        for (row=0; row<hdim; row++) {
            cel = row+hdim*col;
            mel  = absv(A2c[cel]);
            mel00  = A2c[row+hdim*row];
            mel2 = pow(mel,2);
            // filling only positive energy part:
            if (row>col){
                ind = (int)((theham->W[row]-theham->W[col])/5.*(nbins));
                if (ind<0) ind = 0;
                if (ind>=nbins) ind = nbins-1;
                melL = log(mel);
                ae0[ind] += mel;
                ae1[ind] += mel2;
                aeoE[ind] += mel0*mel00;
                aeoEL[ind] += log(absv(mel0*mel00));
                ae0L[ind] += melL;
                cntr[ind]+=1;
                //ae1n += pow(A2c[row+hdim*col],2);
            }
            // filling all energies for spectral function
            ind = (int)((theham->W[row]-theham->W[col]+1.5)/3.*(nbins));
            if (ind<0) ind = 0;
            if (ind>=nbins) ind = nbins-1;
            aew[ind] += mel2; // spectral function includes contribution at zero energy
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
    if (iter ==0)     cout<<"awSz(";
    else if (iter==1) cout<<"awSpSm(";
    else if (iter==2) cout<<"awJ(";
    else    cout<<"aw(";
    cout<<crun+1<<",:)=["<<aew[0];
    for (col=1; col<nbins; col++){
        cout<<","<<aew[col];
    }
    cout<<"];"<<endl;
    delete[] cntr;
    //cout<<"];";
}

//fills in distribution of matrix elements with index in the middle of the band
// moES now is just diagonal matrix element with energy difference!!!
void obs_mel::fillaN(double *beta,double *betaL,double *m,double *mL,double *moeS,double *moeSL){
    int col, row, ind, ind2;
    integer cel;
    //double ae1n;
    int *cntr = new int[nbins];
    for (row=0; row<nbins; row++) {
        beta[row] = 0.;
        betaL[row] = 0.;
        m[row] = 0.;
        mL[row] = 0.;
        moeS[row] = 0.;
        moeSL[row] = 0.;
        cntr[row] = 0;
    }
    //cout<<"matrix elements with energy"<<endl;
    //cout<<"mat=[";
    //ae1n = 0;
    double mel,mel2,mel0,mel00;
    double melL;
    for (col=ez-dt/2; col<ez+dt/2; col++) {
        mel0 = A2c[col+hdim*col];
        for (row=0; row<hdim; row++) {
            cel = row+hdim*col;
            mel  = absv(A2c[cel]);
            mel2 = pow(mel,2);
            //mdiagL = log(absv(A2c[col+hdim*col]));
            //mdiag2 = pow(A2c[col+hdim*col],2);
            melL = log(mel);
            ind = abs(row-col);
            //if (ind<0) ind = 0;
            if (ind>=nbins) ind = nbins-1;
            //if (ind==1)  cout<<mel2<<"/"<<mdiag2<<"="<<mel2/mdiag2<<endl;
            ind2 = (int)floor(mel*nbins);
            if (ind2<0) ind2=0;
            if (ind2>=nbins) ind2 = nbins-1;
            beta[ind2] += 1;
            ind2 = (int)floor(-melL/30*nbins);
            if (ind2<0) ind2=0;
            if (ind2>=nbins) ind2 = nbins-1;
            betaL[ind2] += 1;
            m[ind] +=mel;
            mL[ind] +=melL;
            cntr[ind]++;
            mel00=mel0*A2c[row+hdim*row];
            moeS[ind]+=mel00;
            moeSL[ind]+=log(absv(mel00));
        }
    }
    for (row=0; row<nbins; row++) {
        if (cntr[row]>0){
            m[row] /= (double)cntr[row];
            mL[row] /= (double)cntr[row];
            moeSL[row] = moeSL[row]/(double)cntr[row];
            moeS[row] /= (double)cntr[row];
        }
    }
    delete[] cntr;
}




// measuring distribution of the matrix elements for all operators
void obs_mel::measureAstat2(){
    //loop over spins
    int cind=0;
    int cspin,cspin2;
    // sz in the beginning
    for (iter=0;iter<nobs;iter++){
	    cspin=1;
	    cspin2=1;
	    cind = crun+nruns*iter;
	    if (iter==0)  szA2vec(cspin);    //szmA2(cspin);
	    else if (iter==1) spsmmA2vec(cspin2);
	    else if (iter==2) curr1A2(0);
	    else if (iter==3) currA2();// sp__smMA2();
   	    buildA2c();
  	    fillae(ae0[cind],ae1[cind],aeoE[cind],ae0L[cind],aeoEL[cind],aew[cind]);
        fillaN(beta[cind],betaL[cind],m[cind],mL[cind],moES[cind],moESL[cind]);
	    fillhist2sort(hyst1[cind],hyst2[cind],mdE[cind],cE[crun],mdE2[cind],cE2[crun],mdd[cind],mdc[cind],mdce[cind],mdce2[cind]);
        //fillIPR(iprM+cind,iprMstat[cind]);
	}
    // for (int iter=0;iter<nobs*nruns;iter++){
    //     cout<<iprM[iter]<<",";
    // }
    //cout<<endl;
     /*
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


void obs_mel::findez(){
    for (int i=0; i<hdim; i++) {
        if(theham->W[i]>0){
            ez=i;
            break;
        }
    }
    //cout<<"%found zero energy at="<<ez<<endl;
    if ((ez-dt/2<0)||(ez+dt/2>=hdim)){
        cout<< "obs_mel::error -- too large dt, choosing ez to be hdim/2"<<endl;
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

void obs_mel::measuredata(){
    findez();
    measureAstat2();
    // increase number of current runs by 1
    crun += 1;
}



// writing data to file=filep, which must be initialized!!!
void obs_mel::writedata() {    
    double *darr1, *darr2;
    darr1 = new double[nspin];
    darr2 = new double[nspin];
    double *darr3, *darr4;
    darr3 = new double[nbins];
    darr4 = new double[nbins];
    // averaging first the data and saving it into array to be dumped for all quantities
    for (int i=0; i<nobs; i++) {
        getav(hyst1+nruns*i, crun, nbins, darr3, darr4);
        fput("Mel over nn for operator at edge", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    
    for (int i=0; i<nobs; i++) {
        getav(hyst2+nruns*i, crun, nbins, darr3, darr4);
        fput("Mel over dE for operator at edge", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    getav(cE, crun, band2, darr3, darr4);
    fput("dos scaled", darr3, band2);
    fput("disp", darr4, band2);
    
    getav(cE2, crun, band2, darr3, darr4);
    fput("dos-vs-e", darr3, band2);
    fput("disp", darr4, band2);
    
    for (int i=0; i<nobs; i++) {
        getav(mdE+nruns*i, crun, band2, darr3, darr4);
        fput("Mnn/delta for operator at edge", darr3, band2);
        fput("disp", darr4, band2);
    }

    for (int i=0; i<nobs; i++) {
        getav(mdE2+nruns*i, crun, band2, darr3, darr4);
        fput("Mnn/delta over energy range for operator at edge", darr3, band2);
        fput("disp", darr4, band2);
    }
    
    for (int i=0; i<nobs; i++) {
        getav(mdd+nruns*i, crun, nbins, darr3, darr4);
        fput("Mnn/delta distr", darr3, nbins);
        fput("disp", darr4, nbins);
    }

    for (int i=0; i<nobs; i++) {
        getav(mdc+nruns*i, crun, nbins, darr3,darr4);
        fput("prob of resonance", darr3, nbins);
        fput("disp", darr4, nbins);
    }

    for (int i=0; i<nobs; i++) {
        getav(mdce+nruns*i, crun, band2, darr3, darr4);
        fput("Nres", darr3,band2);
        fput("disp", darr4, band2);
    }

    for (int i=0; i<nobs; i++) {
        getav(mdce2+nruns*i, crun, band2, darr3, darr4);
        fput("Nres2", darr3,band2);
        fput("disp", darr4, band2);
    }

    
    for (int i=0; i<nobs; i++) {
        getav(ae0+nruns*i, crun, nbins, darr3, darr4);
        fput("<M> vs Energy", darr3, nbins);
        fput("disp", darr4, nbins);
    }

    for (int i=0; i<nobs; i++) {
        getav(ae1+nruns*i, crun, nbins, darr3, darr4);
        fput("<M^2> vs Energy", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    for (int i=0; i<nobs; i++) {
        getav(aeoE+nruns*i, crun, nbins, darr3, darr4);
        fput("<M^2/E> vs Energy", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    for (int i=0; i<nobs; i++) {
        getav(ae0L+nruns*i, crun, nbins, darr3, darr4);
        fput("<log M> vs Energy", darr3, nbins);
        fput("disp", darr4, nbins);
    }

    for (int i=0; i<nobs; i++) {
        getav(aeoEL+nruns*i, crun, nbins, darr3, darr4);
        fput("<log M^2/E> vs Energy", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    for (int i=0; i<nobs; i++) {
        getav(aew+nruns*i, crun, nbins, darr3, darr4);
        fput("A(w) ", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    for (int i=0; i<nobs; i++) {
        getav(beta+nruns*i, crun, nbins, darr3, darr4);
        fput("beta vs N", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    for (int i=0; i<nobs; i++) {
        getav(betaL+nruns*i, crun, nbins, darr3, darr4);
        fput("log beta vs N", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    for (int i=0; i<nobs; i++) {
        getav(m+nruns*i, crun, nbins, darr3, darr4);
        fput("m vs N", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    for (int i=0; i<nobs; i++) {
        getav(mL+nruns*i, crun, nbins, darr3, darr4);
        fput("log m vs N", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    for (int i=0; i<nobs; i++) {
        getav(moES+nruns*i, crun, nbins, darr3, darr4);
        fput("<sum M^2/E> vs N", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    
    for (int i=0; i<nobs; i++) {
        getav(moESL+nruns*i, crun, nbins, darr3, darr4);
        fput("log sum M2/E vs N", darr3, nbins);
        fput("disp", darr4, nbins);
    }
    /*
    for (int i=0; i<nobs; i++) {
        getav(iprM+nruns*i, crun, darr3, darr4);
        fput("average IPR", darr3, 1);
        fput("disp", darr4, 1);
    }
    
    for (int i=0; i<nobs; i++) {
        getav(iprMstat+nruns*i, crun, nbins, darr3, darr4);
        fput("histogram of coefficients", darr3, nbins);
        fput("disp", darr4, nbins);
    }*/
    delete[] darr1;
    delete[] darr2;
    delete[] darr3;
    delete[] darr4;
}


