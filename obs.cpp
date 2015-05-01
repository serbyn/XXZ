#include <iostream>
#include <time.h>
#include <stdio.h>
#include <ran.h>
#include <obs.h>
#include <hsp.h>
#include <ham.h>
#include <fstream>
#include <sstream>
#include <memory.h>
#include <globalfunctions.h>

obs::obs(hsp *thehsp1,ham* theham1, Ran* theran1, int nruns1, FILE* filep1)
{
    thehsp = thehsp1;
    theham = theham1;
    theran = theran1;
    //cout << "Constructing observables object ...\n";
    nspin = thehsp->nspin;
    nspin2 = nspin/2;
    hdim0 = thehsp->hdim0;
    hdim = thehsp->hdim;
    hdim2 = thehsp->hdim2;
    hdimh = thehsp->hdimh;
    nruns = nruns1;
    filep = filep1;
    crun = 0;
}

void obs::stringput(string comment, double num){
    sline1<<comment<<", ";
    sline2<<num<<", ";
}

void obs::stringput(string comment, int num){
    sline1<<comment<<", ";
    sline2<<num<<", ";
}

void obs::stringwrite(){
    string str = sline1.str();
    str.erase (str.end()-2, str.end());
    fprintf(filep,"%s \n",&str[0]);
    str = sline2.str();
    str.erase (str.end()-2, str.end());
    fprintf(filep,"%s \n",&str[0]);
    // cleaning strings up!!!
    sline1.str("");
    sline2.str("");
}

// helper functions for writing data to file
void obs::fput(string comment, double* array, int arrl) {
    fprintf(filep,"# %s \n",&comment[0]);
    for(int i=0;i<arrl-1;i++){
        fprintf(filep, "%.15g, ", array[i] );
    }
    fprintf(filep, "%.15g \n ", array[arrl-1] );
    //fprintf(filep,"\n");
}

void obs::fput(string comment,int n) {
    fprintf(filep,"# %s \n",&comment[0]);
    fprintf(filep,"%i \n", n);
    //fprintf(filep,"\n");
}

void obs::fput(string comment,double n) {
    fprintf(filep,"# %s \n",&comment[0]);
    fprintf(filep,"%.15g \n", n);
    //fprintf(filep,"\n");
}

void obs::writebasicdata(){
    stringput("nspin",nspin);
    stringput("nruns", crun);
    stringwrite();
    // write vector quantities
    
}

obs::~obs() {
}
