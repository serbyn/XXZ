#include <iostream>
#include <time.h>
#include <stdio.h>
#ifdef __APPLE__
#include <math.h>
#else
#include "/usr/include/math.h"
#endif
#include "hsp.h"
#include <fstream>
#include <sstream>
#include <memory.h>
#include <globalfunctions.h>
//constructing a spin-half chain of a given size
hsp::hsp( int systemSize )
{
    //system size should be multiple of 2!
    if(systemSize%2 != 0) {cout<<"error. SystemSize not multiple of 2\n"; exit(1);}
    
    nspin =systemSize;
    nspin2 = systemSize/2;
    
    //setting naive Hspace dim
    hdim0 = (integer) pow(2,nspin);
    //setting true Hspace dim
    hdim = (integer)C(nspin,nspin2);
    hdim2 = ((integer)hdim)*((integer)hdim);
    hdimh = (integer) pow(2,nspin2);
    cout <<"%Constructing hilbert space for number of spins="<<nspin<<"; naive Hilbert space dimension="<< hdim0<<"; Sz=0 sector dimension = "<<hdim<< "\n";
    a = new int[nspin];
    hspace = new int[hdim];
    dict = new int[hdim0];
    mult = new int[hdimh];
    multgen = createint(nspin2,hdim+1); // generic array of multiplicities
    makehspace();
}

// Conversion between state and number: binary representation is reversed, in the beginning of array we have smallest numbers
int hsp::array_to_int( int *arr) {
	int j = 0;
    int b = 1;
    for (int i=0;i<nspin;i++) 
    {
        j+= b*arr[i];
        b = b*2;
    }
    return j;
}
// Conversion from number to state
void hsp::int_to_array( int *arr, int numb) {
	int j = numb;
    for(int i = 0;i<nspin;i++) {
        arr[i] = j%2;
        j = j/2;
    }
}

// Make Hilbert space 
void hsp::makehspace(){
    int i,j,k,l,ii,hdimhloc;
    //zeroing hspace
    for (i=0; i<hdim; i++) {
        hspace[i]=0;
    }
    // Determining starting point:
    for (i=0;i<nspin;i++) 
    {
        a[i] = (2*i<nspin)?1:0;
    }
    int snumb = array_to_int(a);
    // Selecting numbers
    j=0;
    for (i = snumb; i<hdim0+1;i++){
        int_to_array(a,i);
        k = 0;
        for (l=0;l<nspin;l++)
            k+=a[l];
        if (k==nspin2){
            hspace[j] = i;
            j+=1;
        }
    }
    //zeroing dict
    for (i=0; i<hdim0; i++) {
        dict[i]=0;
    }
    // Filling the dictionary
    for (i = 0; i<hdim;i++){
        dict[hspace[i]] = i;
    }
    //zeroing mult arrays
    for (i=0; i<hdimh; i++) {
        mult[i]=0;
    }
    for(j=0; j<nspin2; j++){
        for (i=0; i<hdim+1; i++) {
            multgen[j][i]=0;
        }
    }
    // Filling multiplicities array
    k = hspace[0]/hdimh;
    mult[0] = 1;
    ii = 0;
    for (i=1; i<hdim; i++) {
        j = hspace[i]/hdimh;
        if (j!=k) {
            ii = ii+1;
            k = j;
        }
        mult[ii] = mult[ii]+1;
    }
    // Filling generic multiplicities array for the state of rightmost ns-j spins in the beginning
    for (j=0; j<nspin2; j++){
        hdimhloc =  ((int)pow(2,j+1));
        k = hspace[0]/hdimhloc;
        multgen[j][0] = 1;
        ii = 0;
        for (i=1; i<hdim; i++) {
            l = hspace[i]/hdimhloc;
            if (l!=k) {
                ii = ii+1;
                k = l;
            }
            multgen[j][ii] = multgen[j][ii]+1;
        }
        multgen[j][ii+1] = -1;
    }
/*    int norm;
    for(j=0; j<nspin2; j++){
        cout <<"hspace for j="<<j<<"={";
        norm=0;
        for (i=0; i<hdim+1; i++) {
            cout<<multgen[j][i]<<";";
            norm+=multgen[j][i];
        }
        cout<<"}"<<endl;
        cout<<"norm="<<norm<<endl;
    }
 */
    
}

// Print Hspace
void hsp::printhspace(){
    int i,j;
    cout <<"Hspace\n";
    for (i=0;i<hdim;i++) 
    {
        cout << hspace[i]<< "; ";
    }
    cout <<"\n Dictionary \n";
    for (i=0;i<hdim0;i++) 
    {
        cout << dict[i]<< "; ";
    }
    cout <<"\n States in Hspace \n";
    for (j=0;j<hdim;j++)
    {
        int_to_array(a,hspace[j]);
        cout <<"num = "<<hspace[j] <<" array=";
        for (i=0;i<nspin;i++) 
        {
            cout << a[nspin-1-i]<< "; ";
        }
        cout <<"\n";
    }
    cout <<"\n Multiplicity \n";
    for (i=0;i<(int) pow(2,nspin2);i++) 
    {
        cout << mult[i]<< "; ";
    }
    cout <<"\n";
}

integer hsp::gethdim(){
    return hdim;
}

integer hsp::gethdimh(){
    return hdimh;
}

hsp::~hsp() {
    cout << "%Destructing hspace.. bye!\n";
    delete[] a;
    delete[] hspace;
    delete[] dict;
    delete[] mult;
    destroy(multgen, nspin2);
}




