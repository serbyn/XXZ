#ifdef __APPLE__
#include <math.h>
#else
#include "/usr/include/math.h"
#endif
#include <iostream>
#include <ran.h>
#include <hsp.h>
#include <ham.h>
#include <dmt.h>
#include <obs.h>
#include <obs_band.h>
#include <obs_frac.h>
#include <obs_mel.h>
#include <globalfunctions.h>
#include "time.h"
#include <sstream>
#include <fstream>

#include "H5Cpp.h"
#include "hdf5_hl.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;
Ran myran(1);

inline bool file_exists_test (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

// running measurments using ED data
void runmeas_ent (int nspin, int nruns, double hzinp, double Jzinp, string ver)
{   //Global intialization
    // Helper ints
    string measname;
    std::ostringstream out;
    out.str("");

    int i,k,hav, kmax;
    // Number of hz configurations
    hav = nruns;
    double time_start, time_end;
    // creating structures
    hsp sc_hsp(nspin);
    ham sc_ham(&sc_hsp);
    dmt sc_dmt(&sc_hsp);
    // Name of input file
    out.str("");
    out << ver<<"_"<<nspin<<"_"<<Jzinp<<"_"<<hzinp<<".h5";
    string fname;
    string gname;
    fname = out.str();

    // Name of output file
    out.str("");
    out << ver<<"_entspec_"<<nspin<<"_"<<Jzinp<<"_"<<hzinp<<".h5";
    string fname_spec;
    fname_spec = out.str();
    
    int kread[1];
    int par[1]={nspin};
    double par_d[1]={hzinp};
    IntType inttype( PredType::NATIVE_INT );
    inttype.setOrder( H5T_ORDER_LE );
    IntType doubtype( PredType::NATIVE_DOUBLE);
    doubtype.setOrder( H5T_ORDER_LE );
    hsize_t     dimsf[1];              // dataset dimensions
    dimsf[0] = 1;
    H5std_string  FILE_NAME( &fname[0] );
    H5std_string  FILE_NAME_SPEC( &fname_spec[0] );
    H5File *file,*file_spec;
    DataSet datasetN;
    DataSet dataset;
    if (file_exists_test(fname)){
        cout<<"%File  "<<fname<<"  exists ==> doing measurment using data from existing file"<<endl;
        file = new H5File( FILE_NAME, H5F_ACC_RDWR);
        //H5LTread_dataset(file,"/nruns",kread);
        datasetN  = file->openDataSet("/nruns");
        datasetN.read(kread, PredType::NATIVE_INT);
        kmax = kread[0]+1;
        cout<<"%File  "<<fname<<"  contains "<<kmax<<" realizations;"<<endl;
        if (kmax>nruns) kmax = nruns;
        cout<<"%Using "<<kmax<<" of them... beginning measurments..."<<endl;

        file_spec = new H5File( FILE_NAME_SPEC, H5F_ACC_TRUNC );
        DataSpace dataspace( 1, dimsf );
        H5std_string  Bpar_NAME( "L" );
        dataset = file_spec->createDataSet( Bpar_NAME, inttype, dataspace );
        dataset.write(par, PredType::NATIVE_INT );
        Bpar_NAME = "W";
        dataset = file_spec->createDataSet( Bpar_NAME, doubtype, dataspace );
        dataset.write(par_d, PredType::NATIVE_DOUBLE );
        Bpar_NAME = "Jz";
        dataset = file_spec->createDataSet( Bpar_NAME, doubtype, dataspace );
        par_d[0]=Jzinp;
        dataset.write(par_d, PredType::NATIVE_DOUBLE );
        Bpar_NAME = "nruns";
        datasetN = file_spec->createDataSet( Bpar_NAME, inttype, dataspace);
        par_d[0]=kmax;
        datasetN.write(par_d, PredType::NATIVE_DOUBLE );
    }else {
        cout<<"%FILE "<<fname<< " DOES NOT EXIST!!!, TERMINATING"<<endl;
        exit(-1);
    }
    Group group;
    H5std_string  Warr_NAME( "/Warr" );
    H5std_string  Jzz_NAME( "/Jzz" );
    H5std_string  Jperp_NAME( "/Jperp" );
    H5std_string  W_NAME( "/W" );
    H5std_string  A_NAME( "/EntSpec" );
    dimsf[0] = nspin;
    DataSpace dataspace1( 1, dimsf );
    dimsf[0] = sc_hsp.hdim;
    DataSpace dataspace2( 1, dimsf );
    dimsf[0] = sc_dmt.size*sc_ham.size;
    DataSpace dataspace3( 1, dimsf );
    
    // arrays with parameters of Hamiltonian
    double * hz;
    hz = new double[nspin];
    double * Jz;
    Jz = new double[nspin];
    double * Jp;
    Jp = new double[nspin];
    for (i=0;i<nspin;i++)
    {
        hz[i]= hzinp*(2.*(myran.doub())-1.);
        Jz[i] = Jzinp;
        Jp[i] = 1.;
    }
    sc_ham.setparam(hz,Jz,Jp);
    sc_ham.createH();
    
    double *specs;
    specs = new double[sc_dmt.size*sc_ham.size];
    // loop over different realizations of disorder
    for (k=0; k<kmax; k++) {
        out.str("");
        out << "/run_"<<k;
        gname = out.str();
        // reading from hdf5
        time_start = time(0);
        dataset  = file->openDataSet(gname+"/W");
        dataset.read(sc_ham.W, PredType::NATIVE_DOUBLE);
        dataset  = file->openDataSet(gname+"/A");
        dataset.read(sc_ham.A, PredType::NATIVE_DOUBLE);
        time_end = time(0);
        cout <<"%k="<<k<<"; Elapsed for reading  data time is t="<<time_end-time_start<<endl;
        time_start = time(0);
        group = file_spec->createGroup(gname);
        dataset = file_spec->createDataSet( gname+W_NAME, doubtype, dataspace2 );
        dataset.write( sc_ham.W, PredType::NATIVE_DOUBLE );
        for (int v=0; v<sc_ham.size; v++) {
            sc_dmt.constructrho(v, &sc_ham);
            sc_dmt.svd();
            memcpy(specs+v*sc_dmt.size, sc_dmt.W, (sizeof(double))*sc_dmt.size);
        }
        dataset = file_spec->createDataSet( gname+A_NAME, doubtype, dataspace3 );
        dataset.write( specs, PredType::NATIVE_DOUBLE );
        time_end = time(0);
        cout <<"Elapsed for measuring data time is t="<<time_end-time_start<<endl;
    }
    // final output

    
    
    delete[] specs;
    file->close();
    file_spec->close();
    return;
}


int main (int argc, char const *argv[]){ // hz, Jz, name of version
    //initialize random generator!
    unsigned long long int sek;
    sek = time(0);
    myran.RanI(sek);
    //cout<<"Generator is initialized with seed = "<<sek<<endl;
    double time_start, time_end;
    time_start = time(0);
    double hzinp, Jzinp;
    int L,N;
    string ver;
    // different modes depending on the number of arguments, function is called with:
    if (argc == 4){//running full simulation
        hzinp  = atof(argv[1]);
        Jzinp  = atof(argv[2]);
        ver = argv[3];
        //format for calling: runsim(nspins, nruns, hz, Jz, version)
        runmeas_ent(8,10,hzinp,Jzinp,ver);
        //runsim(16,15,hzinp,Jzinp,ver);
    } else if (argc == 6){//running full simulation
        hzinp  = atof(argv[3]);
        Jzinp  = atof(argv[4]);
        L = (int) atof(argv[1]);
        N = (int) atof(argv[2]);
        ver = argv[5];
        //format for calling: runsim(nspins, nruns, hz, Jz, version)
        cout<<"L="<<L<<"; nruns="<<N<<"; W="<<hzinp<<"; V="<<Jzinp <<endl;
        runmeas_ent(L,N,hzinp,Jzinp,ver);
    }
    else{//running in the test mode
        cout<<"%Test mode"<<endl;
        //runsim(12,5,.5,1.,"test");
        runmeas_ent(8,100,.6,1.,"h5test");
    }
    time_end = time(0);
    //cout<<"Total running time was"<<time_end-time_start<<endl;
    return 0;
}


