#ifdef __APPLE__
#include <math.h>
#else
#include "/usr/include/math.h"
#endif
#include <iostream>
#include <ran.h>
#include <hsp.h>
#include <ham.h>
#include <obs.h>
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


// running just ED
void rundiag (int nspin, int nruns, double hzinp, double Jzinp, string ver)
{   //Global intialization
    // Helper ints
    int i,k,hav, kstart;
    // Number of hz configurations
    hav = nruns;
    double time_start, time_end;
    // creating structures
    hsp sc_hsp(nspin);
    ham sc_ham(&sc_hsp);
    // Name of output file
    std::ostringstream out;
    out.str("");
    out << ver<<"_"<<nspin<<"_"<<Jzinp<<"_"<<hzinp<<".h5";
    // pointer for output file
    string fname;
    string gname;
    fname = out.str();
    
    
    //fileout = fopen( &fname[0], "w");
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
    H5File *file;
    DataSet datasetN;
    DataSet dataset;
    if (file_exists_test(fname)){
        cout<<"File  "<<fname<<"  exists ==> appending to existing file"<<endl;
        file = new H5File( FILE_NAME, H5F_ACC_RDWR);
        //H5LTread_dataset(file,"/nruns",kread);
        datasetN  = file->openDataSet("/nruns");
        datasetN.read(kread, PredType::NATIVE_INT);
        kstart = kread[0]+1;
        cout<<"File  "<<fname<<"  contains "<<kstart<<" iterations"<<endl;
    }else {
        file = new H5File( FILE_NAME, H5F_ACC_TRUNC );
        kstart = 0;
        DataSpace dataspace( 1, dimsf );
        H5std_string  Bpar_NAME( "L" );
        dataset = file->createDataSet( Bpar_NAME, inttype, dataspace );
    
        dataset.write(par, PredType::NATIVE_INT );

        Bpar_NAME = "W";
        dataset = file->createDataSet( Bpar_NAME, doubtype, dataspace );
        dataset.write(par_d, PredType::NATIVE_DOUBLE );

        Bpar_NAME = "Jz";
        dataset = file->createDataSet( Bpar_NAME, doubtype, dataspace );
        par_d[0]=Jzinp;
        dataset.write(par_d, PredType::NATIVE_DOUBLE );

        Bpar_NAME = "nruns";
        datasetN = file->createDataSet( Bpar_NAME, inttype, dataspace);
    }
    Group group;
    H5std_string  Warr_NAME( "/Warr" );
    H5std_string  Jzz_NAME( "/Jzz" );
    H5std_string  Jperp_NAME( "/Jperp" );
    H5std_string  W_NAME( "/W" );
    H5std_string  A_NAME( "/A" );
    dimsf[0] = nspin;
    DataSpace dataspace1( 1, dimsf );
    dimsf[0] = sc_hsp.hdim;
    DataSpace dataspace2( 1, dimsf );
    dimsf[0] = sc_hsp.hdim2;
    DataSpace dataspace3( 1, dimsf );

    
    // arrays with parameters of Hamiltonian
    double * hz;
    hz = new double[nspin];
    double * Jz;
    Jz = new double[nspin];
    double * Jp;
    Jp = new double[nspin];
    bool skip=false;
    // loop over different realizations of disorder
    for (k=kstart; k<hav+kstart; k++) {
        par[0]=k;
        datasetN.write(par,PredType::NATIVE_INT);
        out.str("");
        out << "/run_"<<k;
        gname = out.str();
        if (skip) skip = false;
        else group = file->createGroup(gname);
        for (i=0;i<nspin;i++)
        {
            hz[i]= hzinp*(2.*(myran.doub())-1.);
            Jz[i] = Jzinp;
            Jp[i] = 1.;
        }
        // Open BC!!!
        //Jz[nspin-1] = 0.;
        //Jp[nspin-1] = 0.;
        // main Hamiltonian
        sc_ham.setparam(hz,Jz,Jp);
        sc_ham.createH();
        time_start = time(0);
        sc_ham.diagonalize();
        //sc_ham.printW();
        if (sc_ham.getinfo()>0) {
            sc_ham.resetinfo();
            k = k - 1;
            skip = true; // skipping the group creation!
            cout <<"skipping current iteration, due to error with sc_ham"<<endl;
            cout<<"printing values of parameters, hz="<<endl;
            for (i=0;i<nspin;i++)
            {
                cout << hz[i]<<", ";
            }
            cout<<"Jz[0]= "<<Jz[0]<<endl;
            cout<<"Jp[0]= "<<Jp[0]<<endl;
            continue;
        }
        time_end = time(0);
        cout <<"k="<<k<<"; Elapsed for diagonalization of H time is t="<<time_end-time_start<<endl;
        // writing to hdf5
        time_start = time(0);
        dataset = file->createDataSet( gname+Warr_NAME, doubtype, dataspace1 );
        dataset.write( hz, PredType::NATIVE_DOUBLE );
        dataset = file->createDataSet( gname+Jzz_NAME, doubtype, dataspace1 );
        dataset.write( Jz, PredType::NATIVE_DOUBLE );
        dataset = file->createDataSet( gname+Jperp_NAME, doubtype, dataspace1 );
        dataset.write( Jp, PredType::NATIVE_DOUBLE );
        dataset = file->createDataSet( gname+W_NAME, doubtype, dataspace2 );
        dataset.write( sc_ham.W, PredType::NATIVE_DOUBLE );
        dataset = file->createDataSet( gname+A_NAME, doubtype, dataspace3 );
        dataset.write( sc_ham.A, PredType::NATIVE_DOUBLE );
        time_end = time(0);
        cout <<"k="<<k<<"; Elapsed for writing  data time is t="<<time_end-time_start<<endl;

    }
    // final output
    delete[] hz;
    delete[] Jz;
    delete[] Jp;
    file->close();
    return;
}

// running measurments using ED data
void runmeas (int nspin, int nruns, double hzinp, double Jzinp, string ver)
{   //Global intialization
    // Helper ints
    string measname;
    std::ostringstream out;
    out.str("");
    out << ver<<"_mel_"<<nspin<<"_"<<Jzinp<<"_"<<hzinp<<".out";
    measname = out.str();
    FILE* fileout;
    fileout = fopen( &measname[0], "w");
    int i,k,hav, kmax;
    // Number of hz configurations
    hav = nruns;
    double time_start, time_end;
    // creating structures
    hsp sc_hsp(nspin);
    ham sc_ham(&sc_hsp);
    // Name of output file
    out.str("");
    out << ver<<"_"<<nspin<<"_"<<Jzinp<<"_"<<hzinp<<".h5";
    // pointer for output file
    string fname;
    string gname;
    fname = out.str();
    //fileout = fopen( &fname[0], "w");
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
    H5File *file;
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
    }else {
        cout<<"%FILE "<<fname<< " DOES NOT EXIST!!!, TERMINATING"<<endl;
        exit(-1);
    }
//    obs_band sc_obs_band(&sc_hsp, &sc_ham, &myran, kmax, fileout, &sc_dmt);
//    obs_frac sc_obs_frac(&sc_hsp, &sc_ham, &myran, kmax, fileout);
    obs_mel sc_obs_frac(&sc_hsp, &sc_ham, &myran, kmax, fileout);
    Group group;
    H5std_string  Warr_NAME( "/Warr" );
    H5std_string  Jzz_NAME( "/Jzz" );
    H5std_string  Jperp_NAME( "/Jperp" );
    H5std_string  W_NAME( "/W" );
    H5std_string  A_NAME( "/A" );
    dimsf[0] = nspin;
    DataSpace dataspace1( 1, dimsf );
    dimsf[0] = sc_hsp.hdim;
    DataSpace dataspace2( 1, dimsf );
    dimsf[0] = sc_hsp.hdim2;
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
        sc_obs_frac.measuredata();
        time_end = time(0);
        cout <<"Elapsed for measuring data time is t="<<time_end-time_start<<endl;
    }
    // final output
    sc_obs_frac.filep = fileout;
    sc_obs_frac.writebasicdata();
    sc_obs_frac.writedata();
    fclose(fileout);
    delete[] hz;
    delete[] Jz;
    delete[] Jp;
    file->close();
    return;
}


void runsim (int nspin, int nruns, double hzinp, double Jzinp, string ver)
{   //Global intialization
    // Name of output file
    std::ostringstream out;
    out.str("");
    out << ver<<"_mel_"<<nspin<<"_"<<Jzinp<<"_"<<hzinp<<".out";
    // pointer for output file
    FILE* fileout;
    string fname;
    fname = out.str();
    //fileout = fopen( &fname[0], "w");
    // Helper ints
    int i,k,hav;
    // Number of spins
    int info = 0;
    // Number of hz configurations
    hav = nruns;
    double time_start, time_end;
    // creating structures
    hsp sc_hsp(nspin);
    ham sc_ham(&sc_hsp);
    obs_mel sc_obs_band(&sc_hsp, &sc_ham, &myran, nruns, fileout);
    // arrays with parameters of Hamiltonian
    double * hz;
    hz = new double[nspin];
    double * Jz;
    Jz = new double[nspin];
    double * hz0;
    hz0 = new double[nspin];
    double * Jz0;
    Jz0 = new double[nspin];
    double * Jp;
    Jp = new double[nspin];
    // loop over different realizations of disorder
    for (k=0; k<hav; k++) {
        for (i=0;i<nspin;i++)
        {
            hz[i]= hzinp*(2.*(myran.doub())-1.);
            Jz[i] = Jzinp;
            hz0[i] = 0.001*hzinp*(2.*(myran.doub())-1.);
            Jz0[i] = 0;
            Jp[i] = 1.;
        }
        // Open BC!!!
        //Jz[nspin-1] = 0.;
        //Jp[nspin-1] = 0.;
        // main Hamiltonian
        sc_ham.setparam(hz,Jz,Jp);
        sc_ham.createH();
        time_start = time(0);
        sc_ham.diagonalize();
        //sc_ham.printW();
	if (sc_ham.getinfo()>0) {
            sc_ham.resetinfo();
            k = k - 1;
            cout <<"skipping current iteration, due to error with sc_ham"<<endl;
            cout<<"printing values of parameters, hz="<<endl;
            for (i=0;i<nspin;i++)
            {
                cout << hz[i]<<", ";
            }
            cout<<"Jz[0]= "<<Jz[0]<<endl;
            cout<<"Jp[0]= "<<Jp[0]<<endl;
            continue;
        }
        time_end = time(0);
        //cout <<"k="<<k<<"; Elapsed for diagonalization of H time is t="<<time_end-time_start<<endl;
        time_start = time(0);
        sc_obs_band.measuredata();
        if (((nspin==12)&&(k%50==0))||((nspin==14)&&(k%10==0))||(nspin>14)) {
            fileout = fopen( &fname[0], "w");
            sc_obs_band.filep = fileout;
            sc_obs_band.writebasicdata();
            sc_obs_band.writedata();
            fclose(fileout);
        }
        time_end = time(0);
        //cout <<"k="<<k<<"; Elapsed for measuring time is t="<<time_end-time_start<<endl;
    }
    // final output
    fileout = fopen( &fname[0], "w");
    sc_obs_band.filep = fileout;
    sc_obs_band.writebasicdata();
    sc_obs_band.writedata();
    fclose(fileout);
    
    delete[] hz;
    delete[] Jz;
    delete[] Jz0;
    delete[] hz0;
    delete[] Jp;
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
        runsim(8,30000,hzinp,Jzinp,ver);
        runsim(10,10000,hzinp,Jzinp,ver);
        runsim(12,4000,hzinp,Jzinp,ver);
        runsim(14,400,hzinp,Jzinp,ver);
        //runsim(16,15,hzinp,Jzinp,ver);
    } else if (argc == 6){//running full simulation
        hzinp  = atof(argv[3]);
        Jzinp  = atof(argv[4]);
        L = (int) atof(argv[1]);
        N = (int) atof(argv[2]);
        ver = argv[5];
        //format for calling: runsim(nspins, nruns, hz, Jz, version)
        cout<<"L="<<L<<"; nruns="<<N<<"; W="<<hzinp<<"; V="<<Jzinp <<endl;
        rundiag(L,N,hzinp,Jzinp,ver);
    }
    else{//running in the test mode
        cout<<"%Test mode"<<endl;
        //runsim(12,5,.5,1.,"test");
        runsim(8,100,.6,1.,"h5test");
    }
    time_end = time(0);
    //cout<<"Total running time was"<<time_end-time_start<<endl;
    return 0;
}


