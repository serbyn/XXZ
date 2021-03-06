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
#include <globalfunctions.h>
#include "time.h"
#include <sstream>
#include <fstream>


using namespace std;
Ran myran(1);

void runsim (int nspin, int nruns, double hzinp, double Jzinp, string ver)
{   //Global intialization
    // Name of output file
    std::ostringstream out;
    out.str("");
    out << ver<<"_"<<nspin<<"_"<<Jzinp<<"_"<<hzinp<<".out";
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
    dmt sc_dmt(&sc_hsp);
    obs_band sc_obs_band(&sc_hsp, &sc_ham, &myran, nruns, fileout, &sc_dmt);
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
        runsim(L,N,hzinp,Jzinp,ver);
    }
    else{//running in the test mode
        cout<<"%Test mode"<<endl;
        //runsim(10,1,.4,1.,"test");
        runsim(12,5,.5,1.,"test");
        //runsim(8,10,0.5,.1,"test");
    }
    time_end = time(0);
    //cout<<"Total running time was"<<time_end-time_start<<endl;
    return 0;
}


