#ifndef LATTICE_H
#define LATTICE_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <thread>
#include <iostream>
#include <unistd.h>
#include <valarray>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <stdarg.h>
#include <algorithm>

#include "atom.h"
#include "modelstats.h"
#include "MersenneTwister.h"

#define SQUARE(A) (A*A)
#define CUBE(A) (A*A*A)
#define MAXNNDIS 6.0
#define TOL 1e-3

class lattice
{
  public:
        MTRand RNDMNO;
        string Structure,DIRAniso;
        double c2aRatio, STRUCTCUBDIS, STRUCTTETDIS,Kaniscub,Kanistet;
        double J,K,Ut,Uc,U1,U2;
        double xsub,ysub,zsub;
        double **REPDATA,**STRUCTDATA,**MAGCUBICDATA,**MAGTETRADATA;
        int ORGCELLREPLICA, REPNOS, STRDATATYPE,EAATOMTOMS,STRUCTPOINTS,MAGCUBICPOINTS,MAGTETRAPOINTS,POLY;
        lattice();
        virtual ~lattice();
        void setatoms(vector <atom> atoms);
        void AtomSubstitute ( int ATin, int ATout, int COUNTREP );
        void CreateSimDomain (modelstats modelnum);
        void setSIMDOMAINATOMS(int SIMDOMAINATOMS);
        int  getSIMDOMAINATOMS();
        void setSIMDOMAINDATA(int i, atom thisatom);
        void CreateMagNNList (string NAME, int INITJ, double **DATAINITJ, double _c2aRatio);
        void MAXStrNNList (string NAME, double _c2aRatio);
        void DISStrNNList (string NAME, double _c2aRatio);
        void Sleep(float s)
      {
           int sec = int(s*1000000);
           usleep(sec);
      }
      void SwapOneDMatrix ( int i, int j, int SIZE, double **DATA )
{
	double *TEMPDATA;
	TEMPDATA = new double [SIZE];
	if ( !TEMPDATA )
	{
		cout<<"Error allocating memory to TEMPDATA in SwarOneDMatrix ! Exiting ... "<<endl;
		exit(0);
	}
	for ( int VAR = 0; VAR < SIZE; VAR++ )
	{
	    TEMPDATA[VAR] = DATA[i][VAR];
	}
	for ( int VAR = 0; VAR < SIZE; VAR++ )
	{
	    DATA[i][VAR] = DATA[j][VAR];
	}
	for ( int VAR = 0; VAR < SIZE; VAR++ )
	{
	    DATA[j][VAR] = TEMPDATA[VAR];
	}
}
       vector<atom> getatoms();
       atom* getSIMDOMAINDATA();
       double TotalEnergy( modelstats modelnum, double kBTLnp, double gmH);
       double delStrEng (int sigmahext, int i,int SIGMAi, double kBTLnp, double gmH );
       double DeformTrialEq (int sigmahext, int i, double kBTLnp, double gmH );
       double SpinTrial (modelstats modelnum, double kBTLnp,double gmH,int i );
       double SpinTrial_au(modelstats modelnum, double kBTLnp,double gmH,int i );
       double delMagEng_old(int spinhext, int i, int SPINi, int SPINn,double gmH);
       double delMagEng_old_au(int spinhext, int i, int SPINi, int SPINn,double gmH);
       double DeformTrialPd (int sigmahext,int spinhext,int i,double kBT,double kBTLnp,double gmH);
       double IntMagEng ( int spinhext,int i,  int SPINi, int SIGMAi ,double gmH);
       double delTotalEnergy (modelstats modelnum, double kBTLnp,double gmH,int i,int SPINi,int SPINn);
       double delTotalEnergy_au(modelstats modelnum, double kBTLnp,double gmH,int i,int SPINi,int SPINn);
       double NMCalc();
       double NMCalcnew();
       double TDCalc();
       void convertSIMDOMAINDATA();
       void creatsubdomain();
       void creatprob(double gmH);
       void assignkspin();
       double MagCalc(int isau);
       void assignspin();
       void assignsigma();
    protected:

    private:
        vector<atom> _atoms;
        int _SIMDOMAINATOMS;
        atom* _SIMDOMAINDATA;
        int *_SIMDOMAINATOMCOUNT;
        vector < vector <double > > _SIMADDRESS;
        vector < map <double, vector<int>> > _MCNNLIST, _MTNNLIST, _SCNNLIST, _STNNLIST;
        vector < int > _SPINLIST;
        vector < int > _SIGMALIST;
        vector < int > _KSPINLIST;
        vector < double > _CubmagmomentinSD;
        vector < double > _TetmagmomentinSD;
        map<int,vector<int > > _subdomainnnlist;//identify the neighbour domain of each domain
        vector<int > _subdomainnumber;//store domain of each atom
        vector<int > _tetproblist,_cubproblist;//whether atom can interacting with other domains

};

#endif // LATTICE_H
