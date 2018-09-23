#include <cstdio>
#include <cstdlib>
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
#include <sstream>

#include "atom.h"
#include "lattice.h"
#include "modelstats.h"
#include "ReadInput.h"
//#include "include/DrawSpinStates.h"

ofstream outFILE;
ifstream inFILE;

#define DISPLAYGRID
#define SQUARE(A) (A*A)

using namespace std;


//function to call display function for each lattice
//#ifdef DISPLAYGRID
//void outputLattices(vector<lattice> &models,FILE * outfile){
//  for(unsigned int i=0;i<models.size();i++)
//    models[i].display(outfile);
//}
//#endif /* DISPLAYGRID */

void OpenDataFiles ( int in_argc, char** in_argv )
{    string inFILEName, outFILEName;
    if ( in_argc > 1 )
    {
        inFILE.open ( in_argv[1], ios::in );
        if ( !inFILE.is_open () )   {
            cout<<"ERROR : Opening Input File :: "<<in_argv[1]<<endl<<"\nExiting Now ... \n"<<endl;
            exit(0);
        }
        inFILEName = in_argv[1];
		size_t pos;
		int isOUT = 0;
		for ( int i = 2; i < in_argc; i++ )    {
		    if ( ! ( strcmp (in_argv[i], "-o" ) ) ) {
		        isOUT = 1;
		        outFILEName = in_argv[i+1];
		    }
		}
		if ( isOUT == 0 )   {
            pos = inFILEName.find(".inp");
            inFILEName.replace (pos, 4, "");
            outFILEName = inFILEName + ".out";
            outFILE.open ( outFILEName.c_str(), ios::out );
            if ( !outFILE.is_open() )
            {
                cout<<"ERROR : Opening Output File :(\nExiting Now ... \n"<<endl;
                exit(0);
            }
		}
		else    {
		    outFILE.open ( outFILEName.c_str(), ios::out );
            if ( !outFILE.is_open() )   {
                cout<<"ERROR : Opening Output File :(\nExiting Now ... \n"<<endl;
                exit(0);
            }
		}
    }
    else
    {
        cout << "Program Usage ... \n";
        cout << in_argv[0] << "\tINPUT_File \t OPTIONS \nOptions Avialable... \n";
        cout<<"-temp \t\t Tin Tout delT\n";
        cout<<"-o \t\t OUTPUT_File\n";
        cout<<"-mcs \t\t Equilbrium_Steps Production_Steps DataCollection_Steps \n";
        cout<<"-h \t\t External Magnetic Field\n";
        cout<<"-sd \t\t WRITE - Write Simdomain file, READ - Read Simdomain file, Sim_Domain_FileName\n";
        cout<<"-cell \t\t Number of original cell replications\n";
        cout<<"-ieq \t\t Extra number of initial steps for equilibrium\n";
        exit(0);
    }
}




string StrData ( int NdS )
{
    stringstream ss;
    ss << NdS;
    return ss.str();
}


int main(int argc, char *argv[]){

    time_t TimeStart, TimeEnd;
    time(&TimeStart);

    lattice latticeobj;
    modelstats modelnum;
    OpenDataFiles ( argc, argv );
    
    ReadInputFile (latticeobj,modelnum);

    int initEQ = 0;
    latticeobj.CreateSimDomain(modelnum);



    time_t systemtime;
    time(&systemtime);
    cout<<"\nSYSTEM INITIALIZED ... PERFORMING CALCULATIONS ... "<<endl;
    outFILE<<"******************************************************\n";
    outFILE<<"*               SIMULATION DOMAIN & PARAMETERS\n";
    outFILE<<"*\t\t"<<ctime(&systemtime);
    outFILE<<"******************************************************\n";
    outFILE<<"STRCUTURAL PARAMETERS - J, K, U1, Uc, Ut :: "<<latticeobj.J<<"  "<<latticeobj.K<<"  "<<latticeobj.U1<<"  "<<latticeobj.Uc<<"  "<<latticeobj.Ut<<endl;
    outFILE<<"External Magnetic Field :: "<<modelnum.Hext<<endl;
    outFILE<<"c/a Ratio :: "<<latticeobj.c2aRatio<<endl;
    outFILE<<"EQUILBRIUM, PRODUCTION, DATA COLLECTION STEPS :: "<<modelnum.EQUILIBRIUM<<"\t"<<modelnum.PRODUCTION<<"\t"<<modelnum.COLLECTION<<endl;
    outFILE<<"TEMPERATURE OF SYSTEM :: "<<modelnum.Tin<<"\t"<<modelnum.Tout<<"\t"<<modelnum.delT<<endl;
    outFILE<<"******************************************************\n\n";

    latticeobj.CreateMagNNList(string("cubic MagNNlist"),latticeobj.MAGCUBICPOINTS,latticeobj.MAGCUBICDATA,1);
    cout<<"\nFinished Magnetic Cubic....\n";
    latticeobj.CreateMagNNList ("Magnetic Tetra", latticeobj.MAGTETRAPOINTS, latticeobj.MAGTETRADATA, latticeobj.c2aRatio);
    cout<<"Finished Magnetic Tetra...."<<endl<<endl;
    if ( latticeobj.STRDATATYPE == 0 )
    {
        latticeobj.MAXStrNNList ("Structural Cubic", 1);
        cout<<"\nFinished Structural Cubic - Based on Spehere Radius...."<<endl<<endl;
        latticeobj.MAXStrNNList ("Structural Tetra", latticeobj.c2aRatio);
        cout<<"\nFinished Structural Tetra - Based on Spehere Radius...."<<endl<<endl;
    }
    else if ( latticeobj.STRDATATYPE == 1 )
    {
        latticeobj.DISStrNNList ("Structural Cubic", 1);
        cout<<"\nFinished Structural Cubic - Based on Dis info...."<<endl<<endl;
        latticeobj.DISStrNNList ("Structural Tetra", latticeobj.c2aRatio);
        cout<<"\nFinished Structural Tetra - Based on Dis info...."<<endl<<endl;
    }
    else
    {
        cout<<endl<<"ERROR: Wrong value of variable STRDATATYPE :: "<< latticeobj.STRDATATYPE <<"  ... Please check your input file ... "<<endl;
        cout<<"\n\nExiting now .... "<<endl<<endl;
        exit(0);
    }

    int MAGATOMS = 0;
    for ( int k = 0; k < latticeobj.EAATOMTOMS; k++ )
    {
        if ( latticeobj.getatoms()[k].MAGSTATES>0 )
        {
            MAGATOMS += latticeobj.getatoms()[k].getcount();
        }
    }



 for ( int k = 0; k < latticeobj.getSIMDOMAINATOMS(); k++ )
     {

    if ( latticeobj.getSIMDOMAINDATA()[k].MAGSTATES==0 )
    {
                    latticeobj.getSIMDOMAINDATA()[k].setspin(0);
                    //latticeobj.getSIMDOMAINDATA()[k].setkspin(0);
    }
    }


    outFILE<<"T \t m \t m2 \t MSus\t eps \t eps2 \t SSus \t H \t H2 \t Cmag \t CmagT \t Entropy"<<endl;
    double m, m0, m2, eps, eps2, eps0, HEn0, HEn, HEn2, Cmag, CmagT, Sstart = 0, Entropy, Ssum = 0.0, MSus, SSus,delH;
    int COLLCOUNT = 0, PRCOUNTER = 1, EQCOUNTER = 0, CORRCount = 0, TCount = 0;
    double kB = 0.086173303, LandeFac = 2, muB = 0.057883818012, gmH,KgmH,kBT,kBTLnp,sigmahext,spinhext;
    gmH = LandeFac * muB * modelnum.Hext;

    latticeobj.convertSIMDOMAINDATA();
    sigmahext=modelnum.getsigmahext();
    spinhext=modelnum.getspinhext();
    latticeobj.creatsubdomain();
    latticeobj.assignkspin();
    latticeobj.creatprob(gmH);
    latticeobj.assignspin();
    
    int isau=0;
    double mag=latticeobj.MagCalc(isau);
    MTRand RNDMNO2;
    cout<<"MAGATOMS: "<<MAGATOMS<<endl;

   

       for ( int T = modelnum.Tin; T <= modelnum.Tout; T = T + modelnum.delT )
    {
            int flag=0;
            int i;
            vector<int> randomsitelist;
            vector<int>::iterator itt;
            while (flag==0)
            {  i=RNDMNO2.randInt(latticeobj.getSIMDOMAINATOMS()-1);
               itt=find(randomsitelist.begin(),randomsitelist.end(),i);
               if(itt==randomsitelist.end()||randomsitelist.size()==0)
               {
                randomsitelist.push_back(i);
               }
                if (randomsitelist.size()==latticeobj.getSIMDOMAINATOMS())
                {
                    flag=1;
                }
            }
            //latticeobj.creatprob(gmH);
        kBT = kB*T;
        kBTLnp = kBT*log(2);
        EQCOUNTER = 1, CORRCount = 0;

        cout<<"\nStarting Equilibrium for ... "<<T<<endl;

        if ( T == modelnum.Tin ) EQCOUNTER  = -initEQ;
        else            EQCOUNTER  = 0;
        while ( EQCOUNTER <= modelnum.EQUILIBRIUM )
        {
            for ( int ii = 0; ii < latticeobj.getSIMDOMAINATOMS(); ii++ )
             {   
                  int i=randomsitelist[ii];


               double delTrEng=latticeobj.DeformTrialEq(sigmahext,i,kBTLnp,gmH);

               if(isau==0) latticeobj.SpinTrial ( modelnum,  kBTLnp, gmH, i );
            else latticeobj.SpinTrial_au( modelnum,  kBTLnp, gmH, i );

            }

            EQCOUNTER++;
        }
        cout<<"Finished with equilibrium for ..."<<T<<endl;

        /***    Output spin data to a file - Atoms and spin directions ***/
       // if ( ! strcmp ( modelnum.WSPINFiles.c_str(), "ON" ) )  SpinData ( argv[1], T );

        m=0, m2=0, eps=0, eps2=0, HEn = 0, HEn2=0;
        COLLCOUNT = 0, PRCOUNTER = 1;
        HEn0 = latticeobj.TotalEnergy(modelnum,kBTLnp,gmH);
        vector <double> variance;
        double stand_dev=0;
        while ( PRCOUNTER <= modelnum.PRODUCTION)
        {
            for ( int ii = 0; ii < latticeobj.getSIMDOMAINATOMS(); ii++ )
            {    
                int i=randomsitelist[ii];
                HEn0 +=latticeobj.DeformTrialEq(sigmahext,i,kBTLnp,gmH);
             if(isau==0) HEn0+=latticeobj.SpinTrial ( modelnum,  kBTLnp, gmH, i );
             else  HEn0+=latticeobj.SpinTrial_au( modelnum,  kBTLnp, gmH, i );
             
            }
            if ( PRCOUNTER % modelnum.COLLECTION == 0 )
            {
                m0 = latticeobj.MagCalc(isau);
                variance.push_back(m0);
                eps0 = latticeobj.TDCalc ();
                m  += m0;
                m2 += SQUARE(m0);
                eps  += eps0;
                eps2 += SQUARE(eps0);
                HEn += HEn0;
                HEn2 += SQUARE(HEn0);
                COLLCOUNT++;
                CORRCount++;
            }
            PRCOUNTER++;
        }

          m  = m  / COLLCOUNT;
          for(int num=0; num<variance.size();num++)
          {
            stand_dev+=(variance[num]-m)*(variance[num]-m);
          }
          stand_dev/=COLLCOUNT;
          stand_dev=sqrt(stand_dev);
          m2 = m2 / COLLCOUNT;
          eps  = eps  / COLLCOUNT;
          eps2 = eps2 / COLLCOUNT;
          HEn  = HEn  / COLLCOUNT;

        if(abs(eps)<0.3&&isau!=1)
         {
            T=T-modelnum.delT;
            isau=1;
         }
         else{
         outFILE<<setprecision(12)<<T<<"\t"<<m<<"\t"<<stand_dev<<"\t"<<abs(eps)<<"\t"<<eps2<<endl;

         cout<<"\nFinished with production for ..."<<T<<endl;}

          TCount++;
    }

// cooling processs

       for ( int T = modelnum.Tout; T >= modelnum.Tin; T = T - modelnum.delT )
    {
            int flag=0;
            int i;
            vector<int> randomsitelist;
            vector<int>::iterator itt;
            while (flag==0)
            {  i=RNDMNO2.randInt(latticeobj.getSIMDOMAINATOMS()-1);
               itt=find(randomsitelist.begin(),randomsitelist.end(),i);
               if(itt==randomsitelist.end()||randomsitelist.size()==0)
               {
                randomsitelist.push_back(i);
               }
                if (randomsitelist.size()==latticeobj.getSIMDOMAINATOMS())
                {
                    flag=1;
                }
            }
           // latticeobj.creatprob(gmH);
        kBT = kB*T;
        kBTLnp = kBT*log(2);
        EQCOUNTER = 1, CORRCount = 0;

        cout<<"\nStarting Equilibrium for ... "<<T<<endl;

        if ( T == modelnum.Tout ) EQCOUNTER  = -initEQ;
        else            EQCOUNTER  = 0;
        while ( EQCOUNTER <= modelnum.EQUILIBRIUM )
        {
            for ( int ii = 0; ii < latticeobj.getSIMDOMAINATOMS(); ii++ )
             {   
                  int i=randomsitelist[ii];


                double delTrEng=latticeobj.DeformTrialEq(sigmahext,i,kBTLnp,gmH);
                if(isau==0) latticeobj.SpinTrial ( modelnum,  kBTLnp, gmH, i );
               else latticeobj.SpinTrial_au( modelnum,  kBTLnp, gmH, i );

            }
           
            EQCOUNTER++;
        }
        cout<<"Finished with equilibrium for ..."<<T<<endl;

        /***    Output spin data to a file - Atoms and spin directions ***/
       // if ( ! strcmp ( modelnum.WSPINFiles.c_str(), "ON" ) )  SpinData ( argv[1], T );

        m=0, m2=0, eps=0, eps2=0, HEn = 0, HEn2=0;
        COLLCOUNT = 0, PRCOUNTER = 1;
        HEn0 = latticeobj.TotalEnergy(modelnum,kBTLnp,gmH);
        while ( PRCOUNTER <= modelnum.PRODUCTION)
        {
            for ( int ii = 0; ii < latticeobj.getSIMDOMAINATOMS(); ii++ )
            {    //atom &thisatom=latticeobj.getSIMDOMAINDATA()[i];
                int i=randomsitelist[ii];
                HEn0 += latticeobj.DeformTrialEq(sigmahext,i,kBTLnp,gmH);
                if(isau==0) HEn0+=latticeobj.SpinTrial ( modelnum,  kBTLnp, gmH, i );
             else  HEn0+=latticeobj.SpinTrial_au( modelnum,  kBTLnp, gmH, i );
              
            }
            if ( PRCOUNTER % modelnum.COLLECTION == 0 )
            {
                m0 = latticeobj.MagCalc(isau);
                eps0 = latticeobj.TDCalc ();
                m  += m0;
                m2 += SQUARE(m0);
                eps  += eps0;
                eps2 += SQUARE(eps0);
                HEn += HEn0;
                HEn2 += SQUARE(HEn0);
                COLLCOUNT++;
                CORRCount++;
            }
            PRCOUNTER++;
        }

          m  = m  / COLLCOUNT;
          m2 = m2 / COLLCOUNT;
          eps  = eps  / COLLCOUNT;
          eps2 = eps2 / COLLCOUNT;
          HEn  = HEn  / COLLCOUNT;

         if(abs(eps)>0.3&&isau==1){
          T=T+modelnum.delT;
           isau=0;

         }
         else{
         outFILE<<setprecision(12)<<T<<"\t"<<m<<"\t"<<m2<<"\t"<<abs(eps)<<"\t"<<eps2<<endl;
         cout<<"\nFinished with production for ..."<<T<<endl;}

          TCount++;
    }

    outFILE<<"END"<<endl;
    outFILE.close();
    time(&TimeEnd);
    cout<<"Simulations Runtime: "<<int(difftime(TimeEnd, TimeStart) / 3600) <<" Hours, "<<(int(difftime(TimeEnd, TimeStart)) % 3600) / 60 <<" Minutes, "<<int(difftime(TimeEnd, TimeStart)) % 60<<" Seconds...\n"<<endl;
   return 0;
}


