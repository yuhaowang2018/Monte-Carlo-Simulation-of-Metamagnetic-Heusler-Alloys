#ifndef __DRAWSPINSTATES_H
#define __DRAWSPINSTATES_H

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <sstream>
#include <vector>
#include <stdarg.h>
#include <algorithm>

extern int *MAGSTATES, *SPIN, ORGCELLREPLICA, EAATOMTOMS, SIMDOMAINATOMS;
extern string *ATOMSNAME;
extern double **SIMDOMAINDATA, ZPlot;


ofstream wFile;
string Color[] = {"Red", "Blue", "Green", "Navy", "Maroon", "YellowGreen", "Cyan", "NavyBlue"};

string StrTemp(int number)
{
   stringstream ss;
   ss << number;
   return ss.str();
}

void SpinData (char* toFile, int TempR)
{
    string FileName, add;
    FileName = toFile;
    size_t pos;
    pos = FileName.find(".inp");
    add.append ("_T");
    add.append (StrTemp (TempR) );
    add.append (".str");
    FileName.replace (pos,4,add);
    wFile.open ( FileName.c_str());
    if ( !wFile.is_open () )
    {
        cout<<"ERROR : Opening Input File :: "<<FileName<<endl<<"\nExiting Now ... \n"<<endl;
        exit(0);
    }
    wFile << "title Spin_Monents"<< endl;
    wFile << "cell 20.0 20.0 20.0"<<endl;
    wFile << "spgp P" << endl;
    for ( int i = 0; i < EAATOMTOMS; i++ )
    {
        wFile << "sphere  "<< ATOMSNAME[i] << "  0.4  " << Color[i] <<" filter 0.8" <<endl;
    }
    int MAXMAGSTATES = MAGSTATES[0];
    for ( int i = 1; i < EAATOMTOMS; i++ )
    {
        if ( MAGSTATES[1] > MAXMAGSTATES )  MAXMAGSTATES = MAGSTATES[1];
    }
    double **ArrowVec, coorX, coorY, coorZ;
    ArrowVec = new double* [MAXMAGSTATES];
    for ( int i = 0; i < MAXMAGSTATES; i++ )
    {
        ArrowVec[i] = new double [2];
    }
    double AngIncr = 2*M_PI/MAXMAGSTATES;
    for ( int i = 0; i < MAXMAGSTATES; i++ )
    {
        ArrowVec[i][0] = cos(M_PI_2+i*AngIncr);
        ArrowVec[i][1] = sin(M_PI_2+i*AngIncr);
    }
    wFile<<setprecision(5);
    cout<<"\nZPlot :: "<<ZPlot<<endl;
    for ( int i = 0; i < SIMDOMAINATOMS; i++ )
    {
        coorX = SIMDOMAINDATA[i][1];
        coorY = SIMDOMAINDATA[i][2];
        coorZ = SIMDOMAINDATA[i][3];
        if ( coorZ < ZPlot )
        {
            coorX /= ORGCELLREPLICA;
            coorY /= ORGCELLREPLICA;
            coorZ /= ORGCELLREPLICA;
            if ( MAGSTATES[int(SIMDOMAINDATA[i][0])] != 0 )
            {
                wFile << "atom "<<ATOMSNAME[int(SIMDOMAINDATA[i][0])]<<"  "<<SPIN[i]<<"  "<<coorX<<"  "<<coorY<<"  "<<coorZ<<endl;
                wFile << "arrow "<<coorX<<" "<<coorY<<" "<<coorZ<<" "<<ArrowVec[SPIN[i]-1][0]<<" "<<ArrowVec[SPIN[i]-1][1]<<"  0.0000  "<< "2.0  0.2  Black"<<endl;
                if ( coorX < 1e-8 )
                {
                    coorX = 1.0;
                    wFile << "atom "<<ATOMSNAME[int(SIMDOMAINDATA[i][0])]<<"  "<<SPIN[i]<<"  "<<coorX<<"  "<<coorY<<"  "<<coorZ<<endl;
                    wFile << "arrow "<<coorX<<" "<<coorY<<" "<<coorZ<<" "<<ArrowVec[SPIN[i]-1][0]<<" "<<ArrowVec[SPIN[i]-1][1]<<"  0.0000  "<< "2.0  0.2  Black"<<endl;
                }
                else if ( coorY < 1e-8 )
                {
                    coorY = 1.0;
                    wFile << "atom "<<ATOMSNAME[int(SIMDOMAINDATA[i][0])]<<"  "<<SPIN[i]<<"  "<<coorX<<"  "<<coorY<<"  "<<coorZ<<endl;
                    wFile << "arrow "<<coorX<<" "<<coorY<<" "<<coorZ<<" "<<ArrowVec[SPIN[i]-1][0]<<" "<<ArrowVec[SPIN[i]-1][1]<<"  0.0000  "<< "2.0  0.2  Black"<<endl;
                }
            }
            else
            {
                wFile << "atom "<<ATOMSNAME[int(SIMDOMAINDATA[i][0])]<<"  1  "<<coorX<<"  "<<coorY<<"  "<<coorZ<<endl;
                if ( coorX < 1e-8 )
                {
                    coorX = 1.0;
                    wFile << "atom "<<ATOMSNAME[int(SIMDOMAINDATA[i][0])]<<"  1  "<<coorX<<"  "<<coorY<<"  "<<coorZ<<endl;
                }
                else if ( coorY < 1e-8 )
                {
                    coorY = 1.0;
                    wFile << "atom "<<ATOMSNAME[int(SIMDOMAINDATA[i][0])]<<"  1  "<<coorX<<"  "<<coorY<<"  "<<coorZ<<endl;
                }
            }
        }
    }
    wFile<<"nolabels \northographic \nview 0.0 0.0 0.0 \nend"<<endl;
    wFile.close();
}


#endif
