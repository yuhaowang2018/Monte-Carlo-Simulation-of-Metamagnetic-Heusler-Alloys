#ifndef __NEIGHBOURLIST_H
#define __NEIGHBOURLIST_H

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

#define SQUARE(A) (A*A)
#define CUBE(A) (A*A*A)
#define MAXNNDIS 6.0
#define TOL 1e-3

using namespace std;
extern ofstream outFILE;
extern double **SIMDOMAINDATA, STRUCTCUBDIS, STRUCTTETDIS, **STRUCTDATA;
extern int SIMDOMAINATOMS, ORGCELLREPLICA, STRUCTPOINTS;

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

void CreateMagNNList (string NAME, int INITJ, double **DATAINITJ, double c2aRatio, std::vector<vector<vector<double> > > &NNLIST)
{
    /****************************************
    MAGCUBICPOINTS -> INITJ - External
        Intial J Parameter between atoms
    TEMPMAGCUBICPTS -> FINALJ
        Final J points after adding symmetries
    MAGCUBICDATA -> DATAINITJ - External
        Data list of atoms, distance and J values
    TEMPMAGCUBICDATA -> DATAFINALJ
        Data list of atoms, distance and J values
        after adding symmetries
    SIMDOMAINATOMS -> SIMDOMAINATOMS - External
        Total atoms in the simualtion domain
    SIMDOMAINDATA -> SIMDOMAINDATA - External
        Atomic co-ordinate  data of all the atoms
        of simualtion domain
    MAGCUBCOUNT -> COUNTER
        Counter to generate final J data
    ORGCELLREPLICA -> ORGCELLREPLICA - External
        Repetations of the original cell in all directions
    MAGCUBICNNLIST -> NNLIST - External
        Neighbor list generated by the function
    ****************************************/

    /** Check magnetic data **/
    int ATOMi, ATOMj, FINALJ;
    FINALJ = INITJ;
    for ( int i = 0 ; i < INITJ; i++ )
    {
        ATOMi = int(DATAINITJ[i][0]);
        ATOMj = DATAINITJ[i][1];
        if ( ATOMi != ATOMj )
        {
            for ( int j = 0 ; j < INITJ; j++ )
            {
                if ( DATAINITJ[i][2] == DATAINITJ[j][2] )
                {
                    if ( ATOMj == DATAINITJ[j][0] && ATOMi == DATAINITJ[j][1] )
                    {
                        break;
                    }
                }
                if ( j == INITJ - 1 )
                {
                    FINALJ += 1;
                }
            }
        }
    }
    cout<<"MAGCUBICDATA"<<DATAINITJ<<endl;
    cout<<"List of Chains for "<<NAME<<" neighbor list : "<<FINALJ<<endl;

    double **DATAFINALJ;
    int *COUNTJPERAT;
    COUNTJPERAT = new int [EAATOMTOMS];
    if ( !COUNTJPERAT )
	{
		cout<<"Error allocating memory to COUNTJPERAT in CreateCubicNBList ! Exiting ... "<<endl;
		exit(0);
	}
    for ( int i = 0; i < EAATOMTOMS; i++ )
    {
        COUNTJPERAT[i] = 0;
    }
    DATAFINALJ = new double *[FINALJ];
    for ( int i = 0; i < FINALJ; i++ )
    {
        DATAFINALJ[i] = new double [4];
    }
    if ( !DATAFINALJ )
	{
		cout<<"Error allocating memory to DATAFINALJ in CreateCubicNBList ! Exiting ... "<<endl;
		exit(0);
	}

    /** Create new magnetic Cubic parameter lists **/
    int COUNTER = 0;
    for ( int i = 0 ; i < INITJ; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            DATAFINALJ[COUNTER][j] = DATAINITJ[i][j];
        }
        COUNTJPERAT[int(DATAINITJ[i][0])] += 1;
        COUNTER++;
    }
    for ( int i = 0 ; i < INITJ; i++ )
    {
        ATOMi = DATAINITJ[i][0];
        ATOMj = DATAINITJ[i][1];
        if ( ATOMi != ATOMj )
        {
            for ( int j = 0 ; j < INITJ; j++ )
            {
                if ( DATAINITJ[i][2] == DATAINITJ[j][2] )
                {
                    if ( ATOMj == DATAINITJ[j][0] && ATOMi == DATAINITJ[j][1] )
                    {
                        break;
                    }
                }
                if ( j == INITJ - 1 )
                {
                    DATAFINALJ[COUNTER][0] = ATOMj;
                    DATAFINALJ[COUNTER][1] = ATOMi;
                    DATAFINALJ[COUNTER][2] = DATAINITJ[i][2];
                    DATAFINALJ[COUNTER][3] = DATAINITJ[i][3];
                    COUNTJPERAT[ATOMj] += 1;
                    COUNTER += 1;
                }
            }
        }
    }

    cout <<NAME << " J List Before Sorting -> \n";
    for ( int i = 0; i < FINALJ; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            cout << DATAFINALJ[i][j]<<"\t";
        }
        cout<<endl;
    }

    /** Sort the cubic exchange list by
    First -> ATOMi & Second -> Distance **/
    int ATOMi1, ATOMi2, ATOMj1, ATOMj2;
    double DISi, DISj;
    for ( int i = 0; i < FINALJ; i++ )
    {
        for ( int j = i; j < FINALJ; j++ )
        {
            ATOMi1 = DATAFINALJ[i][0];
            ATOMi2 = DATAFINALJ[i][1];
            DISi = DATAFINALJ[i][2];
            ATOMj1 = DATAFINALJ[j][0];
            ATOMj2 = DATAFINALJ[j][1];
            DISj = DATAFINALJ[j][2];
            if ( ATOMj1 < ATOMi1 )
            {
                SwapOneDMatrix ( i, j, 4, DATAFINALJ );
            }
            else if ( ATOMi1 == ATOMj1 )
            {
                if ( ATOMj2 < ATOMi2)
                {
	            SwapOneDMatrix ( i, j, 4, DATAFINALJ );
                }
                else if ( ATOMi2 == ATOMj2 )
                {
                    if ( DISj < DISi )
                    {
                        SwapOneDMatrix ( i, j, 4, DATAFINALJ );
                    }
                }
            }
        }
    }

    cout <<"\n"<<NAME<< " J Exchange List after sorting-> \n";
    for ( int i = 0; i < FINALJ; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            cout << DATAFINALJ[i][j]<<"\t";
        }
        cout<<endl;
    }
    cout<< " Count J for each Atom : ";
    for ( int i = 0; i < EAATOMTOMS; i++ )
    {
        cout << COUNTJPERAT[i] << "\t";
    }
    cout<<"\n";

    /*** Create Cubic Neighbor Lists **/
    cout << "\nCreating "<<NAME<< " Neighbour List ... ";
    cout<<"\n";
    int JCOUNT, JNO;
    double XCOOD, YCOOD, ZCOOD, DISij, ERR, CX, CY, CZ, XHALF, YHALF, ZHALF;
    XHALF = ORGCELLREPLICA/2;
    YHALF = XHALF;
    ZHALF = XHALF*c2aRatio;
    NNLIST.resize(SIMDOMAINATOMS);
    for ( int i = 0; i < SIMDOMAINATOMS; i++ )
    {
        ATOMi1 = int(SIMDOMAINDATA[i][0]);
        JCOUNT = COUNTJPERAT[ATOMi1];
        if ( JCOUNT > 0 )
        {
            NNLIST[i].resize(JCOUNT);
            /** Code to add J parameter at start **/
            for ( int VAR = 0; VAR < FINALJ; VAR++ )
            {
                if ( ATOMi1 == DATAFINALJ[VAR][0] )
                {
                    for ( int k = 0; k < JCOUNT; k++ )
                    {
                        NNLIST[i][k].push_back(DATAFINALJ[VAR+k][3]);
                    }
                    break;
                }
            }
            for ( int ZTRANS = 0; ZTRANS < 2; ZTRANS++ )
            {
                for ( int YTRANS = 0; YTRANS < 2; YTRANS++ )
                {
                    for ( int XTRANS = 0; XTRANS < 2; XTRANS++ )
                    {
                        XCOOD = SIMDOMAINDATA[i][1];
                        YCOOD = SIMDOMAINDATA[i][2];
                        ZCOOD = SIMDOMAINDATA[i][3]*c2aRatio;
                        if ( XCOOD < XHALF )    XCOOD += ORGCELLREPLICA*XTRANS;
                        else                    XCOOD -= ORGCELLREPLICA*XTRANS;
                        if ( YCOOD < YHALF )    YCOOD += ORGCELLREPLICA*YTRANS;
                        else                    YCOOD -= ORGCELLREPLICA*YTRANS;
                        if ( ZCOOD < ZHALF )    ZCOOD += ORGCELLREPLICA*ZTRANS*c2aRatio;
                        else                    ZCOOD -= ORGCELLREPLICA*ZTRANS*c2aRatio;
                        for ( int j = 0; j < SIMDOMAINATOMS; j++ )
                        {
                            if ( j != i )
                            {
                                ATOMj1 = int(SIMDOMAINDATA[j][0]);
                                //cout<<"ATOMj1 "<<ATOMj1<<endl;
                                CX = XCOOD-SIMDOMAINDATA[j][1];
                                CY = YCOOD-SIMDOMAINDATA[j][2];
                                CZ = ZCOOD-SIMDOMAINDATA[j][3]*c2aRatio;
                                DISij = sqrt(SQUARE(CX)+SQUARE(CY)+SQUARE(CZ));

                                if ( DISij < MAXNNDIS )
                                {   //cout<<"DISij "<<DISij<<" ";
                                    JNO = -1;
                                    //cout<<"FINALJ "<<FINALJ<<endl;
                                    for ( int VAR = 0;  VAR < FINALJ; VAR++ )
                                    {
                                        ATOMi2 = int(DATAFINALJ[VAR][0]);
                                        if ( ATOMi1 == ATOMi2 )
                                        {
                                            JNO += 1;
                                            ATOMj2 = int(DATAFINALJ[VAR][1]);
                                            //cout<<"JNO"<<JNO<<endl;
                                            if ( ATOMj1 == ATOMj2 )
                                            {  //cout<<"ATOMj1 "<<ATOMj1<<endl;

                                                ERR = abs(DISij-DATAFINALJ[VAR][2]);

                                                if ( ERR < TOL )
                                                {
                                                    NNLIST[i][JNO].push_back(j);
                                                    //cout<<JNO<<endl;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    int NODES, NoCount = 0;
    for ( int i = 0; i < EAATOMTOMS; i++ )
    {
        NoCount = 0;
        for ( int j = 0; j < SIMDOMAINATOMS; j++ )
        {
            if ( SIMDOMAINDATA[j][0] == i )
            {
                NODES = NNLIST[j].size();
                cout<<"For Atom \""<<ATOMSNAME[i]<<"\" Neighbours number are ";
                if ( NODES == 0 )   cout<<" NONE ";
                for ( int k = 0; k < NODES; k++ )
                {
                    cout<<(NNLIST[j][k].size()-1)<<"\t";
                }
                cout<<endl;
                NoCount++;
                if (NoCount > 1 ) break;
            }
        }
    }

    for ( int i = 0; i < FINALJ; i++ )
    {
        delete [] DATAFINALJ[i];
    }
    delete [] DATAFINALJ;
    delete [] COUNTJPERAT;
}


void MAXStrNNList (string NAME, double c2aRatio, std::vector<vector<vector<double> > > &NNLIST)
{
    /***********************************************************
    SIMDOMAINATOMS -> SIMDOMAINATOMS - External
        Total atoms in the simualtion domain
    SIMDOMAINDATA -> SIMDOMAINDATA - External
        Atomic co-ordinate  data of all the atoms
        of simualtion domain
    MAGCUBCOUNT -> COUNTER
        Counter to generate final J data
    ORGCELLREPLICA -> ORGCELLREPLICA - External
        Repetations of the original cell in all directions
    MAGCUBICNNLIST -> NNLIST - External
        Neighbor list generated by the function
    ***********************************************************/

    double XCOOD, YCOOD, ZCOOD, DISij, CX, CY, CZ, XHALF, YHALF, ZHALF, ATOMICDIS;
    XHALF = ORGCELLREPLICA/2;
    YHALF = XHALF;
    ZHALF = XHALF*c2aRatio;
    NNLIST.resize(SIMDOMAINATOMS);
    if (c2aRatio == 1.0 ) ATOMICDIS = STRUCTCUBDIS;
    else                ATOMICDIS = STRUCTTETDIS;
    cout<<"The ATOMDIS to look for : "<<ATOMICDIS<<endl;
    for ( int i = 0; i < SIMDOMAINATOMS; i++ )
    {
        NNLIST[i].resize(1);
        for ( int ZTRANS = 0; ZTRANS < 2; ZTRANS++ )
        {
            for ( int YTRANS = 0; YTRANS < 2; YTRANS++ )
            {
                for ( int XTRANS = 0; XTRANS < 2; XTRANS++ )
                {
                    XCOOD = SIMDOMAINDATA[i][1];
                    YCOOD = SIMDOMAINDATA[i][2];
                    ZCOOD = SIMDOMAINDATA[i][3]*c2aRatio;
                    if ( XCOOD < XHALF )    XCOOD += ORGCELLREPLICA*XTRANS;
                    else                    XCOOD -= ORGCELLREPLICA*XTRANS;
                    if ( YCOOD < YHALF )    YCOOD += ORGCELLREPLICA*YTRANS;
                    else                    YCOOD -= ORGCELLREPLICA*YTRANS;
                    if ( ZCOOD < ZHALF )    ZCOOD += ORGCELLREPLICA*ZTRANS*c2aRatio;
                    else                    ZCOOD -= ORGCELLREPLICA*ZTRANS*c2aRatio;
                    for ( int j = 0; j < SIMDOMAINATOMS; j++ )
                    {
                        if ( j != i )
                        {
                            CX = XCOOD-SIMDOMAINDATA[j][1];
                            CY = YCOOD-SIMDOMAINDATA[j][2];
                            CZ = ZCOOD-SIMDOMAINDATA[j][3]*c2aRatio;
                            DISij = sqrt(SQUARE(CX)+SQUARE(CY)+SQUARE(CZ));
                            if ( DISij <= ATOMICDIS )
                            {
                                NNLIST[i][0].push_back(j);
                            }
                        }
                    }
                }
            }
        }
    }
    int NODES, NoCount = 0;
    for ( int i = 0; i < EAATOMTOMS; i++ )
    {
        NoCount = 0;
        for ( int j = 0; j < SIMDOMAINATOMS; j++ )
        {
            if ( SIMDOMAINDATA[j][0] == i )
            {
                NODES = NNLIST[j].size();
                cout<<"For Atom \""<<ATOMSNAME[i]<<"\" Neighbours number are ";
                if ( NODES == 0 )   cout<<" NONE ";
                for ( int k = 0; k < NODES; k++ )
                {
                    cout<<(NNLIST[j][k].size()-1)<<"\t";
                }
                cout<<endl;
                NoCount++;
                if (NoCount > 1 ) break;
            }
        }
    }

}

void DISStrNNList (string NAME, double c2aRatio, std::vector<vector<vector<double> > > &NNLIST)
{
    cout << "\nCreating "<<NAME<< " Neighbour List ... ";
    cout<<"\n";
    int ATOMi1, ATOMj1, ATOMi2, ATOMj2, DISCOL;
    double XCOOD, YCOOD, ZCOOD, DISij, ERR, CX, CY, CZ, XHALF, YHALF, ZHALF;
    XHALF = ORGCELLREPLICA/2;
    YHALF = XHALF;
    ZHALF = XHALF*c2aRatio;
    if ( c2aRatio == 1.0 )    DISCOL = 2;
    else                      DISCOL = 3;
    NNLIST.resize(SIMDOMAINATOMS);
    for ( int i = 0; i < SIMDOMAINATOMS; i++ )
    {
        ATOMi1 = int(SIMDOMAINDATA[i][0]);
        NNLIST[i].resize(1);
        for ( int ZTRANS = 0; ZTRANS < 2; ZTRANS++ )
        {
            for ( int YTRANS = 0; YTRANS < 2; YTRANS++ )
            {
                for ( int XTRANS = 0; XTRANS < 2; XTRANS++ )
                {
                    XCOOD = SIMDOMAINDATA[i][1];
                    YCOOD = SIMDOMAINDATA[i][2];
                    ZCOOD = SIMDOMAINDATA[i][3]*c2aRatio;
                    if ( XCOOD < XHALF )    XCOOD += ORGCELLREPLICA*XTRANS;
                    else                    XCOOD -= ORGCELLREPLICA*XTRANS;
                    if ( YCOOD < YHALF )    YCOOD += ORGCELLREPLICA*YTRANS;
                    else                    YCOOD -= ORGCELLREPLICA*YTRANS;
                    if ( ZCOOD < ZHALF )    ZCOOD += ORGCELLREPLICA*ZTRANS*c2aRatio;
                    else                    ZCOOD -= ORGCELLREPLICA*ZTRANS*c2aRatio;
                    for ( int j = 0; j < SIMDOMAINATOMS; j++ )
                    {
                        if ( j != i )
                        {
                            ATOMj1 = int(SIMDOMAINDATA[j][0]);
                            CX = XCOOD-SIMDOMAINDATA[j][1];
                            CY = YCOOD-SIMDOMAINDATA[j][2];
                            CZ = ZCOOD-SIMDOMAINDATA[j][3]*c2aRatio;
                            DISij = sqrt(SQUARE(CX)+SQUARE(CY)+SQUARE(CZ));
                            if ( DISij < MAXNNDIS )
                            {
                                for ( int VAR = 0;  VAR < STRUCTPOINTS; VAR++ )
                                {
                                    ATOMi2 = int(STRUCTDATA[VAR][0]);
                                    if ( ATOMi1 == ATOMi2 )
                                    {
                                        ATOMj2 = int(STRUCTDATA[VAR][1]);
                                        if ( ATOMj1 == ATOMj2 )
                                        {
                                            ERR = abs(DISij-STRUCTDATA[VAR][DISCOL]);
                                            if ( ERR < TOL )
                                            {
                                                NNLIST[i][0].push_back(j);
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    int NODES, NoCount = 0;
    for ( int i = 0; i < EAATOMTOMS; i++ )
    {
        NoCount = 0;
        for ( int j = 0; j < SIMDOMAINATOMS; j++ )
        {
            if ( SIMDOMAINDATA[j][0] == i )
            {
                NODES = NNLIST[j].size();
                cout<<"For Atom \""<<ATOMSNAME[i]<<"\" Neighbours number are ";
                if ( NODES == 0 )   cout<<" NONE ";
                for ( int k = 0; k < NODES; k++ )
                {
                    cout<<(NNLIST[j][k].size()-1)<<"\t";
                }
                cout<<endl;
                NoCount++;
                if (NoCount > 1 ) break;
            }
        }
    }

}

#endif