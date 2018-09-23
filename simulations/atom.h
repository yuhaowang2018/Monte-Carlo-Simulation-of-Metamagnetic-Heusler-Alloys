#ifndef ATOM_H
#define ATOM_H

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
#include <string>
#include <map>
using namespace std;

class atom
{
    public:
        string name;
        int MAGSTATES,MAGCODES;
        double CubMagMoment,TetMagMoment;
        void setnumber(int number);
        int getnumber();
        void setcoordinate(double x, double y,double z);
        void SetMCNN(map <double, vector<int>> MCNN);
        void SetMTNN(map <double, vector<int>> MTNN);
        void SetSCNN(map <double, vector<int>> SCNN);
        void SetSTNN(map <double, vector<int>> STNN);
        map <double, vector<int>> getMCNN();
        map <double, vector<int>> getMTNN();
        map <double, vector<int>> getSCNN();
        map <double, vector<int>> getSTNN();
        double getx();
        double gety();
        double getz();
        void setcount(int count);
        int getcount();
        void setspin(int SPIN);
        int  getspin();
        void setsigma(int SIGMA);
        int  getsigma();
        void setkspin(int KSPIN);
        int  getkspin();
        void setmag(double mag){_mag=mag;}
        double getmag(){return _mag;}
        atom();
        virtual ~atom();


    protected:

    private:
         int _number;
         double _x,_y,_z;
         int _count;
         map <double, vector<int>> _MCNN, _MTNN, _SCNN, _STNN;
         int _SPIN;
         int _SIGMA;
         int _KSPIN;
         double _mag;



};

#endif // ATOM_H
