#ifndef MODELSTATS_H
#define MODELSTATS_H

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

//#include "lattice.h"
#include "atom.h"
using namespace std;


class modelstats
{
    public:
        string WSPINFiles,DomainFileFlag,DomainFileName;
        double Tin,Tout,delT,Hext;
        int EQUILIBRIUM,PRODUCTION,COLLECTION;
        modelstats();
        virtual ~modelstats();


//        double TotalEnergy(lattice latticeobj, double kBTLnp, double KgmH);
        int getspinhext(){
            return _SPINHext;
        }
        int getsigmahext(){
            return _SIGMAHext;
        }
    protected:

    private:
         int _SPINHext, _SIGMAHext;
};

#endif // MODELSTATS_H
