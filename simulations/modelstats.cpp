#include "modelstats.h"
#include "atom.h"
#include "lattice.h"
#include "MersenneTwister.h"

modelstats::modelstats()

{   Hext = 0;
    EQUILIBRIUM = 20000;
    PRODUCTION = 20000;
    COLLECTION = 100;
    WSPINFiles = "OFF";
    DomainFileFlag = "WRITE";
    DomainFileName = "DomainDataFile.sdd";
    _SPINHext=1;
    _SIGMAHext=1;
}

modelstats::~modelstats()
{

}







