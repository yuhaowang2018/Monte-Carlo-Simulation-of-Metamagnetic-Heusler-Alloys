#include "atom.h"

atom::atom()
{
    _number=-1;
    _x=-1;
    _y=-1;
    _z=-1;
    _count=0;
    _SPIN=1;
    _SIGMA=1;
    _KSPIN=1;
}

atom::~atom()
{
    //dtor
}


void atom::setnumber(int number)
{
  _number=number;
}

int atom::getnumber()
{
 return _number;
}

void atom::setcoordinate(double x,double y,double z)
{
  _x=x;
  _y=y;
  _z=z;
}

double atom::getx()
{
 return _x;
}

double atom::gety()
{
 return _y;

}

double atom::getz()
{
  return _z;
}

void atom::setcount(int count)
{
  _count=count;
}

int atom::getcount()
{
  return _count;
}

void atom::SetMCNN(map <double, vector<int>> MCNN){
 _MCNN=MCNN;
}

void atom::SetMTNN(map <double, vector<int>> MTNN){
   _MTNN=MTNN;
}

void atom::SetSCNN(map <double, vector<int>> SCNN){
   _SCNN=SCNN;
}

void atom::SetSTNN(map <double, vector<int>> STNN){
   _STNN=STNN;
}
map <double, vector<int>> atom::getMCNN()
{
  return _MCNN;
}

map <double, vector<int>> atom::getMTNN(){
  return _MTNN;
}

map <double, vector<int>> atom::getSCNN(){
   return _SCNN;
}

map <double, vector<int>> atom::getSTNN(){
   return _STNN;
}

void atom::setspin(int SPIN)
{
  _SPIN=SPIN;
}

int atom::getspin(){
    return _SPIN;
}

void atom::setsigma(int SIGMA)
{
  _SIGMA=SIGMA;

}

int atom::getsigma()
{
   return _SIGMA;
}

void atom::setkspin(int KSPIN)
{
    _KSPIN=KSPIN;
}

int atom::getkspin()
{
   return _KSPIN;
}
