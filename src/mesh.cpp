#include "mesh.h"
#include <iostream>
#include <cassert>
#include <cmath>
using namespace std;
mesh::mesh(double a_,double b_,int bins_)
{
  a=a_;
  b=b_;
  bins=bins_;
  step=(b-a)*1./bins;
}

int mesh::index(double x)
{
  
  int i;
  i=floor((x-a)/step);
  if (i==bins)
    {
      i=i-1;
    }
  
  assert(i>=0);
  assert(i<bins);
  return i;
}
