#ifndef GRADIENTPARTICLES
#define GRADIENTPARTICLES

#include <vector>
#include<iostream>
#include "mpi.h"
using namespace std;
/*
  particle gradient object for 1D data
*/

class packed_data;

class allParticlesGradient1D
{
  
public:
  typedef double value_t;
  typedef vector<double> gradParticles_t;
  
  allParticlesGradient1D();
  allParticlesGradient1D(const vector<int> & Ns);
  void clone( allParticlesGradient1D & g2);
  void pack(packed_data* pack);
  void unpack(packed_data* pack);
  int getPackSize();
  vector<value_t>& operator[](size_t i){return gradient[i];}
  const vector<value_t>& operator[](size_t i) const {return gradient[i];}
  void resize(const vector<int> & Ns);
  void reset();
  void print();
  int size(){return gradient.size();}
  // add the second gradient and returns a reference to itself.
  // the object is modified
  
  allParticlesGradient1D& operator+=(const allParticlesGradient1D & grad2);
  
private:
  vector<vector<value_t> > gradient;
  
};

allParticlesGradient1D* buildAllParticlesGradient1D(string filename);

#endif
