#ifndef RANDOM_H
#define RANDOM_H

#include<vector>
#include<cmath>
using namespace std;

//random number generator implementation
class ranlux
{
  static const int r=24;
  static const int s=10;
  static const int max_int=16777216;
  int c; // the carry bit
  int luxury;// luxury level(0-4), determines how good the generator
  bool skip;
  int n; // number of values generated
  int ni;//current value of the index
  vector<int> p_values; // p(luxury) table
  vector<int> last; // the last integers in the sequence
    
 public:
  ranlux(int seed);
  void rand(vector<double> & vec); //fill the vector with random numbers
  void seed(int seed,int lu); // uses a seed to implement the ranlux algorihm
};

class random1
{
  ranlux prng;
 public:
  random1(int seed);
  void uniform(vector<double> &vec);
  double uniform();
  void gaussian(vector<double> & vec);
  void gaussian(vector<double> & vec,size_t nmax);
};

void gaussian(vector<double> &vec);

#endif
