#include<vector>
#include<cmath>
#include"random.h"
#include"tools.h"
#include"cassert"
#include "ranlxd.h"
using namespace std;

void ranlux::seed(int seed,int lu=3)
{
  
  vector<int> B(31);
  int i=0,j=0,k=0,l=0,m=0;
  last.resize(24);
  // get the ith bith of the seed integer
  for(i=0;i<31;i++)
    {
      B[i]=(seed & ( 1 << i )) >> i;
    }
  j=0;
  k=14%31;
  // create the first 24 24-bit initial numbers
  
  for(m=0;m<24;m++)
    {
  // create 24 bits of the first initializing number
      for(l=0;l<24;l++)
	{
	  B[j]=(B[j]+B[k])%2;
	  j=(j+1)%31;
	  k=(k+1)%31;
	}
      // convert the last 24 bits to an integer and store it
      for (i=0;i<24;i++)
	{
	  
	  last[m]=last[m]+B[loop_index(j-i,31)]*pow(2,i);
	  
	}
      
    }
  // initial settings
  luxury=lu;
  c=0;
  // set the p-value table
  p_values.resize(5);
  p_values[0]=24;
  p_values[1]=48;
  p_values[2]=97;
  p_values[3]=223;
  p_values[4]=389;
  n=0;
  skip=false;
  
}
/*
void ranlux::rand(vector<double> & vec)
{
  // fill the vector with uniform random numbers
  int delta=0;
  int i=0;
  
  while(true)
    {
      
      ni=loop_index(ni+1,24);
      delta=last[loop_index(ni-s,24)]-last[loop_index(ni-r,24)]- c;
      // advance the loop counter
      
      if (delta > 0)
	{
	  last[ni]=delta;
	  c=0;
	}
      else
	{
	  last[ni]=delta + max_int;
	  c=1;
	}
      n=n+1;
      if (n > 24)
	{
	  skip=true;
	  if (n > p_values[luxury])
	    {
	      n=0;
	      skip=false;
	    }
	}
  
      if (!skip)
	{
	  
	  vec[i]=(double)last[ni]*1./max_int;
	  i=i+1;
	  if (i >=vec.size()) {break ;}
	}
      
    }
}
*/

void ranlux::rand( vector<double> &r)
{
  ranlxd(r);
}
ranlux::ranlux(int seed_int) 
{
  
  
  //seed(seed_int);
  rlxd_init(1,seed_int);
}

random1::random1(int seed_int) : prng(seed_int) 
{
  
}

void random1::uniform(vector<double> & vec)
{
  prng.rand(vec);
}

double random1::uniform()
{
  vector<double> tmp(1);
  prng.rand(tmp);
  return tmp[0];
}

void random1::gaussian(vector<double> & vec)
{
  gaussian(vec,vec.size());
}  

void random1::gaussian(vector<double> & vec,size_t nmax)
{
  double tmp1=0;
  double tmp2=0;
  unsigned int i=0;
  
  assert(nmax<= vec.size());
  // uniformly generated random numbers
  ranlxd(vec,nmax);
  
  for(i=0;i<nmax/2;i++)
    {
      tmp1=vec[2*i];
      tmp2=vec[2*i+1];
      // check if the numbers are legal
      assert( (tmp1 > 0) && (tmp1<=1));
      assert( (tmp2 > 0) && (tmp2 <= 1));
      // Box algorithm to generate a couple of random numbers
      vec[2*i]=sqrt(-2*log(tmp1))*cos(2*M_PI*tmp2);
      vec[2*i+1]=sqrt(-2*log(tmp1))*sin(2*M_PI*tmp2);
      
    }
  // add onother random number if the length of the vector is odd
  if (nmax%2 !=0)
    {
      
      tmp1=uniform();
      tmp2=uniform();
      
      vec[nmax-1]=sqrt(-2*log(tmp1))*cos(2*M_PI*tmp2);
      
    }
  
  
  
}

void gaussian(vector<double> &vec)
{
  
  double tmp1=0;
  double tmp2=0;
  vector<double> tmp(2);
  unsigned int i=0;
  
  // uniformly generated random numbers
  ranlxd(vec);
  
  for(i=0;i<vec.size()/2;i++)
    {
      tmp1=vec[2*i];
      tmp2=vec[2*i+1];
      // check if the numbers are legal
      assert( (tmp1 > 0) && (tmp1<=1));
      assert( (tmp2 > 0) && (tmp2 <= 1));
      // Box algorithm to generate a couple of random numbers
      vec[2*i]=sqrt(-2*log(tmp1))*cos(2*M_PI*tmp2);
      vec[2*i+1]=sqrt(-2*log(tmp1))*sin(2*M_PI*tmp2);
      
    }
  // add onother random number if the length of the vector is odd
  if (vec.size()%2 !=0)
    {
      
      ranlxd(tmp);
      
      vec[vec.size()-1]=sqrt(-2*log(tmp[0]))*cos(2*M_PI*tmp[1]);
      
    }
  
}
