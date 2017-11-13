#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include<vector>
#include<sstream>
#include <string>
#include <fstream>
#include <complex>

using namespace std;

string trim(string);
// returns the actual index of a circular array
int loop_index(int i,int len);

void print_vector(vector<double> &vec);
void print_vector(vector<int> &vec);
/**/
class function_array
{
 private:
  vector<double> values;
  double a;// lower range
  double b;// higher range
  int bins;//number of bins in the array
 public:
  // loads an array from a file
  void load_file(string filename);
  double operator()(double value); // get the value of the function
  
};
void clear(vector<double> &v);
string int_to_string(int i);
string real_to_string(double r);
int vec_max_int(vector<int> &vec);
int vec_min_int(vector<int> &vec);
int minIndex(vector<double> &vec);
void print_matrix(vector< vector<double> > &matrix);
template<class T>
inline void zero_element(T &el)
{
  el=0;
}

template<class T>
inline void zero_element_vector(vector<T> &el,int n)
{
  el.resize(n);
  std::fill(el.begin(),el.end(),0);
}

template<class T>
string to_string(const T &e);


vector<double> operator+(const vector<double> &v1,const vector<double> &v2);
vector<double> operator/(const vector<double> &v1,const double & j);
// template<class T>
// void load_vector(const ifstream &f,vector<T> & v)
// {
//   int i;
//   for(i=0;i<v.size();i++)
//     {
//       f>>v[i];
//     }
  
// }
vector<double> operator*(const vector<double> &v1,const vector<double> &v2);
vector<double> operator-(const vector<double> &v1,const vector<double> &v2);

vector<double> operator/(const vector<double> &v1,const vector<double> &v2);


void real_part(double& v,double &p);
void real_part(complex<double>  &v,double &p);

double real_part(complex<double> & v);
double real_part(double & v);

// count the number of elements in the matrix
template<class T>
int get_counts(const vector< vector<T> > &v)
{
  int i,j,c;
  c=0;
  
  for(i=0;i<v.size();i++)
    {
      c+=v[i].size();
    }
  return c;
}
// read a matrix from a file
double string_to_double(const string & s);

vector< vector<double> > readMatrixFromFile(const string &filename);


// a vector of fixed length
template<class T,int N>
class tinyVector
{
public:
  tinyVector()
  {
    a[0]=0;
    a[1]=0;
  }
  
  T& operator[](int i)
  {
    return a[i];
  }

  const T& operator[](int i) const
  {
    return a[i];
  }
  
  tinyVector<T,N> & operator=(T o)
  {
    int i;
    for(i=0;i<N;i++)
      {
	a[i]=o;
      }
  }
  
private:
  
  T a[N];
};

#include "tools.hpp"

template<class T>
void cleanMatrix(T &M)
{
  int i,j;
  
  for(i=0;i<M.size();i++)
    {
      for(j=0;j<M[j].size();j++)
	{
	  M[i][j]=0;
	}
      
    }
}

void printColumnArray(double* A,int n,int m);

namespace tools
{
  void print(vector<double> &vector);
}

#endif
