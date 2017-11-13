#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <complex>
#include <fstream>
#include "tools.h"

using namespace std;
// trim a strip
string trim(string str){
  int pos=0;
  // just remove white spaces at both ends
  pos=str.find_last_not_of(" \t\n");
  
  if( string::npos !=  pos )
    {
      str = str.substr( 0, pos+1 );
    }
  
  pos=str.find_first_not_of(" \t\n");
  if (string::npos != pos)
    {
      str=str.substr(pos,str.length());
    }
  return str;
}

int loop_index(int i/*plain index*/,int len/*length of the array*/)
{
  i=i%len;
  if (i<0) {i=i+len;}
  return i;
}

void print_vector(vector<double> &vec)
{
  int i=0;
  for (i=0;i<vec.size();i++)
    {
      cout<<vec[i]<<", ";
    }
  cout<<endl;
}

void print_vector(vector<int> &vec)
{
  int i=0;
  for (i=0;i<vec.size();i++)
    {
      cout<<vec[i]<<", ";
    }
  cout<<endl;
}
  
/*
------------discretized array of a function
*/

// convert a string to integer
int string_to_int(string name)
{
  return atoi(name.c_str());
}

// converts an integer to a string
string int_to_string(int i)
{
  stringstream s;
  s << i;
  return s.str();
}

double string_to_double(const string & s)
{
  return atof(s.c_str());
}
string real_to_string(double r)
{
  stringstream s;
  s << r;
  return s.str();
}

int vec_max_int(vector <int> &vec)
{
  unsigned int i;
  int max;
  int i_max=0;
  max=vec[0];
  for(i=0;i<vec.size();i++)
    {
      if (vec[i] > max)
	{
	  i_max=i;
	  max=vec[i];
	}
    }
  return i_max;
}

int vec_min_int(vector <int> &vec)
{
  unsigned int i;
  int min;
  int i_min=0;
  min=vec[0];
  for(i=0;i<vec.size();i++)
    {
      if (vec[i] < min)
	{
	  i_min=i;
	  min=vec[i];
	}
      
    }
  return i_min;
}

template<class T>
string to_string(const T &e)
  {
  stringstream s;
  s << e;
  return s.str();
};

template<>
string to_string<vector<double> >(const vector<double> &v)
{
  int i;
  string s="";
  for(i=0;i<v.size();i++)
    {
      s=s + string(" ") + to_string(v[i]);
    }
  return s;
};

template<>
string to_string<vector<int> >(const vector<int> &v)
{
  int i;
  string s="";
  for(i=0;i<v.size();i++)
    {
      s=s + string(" ") + to_string(v[i]);
    }
  return s;
};



vector<double> operator+(const vector<double> &v1,const vector<double> &v2)
{
  int i;
  vector<double> v=v1;
  
  for(i=0;i<v1.size();i++)
    {
      v[i]+=v2[i];
    }
  return v;
}

vector<double> operator/(const vector<double> &v1,const double &j)
{
  int i;
  vector<double> v=v1;
  for(i=0;i<v1.size();i++)
    {
      v[i]/=j;
    }
  return v;
}
vector<double> operator*(const vector<double> &v1,const vector<double> &v2)
{
  int i;
  vector<double> v=v1;
  for(i=0;i<v1.size();i++)
    {
      v[i]*=v2[i];
    }
  return v;
}
vector<double> operator/(const vector<double> &v1,const vector<double> &v2)
{
  int i;
  vector<double> v=v1;
  for(i=0;i<v1.size();i++)
    {
      v[i]/=v2[i];
    }
  return v;
}

vector<double> operator-(const vector<double> &v1,const vector<double> &v2)
{
  int i;
  vector<double> v=v1;
  for(i=0;i<v1.size();i++)
    {
      v[i]-=v2[i];
    }
  return v;
}

inline int get_size(vector<double> &v)
{
  return v.size();
}

inline int get_size(double &v)
{
  return 1;
}

double real_part(complex<double> & v)
 {
   return real(v);
 }
inline double imaginary_part(complex<double> &v)
{
  return imag(v);
}
double real_part(double & v)
 {
   return v;
 }

void real_part(complex<double>  &v,double &p)
{
  p=real(v);
};

void real_part(double  &v,double &p)
{
  
};

vector<vector <double> > readMatrixFromFile(const string & inputString)
{
  stringstream matrixStream;
  stringstream lineStream;
  string line;
  string value;
  vector<double> tmpVec;
  vector<vector <double> > M;
  int m;
  ifstream f;
  
  f.open(inputString.c_str());
  
  while(getline(f,line) )
    {
      lineStream.str("");
      lineStream.clear();
      lineStream<<line;
      
     while(getline(lineStream,value,' ') )
   	{
   	  tmpVec.push_back(string_to_double(value));
  	}
     }
  
   // sets equal to the size of the temporary vector
   m=tmpVec.size();
   M.resize(m);
   
   f.clear();
   f.seekg(0);
   
   int i=0;
   while(getline(f,line) )
     {
       lineStream.str("");
       lineStream.clear();
       lineStream<<line;
       i=0;
       while(getline(lineStream,value,' ') )
	 {
	   if (i >= m)
	     {
	       cout << "Bad matrix."<<endl;
	       exit(1);
	     }
	   
	   M[i].push_back(string_to_double(value));
	   i++;
  	}
     }
   
   // return a matrix of all value types
   return M;
}


void print_matrix(vector< vector<double> > &matrix)
{
  int i,j;
  
  for(i=0;i<matrix.size();i++)
    {
      for(j=0;j<matrix[i].size();j++)
	{
	  cout << matrix[i][j]<<" ";
	}
      cout << endl;
    }
}

void printColumnArray(double* A,int n,int m)
{
  for(int i=0;i<n;i++)
    {
      for(int j=0;j<m;j++)
	{
	  cout << A[j*n+i]<<" ";
	}
      cout << endl;
    }
}


// retruns the index of the minimum value of the array
int minIndex(vector <double> &vec)
{
  int iMin;
  iMin=0;
  
  for(int i=1;i<vec.size();i++)
    {
      if (vec[i]<vec[iMin])
	{
	  iMin=i;
	}
    }
  return iMin;
}

void tools::print(vector<double> &vec)
{
  int i=0;
  for (i=0;i<vec.size()-1;i++)
    {
      printf("%f ,",vec[i]);
    }
  printf("%f\n",vec[i]);
}
