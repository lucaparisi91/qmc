#ifndef OPTIMIZETOOLS_H
#define OPTIMIZETOOLS_H

#include <vector>
#include <cstdio>

using namespace std;
template<class T>
class rawMatrix
{
public:
  rawMatrix();
  void resize(unsigned int n,unsigned int m);
  const T& operator()(int i,int j) const;
  T& operator()(int i,int j);
  T* raw(){return &M[0];}
  int getN(){return n;} const
  int getM(){return m;} const
  void print() const;
  
private:
  inline unsigned int index(int i,int j) const{return n*j+i;} 
  unsigned int n,m;
  vector<T> M;
 
};

class linearMethodStepEstimator
{
public:
  linearMethodStepEstimator(){nP=0;};
  void buildMatrix(const vector<double> &accumulatedData,int nP);
  int getStep(vector<double>& parameters);
  void print() const {printf("H\n");H.print();printf("S\n");S.print();};
  void addDiagonal(double element);
  
private:
  int nP;
  rawMatrix<double> H;
  rawMatrix<double> S;
};

#endif
