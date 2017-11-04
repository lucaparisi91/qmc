/*
  Implements the linear optimization method.
  Returns the parameter variation from a vector of estimators of energy and derivatives of the energy
  - implemented for just one variational parameter( to be extened)
 */

#ifndef OPTIMIZATIONMATRIXLINEARMETHOD_H
#define OPTIMIZATIONMATRIXLINEARMETHOD_H

#include <vector>
using namespace std;

class optimizationMethodLinear
{
public:
  optimizationMethodLinear();
  void buildHMatrix(vector<double> &obs);
  void buildSMatrix(vector<double> &obs);
  double getParameterVariation();
  double getParameterVariation(vector<double> &obs);
private:
  double * H;
  double * S;//overlap matrix
  int N;
  double* leftEigenVectors;
  double* rightEigenVectors;
  double* work;
  double* alphaI;
  double* alphaR;
  double* beta;
  
};

#endif
