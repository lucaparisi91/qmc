#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include "optimizePlan.h"
#include "optimizeTools.h"
#include "accumulator.h"

// linear optimization mwthod
template<class wave_t,class walker_t>
class linearMethodOptimize
{
public:
  typedef typename walker_t::grad_t grad_t;
  
  void setPlan(optimizePlan &planIn,const wave_t* wave);
  void accumulate(walker_t* w);
  int getNParams(){return waves.size();} const
  
  void print() const;
  void setGradient(grad_t * g){gradTmp=g;};
  void transfer(int root){rawAccumulators.transfer(root);};
  
  optimizePlan& getPlan(){return plan;};
  
  int getStep(vector<double> &step);
  
  void addDiagonal(double element){stepEstimator.addDiagonal(element);}
  
private:
  
  grad_t* gradTmp;
  optimizePlan plan;
  vector<wave_t*> waves;
  vector<double> walkerMeasurements;
  vectorAccumulatorVariance<double> rawAccumulators;
  linearMethodStepEstimator stepEstimator;
  
};

template<class wave_t,class walker_t>
class scalarAccumulatorVarianceCorrelated
{
public:
  void setWavefunctions(vector<wave_t *> waves);
  void accumulate(double value,double waveValue);
  void getMean(vector<double> &mean);
  void getVariances(vector<double> &variance);
private:
  vector<wave_t *> waves;
  vectorAccumulatorVariance<double> accumulator;
  
};

#include "optimize.hpp"


#endif
