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
  void setParameters(vector<double> params);
  void accumulate(walker_t* w);
  int getNParams(){return waves.size();} const
  
  void print() const;
  void setGradient(grad_t * g){gradTmp=g;};
  void transfer(int root){rawAccumulators.transfer(root);};
  
  optimizePlan& getPlan(){return plan;};
  
  int getStep(vector<double> &step);
  int getStep(vector<double> &step,double shift);
  void getMean(vector<double> &vecOut) const{rawAccumulators.getMean(vecOut);}
  void reset(){rawAccumulators.reset();}
  void getMeanError(vector<double> &meanOut,vector<double> &errorOut) const{rawAccumulators.getMeanError(meanOut,errorOut);};
    
private:
  
  grad_t* gradTmp;
  optimizePlan plan;
  vector<wave_t*> waves;
  vector<double> walkerMeasurements;
  vectorAccumulatorVariance<double> rawAccumulators;
  linearMethodStepEstimator stepEstimator;
  
};

#include "optimize.hpp"


#endif
