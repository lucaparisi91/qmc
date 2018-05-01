#ifndef PIGSPROPAGATOR_H
#define PIGSPROPAGATOR_H

#include "propagatorFunctions.h"

class pairParticleApproximationChain
{
public:
  
  typedef  freePropagator1D freeParticlePropagator_t;
  typedef pairApproximationPropagator1DUnitary pairParticlePropagator_t;
  typedef system_t::geometry_t geometry_t;
  typedef qmcDriver_t::wave_t wave_t;
  typedef system_t::pos_t pos_t;
  typedef qmcDriver_t::configurations_t configurations_t;
  
  //------------------initialization-----------------------------------
  double logEvaluateLinks( configurations_t &configurations);
  double logEvaluateSingleLinks( configurations_t &configurations);
  double logEvaluatePairLinks( configurations_t &configurations);
  void setWave(wave_t *wave_){wave=wave_;};
  
  void setFreePropagator(freeParticlePropagator_t * gr){freeParticlePropagatorFunction=gr;};
  void setPairPropagator(pairParticlePropagator_t * gr){pairParticlePropagatorFunction=gr;};
  
  //--------------------- evaluation----------------------------
  double logEvaluate(configurations_t & configurations);
  double logEvaluatePairLinks(int iParticle,int iTimeSliceBegin,int iTimeSliceEnd,configurations_t & configurations);
  double logEvaluateSingleLinks(int iParticle,int iTimeSliceBegin,int iTimeSliceEnd,configurations_t & configurations);
  double logEvaluateTail(configurations_t & configurations){return wave->log_evaluate(&configurations[0]);}
  double logEvaluateHead(configurations_t & configurations){return wave->log_evaluate(&configurations[configurations.size()-1]);}
  
  double logEvaluateSingleLinks(int iParticle,int jParticle,int iTimeSliceeBegin,int iTimeSliceEnd,configurations_t & configurations);
  double logEvaluatePairLinks(int iParticle,int jParticle,int iTimeSliceBegin,int iTimeSliceEnd,configurations_t & configurations);
  
  double logEvaluateTimeDerivative(configurations_t & configurations);
  
  void setGeometry(geometry_t* geo_){geo=geo_;}
  geometry_t* getGeometry(){return geo;}
  void setTimeStep(double step){timeStep=step;}
  
private:
  
  double timeStep;
  wave_t * wave;
  freeParticlePropagator_t *freeParticlePropagatorFunction;
  pairParticlePropagator_t *pairParticlePropagatorFunction;
  geometry_t * geo;
};
class linksTable
{
  double single(int iParticle,int iTimeSlice){return linksOneBody[M*(iParticle)+iTimeSlice];}
private:
  int M;
  vector<double> linksOneBody;
  vector<double> linksTwoBody;
  
};
#endif
