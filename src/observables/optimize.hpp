#include "optimize.h"
#include "lapacke.h"
#include "../tools.h"
#include <cmath>

template<class wave_t,class walker_t>
void linearMethodOptimize<wave_t,walker_t>::setPlan(optimizePlan &planIn,const wave_t* wave_)
{  
  int k,iWave,iParam;
  plan=planIn;
  optimizePlan::iterator it;
  k=0;
  for(it=plan.begin();it!=plan.end();++it)
    {
      for(int i=0;i<it->second.size();i++)
	{
	  iWave=it->second[i].first;
	  iParam=it->second[i].second;
	  // creates new wavefunctions \psi(p_i + delta)
	  waves.push_back(new wave_t(wave_));
	  waves[k]->incrementParameter(iParam,iWave,plan.getDelta());
	  k+=1;
	}
    }
  
  k=waves.size();
  int nAccumulators=k*3+2*(k*(k+1))/2+k*k+1;
  walkerMeasurements.resize(nAccumulators);
  rawAccumulators.resize(nAccumulators);
  
}

/*
WARNING: assumes the wavefunction and energy have been updated
*/

template<class wave_t,class walker_t>
void linearMethodOptimize<wave_t,walker_t>::accumulate(walker_t* w)
{
  double psi;
  double psiDelta;
  double psiPrime;
  double eDelta;
  double ePrime;
  double dummy;
  double delta;
  vector<double> parameters;
  int k;
  
  delta=plan.getDelta();
  psi=w->wavefunction_value;
  
  // loop on parameters
  for(int i=0;i<waves.size();i++)
    {
      waves[i]->laplacianGradient(*w->state,eDelta,dummy,*gradTmp);
      waves[i]->getParameters(plan,parameters);
      //tools::print(parameters);
      eDelta+=w->ev;
      psiDelta=waves[i]->log_evaluate(w->state);
      psiPrime=(exp(psiDelta-psi)-1)/delta;
      ePrime=(exp(psiDelta-psi)*eDelta-w->e)/delta-psiPrime*w->e;
      walkerMeasurements[3*i]=psiPrime;
      walkerMeasurements[3*i+1]=psiPrime*w->e;
      walkerMeasurements[3*i+2]=ePrime;
      
    }
  
  k=waves.size()*3;
  
  // loop on couple of parameters(symmetric)
  for(int i=0;i<waves.size();i++)
    {
      for(int j=0;j<=i;j++)
	{
	  walkerMeasurements[k]=walkerMeasurements[3*i]*walkerMeasurements[3*j];
	  walkerMeasurements[k+1]=walkerMeasurements[k]*w->e;
	  
	  k+=2;
	}
    }
  
  // loop on couple of parameters(non symm-part)
  for(int i=0;i<waves.size();i++)
    {
      for(int j=0;j<waves.size();j++)
	{
	  walkerMeasurements[k]=walkerMeasurements[3*i]*walkerMeasurements[3*j+2];
	    k+=1;
	}
    }
  
  walkerMeasurements[k]=w->e;
  
  //tools::print(walkerMeasurements);
  rawAccumulators.accumulate(walkerMeasurements);
  
}

template<class wave_t,class walker_t>
const void linearMethodOptimize<wave_t,walker_t>::print() const
{
  
  printf("Raw measurements for optimization\n");
  vector<double> averages;
  vector<double> variances;

  double mean,error;
  
  rawAccumulators.getMean(averages);
  rawAccumulators.getVariances(variances);
  assert(averages.size()==variances.size());
  assert(averages.size()>0);
  
  for(int i=0;i<averages.size();i++)
    {
      mean=averages[i];
      assert(mean!=0);
      error=sqrt(abs(variances[i])/averages.size())/sqrt(rawAccumulators.getNmeasurements());
      
      printf("%f+-%f(%f%%)\n",mean,error,error/abs(mean)*100 );
    }
  printf("\n");
}
 
template<class wave_t,class walker_t>
int linearMethodOptimize<wave_t,walker_t>::getStep(vector<double> &params)
{
  vector<double> mean;
  rawAccumulators.getMean(mean);
  stepEstimator.buildMatrix(mean,getNParams());
  
  return stepEstimator.getStep(params);
}

template<class wave_t,class walker_t>
int linearMethodOptimize<wave_t,walker_t>::getStep(vector<double> &params,double shift)
{
  vector<double> mean;
  rawAccumulators.getMean(mean);
  stepEstimator.buildMatrix(mean,getNParams());
  stepEstimator.addDiagonal(shift);
  
  return stepEstimator.getStep(params);
}

template<class wave_t,class walker_t>
void linearMethodOptimize<wave_t,walker_t>::setParameters(vector<double> params)
{
  
  assert(params.size()==waves.size());
  for(int i=0;i<params.size();i++)
    {
      params[i]+=plan.getDelta();
    }
  
  for(int i=0;i<waves.size();i++)
    {
      waves[i]->setParameters(plan,params);
    }
  
}

