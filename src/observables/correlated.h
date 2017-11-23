#ifndef CORRELATED_H
#define CORRELATED_H

#include <vector>
using namespace std;

template<class wave_t,class walker_t>
class correlatedEstimatorEnergyDifference
{
  typedef typename walker_t::grad_t grad_t;
public:
  
  void setWavefunctions(vector<wave_t *> wavesIn){waves=wavesIn;resize(waves.size());};
  void addWavefunction(wave_t* waveIn){waves.push_back(waveIn);resize(waves.size());}
  void resize(int n){scratch.resize(n+1);scratchWeights.resize(n+1);waves.resize(n);accumulator.resize(n+1);}
  
  void accumulate(walker_t * walker)
  {
    double e;
    double dummy;
    
    for(int i=0;i<waves.size();i++)
      {
	waves[i]->laplacianGradient(*walker->state,e,dummy,*gradTmp);
	e+=walker->ev;
	
	scratchWeights[i]=exp(waves[i]->log_evaluate(walker->state)-walker->wavefunction_value);
	scratchWeights[i]*=scratchWeights[i];
	scratch[i]=e;
      }
    
    scratch[waves.size()]=walker->e;
    scratchWeights[waves.size()]=1.;
    accumulator.accumulate(scratch,scratchWeights);
    
  }
  
  void getMean(vector<double> &mean) const {accumulator.getMean(mean);} ;
  
  void getVariances(vector<double> &variance) const {accumulator.getVariances(variance);}
  vector<wave_t *> & getWaves(){return waves;};
  
  int getNmeasurements()
  {
    return accumulator.getNmeasurements();
  }
  
  void setGradient(grad_t * g){gradTmp=g;};
  
  void reset(){accumulator.reset();}
  void transfer(int root){accumulator.transfer(root);}
  
  void print()
  {
    
    vector<double> energies;
    vector<double> variances;
    int nMeasurements;
    double mean,error;
    accumulator.getMean(energies);
    accumulator.getVariances(variances);
    
    nMeasurements=accumulator.getNmeasurements();
    assert(nMeasurements>0);
    
    for(int i=0;i<energies.size();i++)
      {
	mean=energies[i];
	error=sqrt(abs(variances[i]))/sqrt(nMeasurements);
	if (mean!=0)
	  {
	    printf("%f+-%f(%% %f)\n",mean,error,error/mean*100);
	  }
	else
	  {
	    printf("%f+-%f\n",mean,error);
	  }
	
      }
  }
  
  int getMinCorrelatedEnergy()
  {
    vector<double> mean;
    accumulator.getMean(mean);
    assert(mean.size()>0);
    double eMin,eTrial;
    int iMin;
    eMin=0;
    eTrial=mean[0]-mean[mean.size()-1];
    iMin=0;
    
    for(int i=0;i<waves.size()-1;i++)
      {
	eTrial=mean[i]-mean[mean.size()-1];
	if (eTrial<eMin)
	  {
	    eMin=eTrial;
	    iMin=i;
	  }
      }
    
    return iMin;
  }
private:
  
  grad_t *gradTmp;
  vector<wave_t *> waves;
  vector<double> scratchWeights;
  vector<double> scratch;
  vectorAccumulatorVariance<double> accumulator;
};

#endif
