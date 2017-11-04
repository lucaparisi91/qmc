#include "optimizationObservablesLinearMethod.h"

template<class walker_t,class wave_t>
void optimizationObservablesLinearMethod<walker_t,wave_t>::make_measurement(walker_t* w,wave_t* wave)
{
    double e2,ef2,psiPrime,psi,psi2,ePrime;
    psi=w->wavefunction_value;
    
    wave2->laplacianGradient(*w->state,e2,ef2,gradTmp);
    e2+=w->ev;
    
    psi2=wave2->log_evaluate(w->state);
    // wavefunction parameter derivative
    psiPrime=(exp(psi2-psi)-1)/delta;
    // energy parameter derivative
    ePrime=(exp(psi2-psi)*e2-w->e)/delta-psiPrime*w->e;
    
    optObs[0]=w->e;
    optObs[1]=psiPrime*w->e;
    optObs[2]=psiPrime;
    optObs[3]=ePrime;
    optObs[4]=psiPrime*psiPrime*w->e;
    optObs[5]=psiPrime*ePrime;
    optObs[6]=psiPrime*psiPrime;
    
    this->ms->add(optObs,0);
    
  }

template<class walker_t,class wave_t>
double optimizationObservablesLinearMethod<walker_t,wave_t>::getEnergy()
  {
    return this->ms->getMean()[0];
  }


template<class walker_t,class wave_t>
vector<double>  optimizationObservablesLinearMethod<walker_t,wave_t>::getMean()
{
  return this->ms->getMean();
}

template<class walker_t,class wave_t>
void optimizationObservablesLinearMethod<walker_t,wave_t>::addParameter(int i,int j)
  {  
    parametersToOptimize.push_back(pair<int,int>(i,j));
  }

template<class walker_t,class wave_t>
void optimizationObservablesLinearMethod<walker_t,wave_t>::setParameter(double p,wave_t * waveTmp) const
  {
    int i,l1,l2;
    
    for(i=0;i<parametersToOptimize.size();i++)
      {
	
	l1=parametersToOptimize[i].first;
	l2=parametersToOptimize[i].second;
	waveTmp->setParameter(p,l2,l1);
	
      }
    
  }


template<class walker_t,class wave_t>
double optimizationObservablesLinearMethod<walker_t,wave_t>::getParameter()
{
  int l1,l2;
  l1=parametersToOptimize[0].first;
  l2=parametersToOptimize[0].second;
  return wave->getParameter(l2,l1);
}

template<class walker_t,class wave_t>
void optimizationObservablesLinearMethod<walker_t,wave_t>::init()
  {
    double p2;
    int l1,l2;
    
    // increase by delta the parameter of the second wavefunction
    if (parametersToOptimize.size()==0)
      {
	cout << " Nothing to optimize !"<<endl;
	exit(1);
      }
    
    p2=getParameter() + delta; 
    setParameter(p2,wave2);
    
  }

#include "../vmc.h"

template class  optimizationObservablesLinearMethod< vmc_walker< D1_t<pbc1d > >,total_wavefunction< vmc<D1_t<pbc1d > > > >;
template class  optimizationObservablesLinearMethod< vmc_walker< D1_t<noPbcD1 > >,total_wavefunction< vmc<D1_t<noPbcD1 > > > >;
