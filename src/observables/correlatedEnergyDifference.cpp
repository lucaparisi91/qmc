#include "correlatedEnergyDifference.h"
#include "../vmc.h"
#include <iostream>
using namespace  std;
template<class walker_t,class wave_t>
correlatedEnergyDifference<walker_t,wave_t>::correlatedEnergyDifference(measure_vector* ms_,vector<int> & ns,wave_t* wave_,int nWaves) : measurement<walker_t,wave_t,measure_vector>(ms_)
{
  for(int i=0;i<nWaves;i++)
    {
      waves.push_back(new wave_t(wave_));
    }
  
  gradTmp.resize(ns);
  obs.resize(2*nWaves+1);
}

template<class walker_t,class wave_t>
void correlatedEnergyDifference<walker_t,wave_t>::make_measurement(walker_t* w,wave_t* wave)
{
  double ef2,psi2;
  
  for(int i=0;i<waves.size();i++)
    {
      
      waves[i]->laplacianGradient(*w->state,obs[2*i],ef2,gradTmp);
      psi2=exp(2*(waves[i]->log_evaluate(w->state)-w->wavefunction_value  ));
      obs[2*i]=(obs[2*i]+w->ev)*psi2;
      obs[2*i+1]=psi2;
    }
  obs[2*waves.size()]=w->e;
  
  this->ms->add(obs,0);
}

template class correlatedEnergyDifference< vmc_walker< D1_t<pbc1d > >,total_wavefunction< vmc<D1_t<pbc1d > > > >;
template class correlatedEnergyDifference< vmc_walker< D1_t<noPbcD1 > >,total_wavefunction< vmc<D1_t<noPbcD1 > > > >;
