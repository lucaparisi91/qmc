#ifndef CORRELATEDENERGYDIFFERENCE_H
#define CORRELATEDENERGYDIFFERENCE_H
#include "../measures.h"

template<class walker_t,class wave_t>
class correlatedEnergyDifference : public measurement<walker_t,wave_t,measure_vector >
{
public:
  typedef typename wave_t::grad_t grad_t;
  
  correlatedEnergyDifference(measure_vector* ms_,vector<int> & ns,wave_t* wave_,int nWaves);
  
  void make_measurement(walker_t* w,wave_t * wave);

  // returns the difference between the current wavefunction and the stored wavefunctions
  
  vector<double> getMean()
  {
    vector<double> energies;
    
    obs=this->ms->getMean();
    energies.resize(waves.size());
    
    for(int i=0;i<waves.size();i++)
      {
	energies[i]=obs[2*i]/obs[2*i+1];
      }
    
    return energies;
  }
  
  
  wave_t* getWave(int i){return waves[i];}
  
private:
  grad_t gradTmp;
  vector<wave_t *> waves;
  vector<double> obs;

};


template<class walker_t,class wave_t>
correlatedEnergyDifference<walker_t,wave_t>* buildCorrelatedEnergy(xml_input* main_input,wave_t * wave)
{
  
  vector<int> ns;
  correlatedEnergyDifference<walker_t,wave_t>* objPtr;
  
  // gets the total number of particles(required to know the total number of particles)
  
  main_input->reset()->get_child("system")->get_first_child();
  
  while(main_input->check() )
    {
      if (main_input->get_name()=="particles")
	{
	  main_input->get_attribute("n");
	  ns.push_back(main_input->get_int());
	}
      
      main_input->get_next();
    }
  
  // create correlated energy differences
  objPtr=new correlatedEnergyDifference<walker_t,wave_t>(new measure_vector(2*3+1,"optimizationCorr",false),ns,wave,3);
  
  main_input->reset();
  
  //allocate mememory and initialize the auxiliary wavefunction
  
  return objPtr;
}

#endif

