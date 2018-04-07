#ifndef PIGS_H
#define PIGS_H

#include "../random.h"
#include "../wavefunction.h"
#include "configurations.h"
#include "measuresPIGS.h"
#include "PIGSpropagator.h"

class pigs_t
{
  
};

template<class comp>
class contactInteractions1DProp
{
public:
  typedef freePropagator1D freeParticlePropagator_t;
  typedef pairApproximationPropagator1DUnitary pairParticlePropagator_t;
  
  typedef typename comp::geometry_t geometry_t;
  typedef typename comp::configurations_t configurations_t;
  typedef typename comp::wave_t wave_t;
  typedef typename comp::pos_t pos_t;
  
};



template<class comp>
class pigsDriver : public qmc<comp>
{
  
public:
  typedef pigs_t qmcKind;
  typedef typename comp::rand_t rand_t;
  typedef typename comp::geometry_t geometry_t;
  typedef double value_t;
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename all_particles_t::particles_t particles_t;
  typedef typename total_wavefunction<pigsDriver<comp> >::wave_t swave_t;
  typedef total_wavefunction<pigsDriver<comp> > wave_t;
  typedef typename comp::configurations_t configurations_t;
  typedef typename comp::pos_t pos_t;
  typedef configurations_t walker_t;
  
  typedef pairParticleApproximationChain< contactInteractions1DProp<pigsDriver<comp>   > > propagator_t;
  
  void setWavefunction()
  {
    vector<swave_t*> waves;
    load_wavefunctions<pigsDriver<comp> >(this->main_input,waves,this);
    wave=new wave_t(this);
    wave->link_wavefunctions(this->main_input,waves,"main");
    wave->print_jastrows();
  }

  pigsDriver()
  {
    printf("Setting the wavefunction....\n");
    setWavefunction();

  };
  
  void run()
  {
    printf("PIGS----------\n");
    load();
    
    step();
    
  };
  
  void step()
  {
    
    pigsMoverO.move(configurations);
    
  }
  
  void measure()
  {
    estimates.make_measurements(&configurations,wave);
    estimates.increment();
  }
  
  void out()
  {
    estimates.out();
  }
  
  void load()
  {
    printf("...loading configurations\n");
    buildConfigurationsPIGS<pigsDriver<comp> > configBuilder;
    printf("...setting measurements\n");
    configBuilder.initConfigurations(this);
    buildPIGSMeasures("input.xml",estimates);
  }
  
  void setWavefunction(wave_t * wave_){wave=wave_;}
  configurations_t & getConfigurations(){return configurations;}

  void getNBeads(){return configurations.size();}

private:

  wave_t* wave;
  configurations_t configurations;
  measures<pigsDriver<comp> > estimates;
  pigsMover< pigsDriver<comp> > pigsMoverO;
};

#endif
