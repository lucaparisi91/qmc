#ifndef PIGS_H
#define PIGS_H

#include "../qmc.h"
#include "../random.h"
#include "../wavefunction.h"
#include "configurations.h"
#include "pigsTraits.h"

class pigs_t{};

//advanced definitions
class pairParticleApproximationChain;
class pigsMover;

class pigsDriver : public qmc<system_t>
{
  
public:
  typedef pigs_t qmcKind;
  typedef typename system_t::rand_t rand_t;
  typedef typename system_t::geometry_t geometry_t;
  typedef double value_t;
  typedef system_t::all_particles_t all_particles_t;
  typedef system_t::particles_t particles_t;
  typedef total_wavefunction<pigsDriver >::wave_t swave_t;
  typedef total_wavefunction<pigsDriver > wave_t;
  typedef configurationsPIGS configurations_t;
  typedef system_t::pos_t pos_t;
  typedef configurations_t walker_t;
  typedef pairParticleApproximationChain propagator_t;
  typedef pigsMover pigsMover_t;
  
  void setWavefunction();
  
  pigsDriver()
  {
    printf("Setting the wavefunction....\n");
    setWavefunction();
    
  };
  
  void run();
  
  void load();
  
  void step();
  
  wave_t* getWave(){return wave;}
  
  void measure();
  
  void out();
  
  void setWavefunction(wave_t * wave_){wave=wave_;}
  configurations_t & getConfigurations(){return configurations;}
  
  int getNBeads(){return configurations.size();}

  pigsMover_t * getMover(){ return pigsMoverO;}
  
private:
  
  wave_t* wave;
  configurations_t configurations;
  measures<pigsDriver > estimates;
  pigsMover_t *pigsMoverO;
};

typedef pigsDriver qmcDriver_t;

#endif
