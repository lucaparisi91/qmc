#ifndef PIGS_H
#define PIGS_H

#include "../random.h"
#include "../wavefunction.h"

template<class comp>
class pigsDriver : public qmc<comp>
{
public:
  typedef typename comp::rand_t rand_t;
  typedef double value_t;
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename all_particles_t::particles_t particles_t;
  
  typedef typename total_wavefunction<pigsDriver<comp> >::wave_t wave_t;
  
  pigsDriver(){};
  void run()
  {
    printf("PIGS----------\n");
    printf("PIGS----------...running");
  };
  
  void load();
private:
  
  wave_t* wave;
  
};


#endif

