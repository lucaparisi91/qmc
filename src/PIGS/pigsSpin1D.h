#include <cstdio>
#include "../particles/orbital.h"
#include "../particles/gradientParticles.h"
#include "configurations.h"
#include "../wavefunction.h"
#include "../random.h"
#include "pigs.h"
#include "../geometry.h"




class PIGSSpinSystem_t
{
public:
  typedef random1 rand_t;
  
  typedef orbitals<orbitals<spinOrbital1D> > all_particles_t;
  typedef all_particles_t::particles_t particles_t;
  typedef allParticlesGradient1D grad_t;
  typedef geometry<pbc1d> geometry_t;
  typedef configurationsPIGS_t<all_particles_t> configurations_t;
  typedef particles_t::particles_t::pos_t pos_t;
  
};

typedef  pigsDriver<PIGSSpinSystem_t> pigsSpin1DDriver_t;
