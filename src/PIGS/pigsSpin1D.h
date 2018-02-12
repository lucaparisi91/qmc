#include <cstdio>
#include "../particles/orbital.h"
#include "../particles/gradientParticles.h"
#include "configurations.h"
#include "../wavefunction.h"
#include "../random.h"
#include "pigs.h"
#include "../geometry.h"

typedef configurations_t<orbitals<orbitals<spinOrbital1D> > > configurationsSpin_t;

class PIGSSpinSystem_t
{
public:
  typedef random1 rand_t;
  typedef configurationsSpin_t walker_t;
  typedef orbitals<orbitals<spinOrbital1D> > all_particles_t;
  typedef allParticlesGradient1D grad_t;
  typedef geometry<pbc1d> geometry_t;
  
};

typedef  pigsDriver<PIGSSpinSystem_t> pigsSpin1DDriver_t;
