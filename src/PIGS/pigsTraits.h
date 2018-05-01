#ifndef PIGS_TRAITS_H
#define PIGS_TRAITS_H

#include "../geometry.h"
#include "../random.h"
#include "../particles/orbital.h"
class spinor1D_t
{
public:
  typedef random1 rand_t;
  typedef geometry<pbc1d> geometry_t;
  typedef orbitals<orbitals<spinOrbital1D> > all_particles_t;
  typedef all_particles_t::particles_t particles_t;
  typedef double pos_t;
};

typedef spinor1D_t  system_t;


#endif
