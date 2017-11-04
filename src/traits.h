#ifndef TRAITS_H
#define TRAITS_H
#include "vector"
#include <complex>
using namespace std;
template<class> class particles;
template<class> class all_particles;
template<class> class geometry;
class spinorsRealDmc;

class pbc1d;
class random1;
class spinors;

class empty_t{};

class keep_history_t{};

typedef random1 rand_t;
typedef double real_t;
template<class T>
class traits
{
};

class dmc_t{};
class vmc_t{};

// Definees the dimensionality of the system
template<class bc>
class D1_t
{
public:
  typedef double position_t;
  typedef double value_t;
  typedef particles<value_t> particles_t ;
  typedef all_particles<particles_t> all_particles_t;
  typedef geometry<bc> geometry_t;
};

// the spinor templated class
template<class bc>
class spinor
{
public:
  typedef complex<double> value_t;
  typedef geometry<bc> geometry_t;
  typedef spinors particles_t;
  typedef all_particles<particles_t> all_particles_t;
};

template<class bc>
class spinorReal
{
public:
  typedef double value_t;
  typedef geometry<bc> geometry_t;
  typedef spinorsRealDmc particles_t;
  typedef all_particles<particles_t> all_particles_t;
};

typedef spinorReal<pbc1d> spinor1D;

#include "jastrow_traits.h"


class VMCMitasSpinSampling
{
  
};

class VMCNoSpinSampling
{
  
};

template<class T>
class VMCMoveTraits
{
public:
  typedef VMCNoSpinSampling move_t;
};

template<>
class VMCMoveTraits<spinor1D>
{
 public:
  //typedef VMCNoSpinSampling move_t; // makes no move when making the spin sampling
  typedef VMCMitasSpinSampling move_t;
};

#endif
