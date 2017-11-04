#ifndef TRAITS_H
#define TRAITS_H
#include "vector"
#include <complex>

using namespace std;
// define the dimensionality of the system
class jastrow_delta;
class spinors;
class jastrow_gaussian;
class random1;
template<class T> class qmc;
template<class T> class dmc;
template<class T> class vmc;
class jastrow_spline;
class pbc1d;
template<class T> class dmc_walker;
template<class T> class vmc_walker;
template<class pbct> class geometry;
template<class comp> class potential;
template<class T> class wavefunction;
template<class T> class total_wavefunction;
class jastrow_delta_bound_state;
class jastrow_delta_phonons;
class jastrow_delta_in_trap;
template < class> class particles;
template <class > class all_particles;
class noPbcD1;
typedef double real_t;
template <class comp> class speciesHarmonicPotential;
template <class comp> class empty_potential;
template <class,class,class> class jastrow_spinor_free_sampling;
template <class,class,class> class jastrow_spinor;
class dmc_kind
{
  
};
class vmc_kind
{
  
};
class D1_t
{
  
};
class empty_t
{};
class D2_t
{
  
};

class spinor1D
{
  
};
class spinor1DVMCMeasure
{
  
};
class spinor1DDMCMeasure
{
  
};
class D1VMCMeasure
{
  
};
class D1DMCMeasure
{
  
};
class D3_t
{
  
};


template<class X>
class traits
{
  
};

template<class X>
class particleTraits
{
  
};

template<>
class particleTraits<D1_t>
{
public:
  typedef particles<double> particles_t;
  typedef all_particles<particles_t> all_particles_t;  
};

template<>
class particleTraits<spinor1D>
{
public:
  typedef spinors particles_t;
  typedef all_particles<particles_t> all_particles_t;  
};





template<>
struct traits<D1_t>
{
  typedef real_t position_t;
  typedef real_t value_t;
  typedef int position_index_t;
  typedef geometry<pbc1d> geometry_t;
};

template<>
struct traits<spinor1D>
{
  typedef real_t position_t;
  typedef complex<double> value_t;
  typedef int position_index_t;
  typedef geometry<pbc1d> geometry_t;
};




template<class walker_t,class wave_t>
class measurementInterface;
template<class wave_t>
class rqmcMeasurementInterface;
class rqmc;
class rqmc_walker;
template<class X>
class qmcTraits
{
public:
  typedef qmc<X> qmc_t;
};
// template<>
// struct traits<rqmc>
// {
// public:
//   typedef rqmc_walker walker_t;
//   typedef total_wavefunction wave_t;
//   typedef rqmcMeasurementInterface<wave_t> measurementInterface_t;
//   typedef vector<walker_t *> measure_obj_t; 
// };
// inheriths some properties

// template<>
// struct traits<wave_traits>
// {
// public:
//   typedef wavefunction swave_t;
//   typedef total_wavefunction wave_t;
// };

// uses a vector of spinors in 1D
class spinorDmc1D
{
  
};

class dmc_t
{
};

class vmc_t
{
};

template<class X>
class dmcTraits
{
public:
  typedef X qmc_comp;
  typedef dmc<X> qmc_t;
  typedef dmc_walker<X> walker_t;
};

template<class X>
class vmcTraits
{
public:
  typedef X qmc_comp;
  typedef vmc<X> comp;
  typedef vmc_walker<X> walker_t;
};

// template<class X>
// class waveTraits
// {
// public:
//   typedef total_wavefunction<X> wave_t;
//   typedef wavefunction<X> swave_t;
// };

template<class X>
class valueTraits
{
  
};

template<>
class valueTraits<spinor1D>
{
public:
  typedef complex<double> value_t;
};

template<>
class valueTraits<D1_t>
{
public:
  typedef  double value_t;
};

template<class X>
class measuresTraits
{
public:
  typedef typename X::wave_t wave_t;
  typedef typename X::walker_t walker_t;
  typedef measurementInterface<walker_t,wave_t> measurementInterface_t;
  typedef walker_t measure_obj_t;
  
};

template<class T>
class randTraits
{
public:
  typedef random1 rand_t;
};

template<class T>
class potentialTraits
{
public:
  typedef potential<T> potential_t;
};





class keep_history_t
{
  
};

class no_spline
{
  
};

class with_spline
{
  
};
#endif
