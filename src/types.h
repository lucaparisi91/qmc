#ifndef TYPES_H
#define TYPES_H

#include<vector>
using namespace std;


typedef double real_t;
// there are differences between the two wavefunction sets
struct wavefunction_diff_sets_t{};
struct wavefunction_same_sets_t{};
// types defining the number of dimensions in the set
struct 1D_t{};
struct 2D_t{};
struct 3D_t{};


// one body jastrow wavefunction
struct one_body_jastrow_t{};
// two body jastrow wavefunction
struct two_body_jastrow_t{};

// template for the drift force
template<class T>
struct all_drift_force_t<1D_t,T>
{
  typedef T value_t;
  typedef vector<value_t> type;
};

// particle position type
template<class d>
particle_position_t
{
  
  typdef vector<real_t> type;
  
};

template<>
particle_position_t<1D_t>
{
  typdef real_t type;
};

// type of the drift force
template<class X>
struct all_drift_force_t
{
  typdef vector<X> type;
};

template<class d,class T>
// general traits of the wavefunction
struct wavefunction_traits
{
  typename typedef d dimension_t;
  typename typedef T value_t;
  typedef drift_force_t<d,value_t>::type drift_force_t;
  typedef all_drift_force_t<drift_force_t>::type all_drift_force_t;
  typedef particle_position_t<d>::type pos_t;//position of the particle
  typedef particles<typename d> particles_t;
  typedef all_particles<typename d> all_particles_t;
  typedef geometry<typename d> geometry_t; // geometry object type
  
};
// implements a periodic boundary condition
struct periodic_boundary_condition_t
{
  
};

#endif
