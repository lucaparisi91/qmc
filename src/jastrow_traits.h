
class jastrow_delta;
class jastrow_delta_phonons;
class jastrow_gaussian;
class jastrow_delta_bound_state;
class jastrow_delta_bound_state_no_pbc;
class jastrow_delta_bound_state_no_pbc2;
class jastrow_delta_bound_state_no_pbc3;
class jastrow_delta_in_trap;
class jastrow_delta_in_trap_exponential;
class jastrow_spline;
class jastrowSpinOrbital;
template<class,class,class> class jastrow_spinor;
template<class,class,class> class jastrow_spinor_free_sampling;
template<>
struct traits<jastrow_delta>
{
  typedef double value_t;
  typedef double position_t;
};

template<>
struct traits<jastrow_delta_phonons>
{
  typedef double value_t;
  typedef double position_t;
};

template<>
struct traits<jastrowSpinOrbital>
{
  typedef double value_t;
  typedef double position_t;
};

template<>
struct traits<jastrow_gaussian>
{
  typedef double value_t;
  typedef double position_t;
  
  
};

template<>
struct traits<jastrow_spinor<jastrow_delta,jastrow_delta,jastrow_delta> >
{
  typedef double value_t;
  typedef double position_t;
};
template<>
struct traits<jastrow_spinor_free_sampling<jastrow_delta,jastrow_delta,jastrow_delta> >
{
  typedef double value_t;
  typedef double position_t;
};

template<>
struct traits<jastrow_delta_bound_state>
{
  typedef real_t value_t;
  typedef real_t position_t;
};

template<>
struct traits<jastrow_delta_bound_state_no_pbc>
{
  typedef real_t value_t;
  typedef real_t position_t;
};

template<>
struct traits<jastrow_delta_bound_state_no_pbc2>
{
  typedef real_t value_t;
  typedef real_t position_t;
};

template<>
struct traits<jastrow_delta_bound_state_no_pbc3>
{
  typedef real_t value_t;
  typedef real_t position_t;
};

template<>
struct traits<jastrow_delta_in_trap>
{
  typedef real_t value_t;
  typedef real_t position_t;
};

template<>
struct traits<jastrow_delta_in_trap_exponential>
{
  typedef real_t value_t;
  typedef real_t position_t;
};


template<>
struct traits<jastrow_spline>
{
  typedef real_t value_t;
  typedef real_t position_t;
};

