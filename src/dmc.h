#ifndef DMC_H
#define DMC_H
#include <vector>
#include <fstream>
#include "qmc.h"
#include "omp.h"
#include "mpi.h"
#include "measures.h"
#include "traits.h"
#include "qmcDriver/moves.h"
#include "wavefunction.h"
#include "particles/gradientParticles.h"
#include "qmcDriver/qProbability.h"
#include "particles/orbital.h"

using namespace std;

template<class T> class walkers;
template<class> class dmc_walker;
template<class> class rabi_walker;
class gatherer;
class packed_data;
class dock;
template <class comp > class measure_dynamic;

bool metropolis(double log_ratio,random1* rand);
// return

template<class comp>
class dmc : public qmc<comp>
{
 public:
  typedef random1 rand_t;
  typedef typename comp::spinSampling_t spinSampling_t;
  typedef typename walkerTraits<comp,spinSampling_t>::dmc_walker_t walker_t;
  typedef typename comp::value_t value_t;
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename comp::particles_t particles_t;
  typedef wavefunction< dmc<comp> > swave_t;
  typedef total_wavefunction< dmc<comp> > wave_t;
  typedef qmcMover<dmc<comp> > qmcMover_t;
  typedef allParticlesGradient1D grad_t;
  typedef comp qmcSystem_t;
  
  typedef potential< dmc<comp> > potential_t;
  typedef dmc_t qmcKind;
  
  dmc();
  int i_walker; // index of the current walker
  xml_input* xml_walkers;
  xml_input* walkers_load_xml;
  dock* dmc_dock;
  int total_walkers,n_walkers_c;
  int mean_walkers;//mean value of walkers
  double delta_tau_or;
  char * inBuffer;
  double delta_walkers;//uncertainity on the walkers
  walkers<comp> *ws;
  //measures* m;
  double e_t;// population control energy
  void run();// run the algorithm
  void step();// performs a single mc step
  void warmup_step();
  void out();// perform outputs to file (and screen if VERBOSE)
  void load();
  void reduce();
  void save();
  void saveParallel(const string & filename);
  void send_mpi();
  void receive_mpi();
  
  void receive_mpi_tmp();
  void send_walkers();
  void collect();
  void ready();
  void population_control();
  int detailed_balance; // whatever to perform an accept-reject step in the diffusion monte carlo
  bool smart_vmc;
  
  void send_receive_determine();
  void pack();
  void unpack();
  void send_packs();
  void recv_packs();
  void move2order(all_particles_t* p,all_particles_t* ptmp);
  
  vector<swave_t*> waves;
  wave_t * wave;
  qmcMover_t* qmcMoverO; // perform updates of the position on the walkers
  
  potential_t* potential_obj;
private:
  vector<double> tmp;
  vector<double> tmp1;
  
};

template<class comp>
class walker
{
public:
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename comp::value_t value_t;
  typedef double waveValue_t;
  typedef allParticlesGradient1D gradient_t;
  
  double wavefunction_value;// wavefunct
  all_particles_t* state;//the state of the physical system
  all_particles_t* state_tmp;
  double e;// the energy
  double ev;//the potential energy
  double e_f;
  double e_old;// a tmp for computing the descendants
  //vector<double> drift_force; // allows to compute the drift force of the system
  double e_difference;
  double weight;
  bool accept;
  
  
  void log(); 
 public:
  // makes the diffusion step and mark for the branching step
  waveValue_t getWavefunctionValue(){return wavefunction_value;};
  void send(gatherer *ws_global);
  void send_mpi(int dest,int tag,vector<MPI_Request> &requests );
  void receive_mpi(int dest,int tag,vector<MPI_Request> &requests);
  template<class qmc_t>
  void save(qmc_t* qmc_obj);
  template<class qmc_t>
  void set(qmc_t* qmc_obj);
  void set_particles(all_particles_t* p);
  void receive(gatherer *ws_global);
  void printGeneral();
  void pack(packed_data* walker_pack);
  void unpack(packed_data* walker_pack);
  template<class walker2_t>
  void clone(walker2_t & w);

  template<class comp1>
  friend ostream& operator<<(ostream& out, const walker<comp1> & walker);

  template<class comp1>
  friend istream& operator>>(istream& in, walker<comp1> & walker);

  
  gradient_t & getParticlesGradient()
  {
    return particlesGradient;
  }
  
  int get_pack_size();
protected:
  template<class qmc_type>
  walker(qmc_type* qmc_obj);
protected:
  walker<comp>& operator=(walker<comp> & w){};
  gradient_t particlesGradient; // stores the gradient of the wavefunction
  
  gradient_t particlesGradientBackup;
  
};

template<class comp>
class dmc_walker : public walker<comp>
{
public:
  int descendants;
  typedef dmc<comp> qmc_t;
  template<class qmc_type> void update(qmc_type* qmc_obj);
  template<class qmc_type>
  dmc_walker(qmc_type* dmc_obj_);
  int get_pack_size();
  double get_weight(const double e_ref,dmc<comp>*  dmc_obj);
  
  template<class qmc_type>
  void set(qmc_type* qmc_obj);
  void pack(packed_data* walker_pack);
  void unpack(packed_data* walker_pack);
  template<class qmc_type>
  void set_descendants(double &e_t,qmc_type* dmc_obj);
  template<class measure_t,class wave_t>
  void make_measurements(measure_t* m,wave_t* wave,double current_step);
  dmc_walker<comp>& operator=(dmc_walker<comp> & w);
  void print();
  vector< measure_dynamic<dmc<comp> > * > md;
};

template<class comp>
class rabi_walker : public dmc_walker<comp>
{
public:
  typedef dmc<comp> qmc_t;
  typedef typename dmc_walker<comp>::all_particles_t all_particles_t;
  typedef typename qmc_t::particles_t particles_t;
  rabi_walker(qmc_t* dmc_obj_) : dmc_walker<comp>(dmc_obj_)
  {
    
  }
  
  
  void pack(packed_data* walker_pack)
  {
    dmc_walker<comp>::pack(walker_pack);
    walker_pack->pack(spinFlipRatios);  
  };
  
  void unpack(packed_data* walker_pack)
  {
    dmc_walker<comp>::unpack(walker_pack);
    walker_pack->unpack(spinFlipRatios);
  }
  
  template<class measure_t,class wave_t>
  void make_measurements(measure_t* m,wave_t* wave,double current_step)
  {
    m->make_measurements(this,wave);
  }
  
  virtual void set(qmc_t * qmcO)
  {
    dmc_walker<comp>::set(qmcO);
    spinFlipRatios.resize((*this->state)[0].size(),0);
    qmcO->wave->spinFlip(*this->state,spinFlipRatios);
  }
  
  virtual rabi_walker<comp> & operator=(rabi_walker<comp> & w)
  {
    dmc_walker<comp>::operator=(w);
    spinFlipRatios=w.spinFlipRatios;
  }
  
  void updateAllatOnce(qmc_t* dmc_obj)
  { 
    dmc_obj->wave->spinFlip(*this->state,spinFlipRatios);
    dmc_obj->qmcMoverO->moveSpinRabi((*this->state)[0],spinFlipRatios);
    dmc_walker<comp>::update(dmc_obj);
  }
  
  void update(qmc_t* dmc_obj)
  {
    updateOneByOne(dmc_obj);
  }
  
  void updateOneByOne(qmc_t* dmc_obj)
  {
    dmc_walker<comp>::update(dmc_obj);
    
    double spinFlipRatio;
    particles_t & p1=(*this->state)[0];
    
    for(int i=0;i<p1.size();i++)
      {
	spinFlipRatio=dmc_obj->wave->spinFlip(i,*this->state);
	dmc_obj->qmcMoverO->moveSpinRabi(p1[i],spinFlipRatio);
      }
    
    
  }
  
private:
  
  vector<double> spinFlipRatios;
  
};


template<class comp>
class dmc_walker_fixed_node : dmc_walker<comp>
{
public:
  
  typedef dmc<comp> qmc_t;
  
  dmc_walker_fixed_node(qmc_t* dmc_obj_) : dmc_walker<comp>(dmc_obj_){};
  dmc_walker_fixed_node<comp>& operator=(dmc_walker_fixed_node<comp> & w);
  
  template<class qmc_type>
  void set_descendants(double &e_ref,qmc_type* dmc_obj);
  
  template<class qmc_type>
  void update(qmc_type* dmc_obj);
  
  int sign;
  int sign_current;
  
};


template<class comp>
class dmc_walker_fixed_phase : dmc_walker<comp>
{
public:
  typedef dmc<comp> qmc_t;
  
  dmc_walker_fixed_phase(qmc_t* dmc_obj_) : dmc_walker<comp>(dmc_obj_){};
  dmc_walker_fixed_phase<comp>& operator=(dmc_walker_fixed_phase<comp> & w);
  template<class qmc_type>
  void set_descendants(double &e_ref,qmc_type* dmc_obj);

  template<class qmc_type>
  void update(qmc_type* dmc_obj);
  
  double phase;
  double phase_current;
  
};

template<class comp>
class walkers
{
public:
  typedef typename dmc<comp>::walker_t walker_t;
  
  // rename the types of a certain class of measurements
  
  typedef measures<dmc<comp> > measures_t;
  
  vector<double> work;
  int n; // number of walkers
  measures_t *m;// collects all required measurements
  vector<walker_t *> ws; // vector of all walkers
  template<class qmc_type>
  walkers(qmc_type* dmc_obj);// initiate the walkers object
  template<class qmc_t> void branch(int &i,qmc_t* dmc_obj);// kills and copy as marked
  template<class qmc_t>
  void generate_random(int n_to_add,qmc_t* dmc_obj);
  //void save();
  //void send(gatherer* g);
  void send_mpi();
  //void receive(gatherer* g);
  void print();
  void log();
  template<class qmc_t>
  void generate_all_to(int n_to_add,double pos,qmc_t* dmc_obj);
  template<class qmc_t>
  void generate_random(int n_to_add,double lBox,qmc_t* dmc_obj);
  template<class qmc_t>
  void generate_uniform(int n_to_add,double lBox,qmc_t* dmc_obj);
  template<class qmc_t>
  void generate_gaussian(int n_to_add,double sigma,double position,qmc_t* dmc_obj);
  void save(const string &filename);
  string saveToString(const string &filename);
  void saveAppend(const string & filename);
  
  template<class comp1>
  friend ostream& operator<<(ostream& out,const walkers<comp1> & walkers);
  
  template<class comp1>
  friend istream& operator<<(istream& in,walkers<comp1> & walkers);

  
  
  
  
  
};

// store temporary variables when performing comunication among walkers

/* class gatherer */
/* { */
  
/*   public: */
/*   double e_t; */
/*   double current_step; */
/*   double n_metropolis; */
/*   double success_metropolis; */
/*   omp_lock_t gethering_lock; */
/*   int n_walkers; */
/*   int n; */
/*   bool saving; */
/*   ifstream load_file; */
/*   ofstream save_file; */
/*   vector<walker *> ws; */
/*   gatherer(); */
/* }; */

// includes the definition of templated functions [ required by the compiler]

template<class comp>
class moveQMC
{
public:
  template<class all_particles_t,class rand_t>
  void move(all_particles_t* p,rand_t* rand_o,const double &delta_tau)
  {
    moveGaussian(p,rand_o,delta_tau);
  }
};
// spin sampling according to mitas
template<>
template<class all_particles_t,class rand_t>
void moveQMC<VMCMitasSpinSampling>::move(all_particles_t* p,rand_t* rand_o,const double &delta_tau)
  {
    moveGaussian(p,rand_o,delta_tau);
    GaussianMoveSpin(p,rand_o,delta_tau);
  }

#include "walker.hpp"
#include "moves.hpp"
#include "dmc.hpp"

#endif
