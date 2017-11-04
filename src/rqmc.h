#include "measures.h"
#include "qmc.h"
#include "dmc.h"
/*
  *****************Reptation Quantum Monte Carlo******

*/

class qmc;
class rqmc_walker : public walker
{
public:
  template<class qmc_type> void update(qmc_type* mc_obj,rqmc_walker * wOld);
  rqmc_walker(qmc* qmc_ob);
  template<class measure_kind,class wavefunction_type>
  void make_measurements(measure_kind* m,wavefunction_type* wave,double current_step);
  void print();
  
};

class path
{
public:
  path(qmc* qmc_ob,int n_);
  void load(qmc* qmc_ob); // load a path from the walkers.xml file
  vector<rqmc_walker *> pathX;
  void save();
  inline int get_n(){return n;};
  inline void turn_switch(){pathSwitch=pathSwitch*(-1);};
  inline int get_switch(){return pathSwitch;};
private:
  int n;
  int pathSwitch;
};
// RQMC measurements

class rqmc : public qmc
{
public:
  rqmc();
  void step() ; // performs a variational like step
  void load(); // load initial VMC configuration
  void run(); // begin the calculation
  int n_metropolis;
  int success_metropolis;
private:
  path* rPath;
  measures<rqmc> *m; // gather measurements at the end of each path
  int time_slices; // the number of steps in a path
  int nBlocks;
  int warmupVMC;
  int deltaY;
  int reverse;  
};
