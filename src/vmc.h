#ifndef VMC_H
#define VMC_H

#include <vector>
#include <string>
#include "wavefunction.h"
#include "observables/optimize.h"
#include "observables/correlated.h"

//**************** Variational Monte Carlo [ using the Metropolis Algorithm]

template<class comp> class vmc;

template<class comp>
class vmc_walker : public walker<comp>
{
public:
  
  double get_weight(const double e_ref,vmc<comp>*  dmc_obj);
  void print();
  template<class qmc_t>
  vmc_walker(qmc_t* qmc_obj) : walker<comp>(qmc_obj) {};
  template<class qmc_type> void update(qmc_type* vmc_obj);
  template<class measure_kind,class wavefunction_type>
  void make_measurements(measure_kind* m,wavefunction_type* wave,double current_step);
  typedef allParticlesGradient1D grad_t;
  
};

#include "observables/optimizationMatrixLinearMethod.h"
#include "observables/correlatedEnergyDifference.h"

template<class comp>
class vmc : public qmc<comp>
{
public:
  double e_t;
  typedef typename comp::value_t value_t;
  typedef wavefunction< vmc<comp> > swave_t;
  typedef total_wavefunction< vmc<comp> > wave_t;
  typedef vmc_walker<comp> walker_t;
  typedef typename comp::all_particles_t all_particles_t;
  typedef measures< vmc<comp> > measures_t;
  typedef typename VMCMoveTraits<comp>::move_t move_t;
  typedef vmc_t qmcKind;
  typedef random1 rand_t;
  typedef allParticlesGradient1D grad_t;
  typedef typename all_particles_t::particles_t particles_t;
  typedef qmcMover1Order<vmc<comp> > qmcMover_t;
  typedef linearMethodOptimize<wave_t,walker_t>  mesOpt_t;
  
  typedef correlatedEstimatorEnergyDifference<wave_t,walker_t> mEnergyCorrelated_t;
  
  int warmupOptimizeSteps;
  vmc();
  
  void step();
  void optimizationStep();
  void correlatedEnergyOptimizationStep();
  void save();
  void saveAddWalker(); // prepend the current walkers [ to be used as input in a dmc calculation]
  void run();
  void runOptimize();
  void out();
  void optimizationOut();
  void optimizationAfterCorrOut();
  void load();
  void runStandard();
  void warmup_step();
  void stabilizationStep();
  void stabilizationOut();
  
  // step in which gather information about optimization
  //void optimizationStep();
  
  walker_t* w;
  measures_t* m;
  bool inputDmc;
  qmcMover_t moveEngine;
  wave_t* wave;
  vector<swave_t*> waves;
  
  void setOptimize(bool optimize_) {optimize=optimize_;};
  bool isOptimize() const {return optimize;}
  void chooseNextOptimizationParameters();
private:
  mEnergyCorrelated_t *mEnergyCorrelated;
  bool optimize;
  mesOpt_t mO;
  int unCorrelationSteps;
  int correlatedEnergySteps;
  int statusGeneralizedEigenValue;
  int statusCorrelatedMeasurements;
  vector<double> optParameters;
  vector< vector<double> > parametersProposal;
  int indexMinEnergyProposal;
  enum mode{absErrMode=1,countStepsMode=0};
  mode optimizationMode;
  double absErrorLimit;
  int statusCorrelated;
};

#include "walker_vmc.hpp"
#include "vmc.hpp"

#endif
