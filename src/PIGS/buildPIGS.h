#ifndef BUILD_PIGS
#define BUILD_PIGS

#include "pigs.h"
#include "pigsMover.h"
#include "PIGSpropagator.h"
#include "measuresPIGS.h"

class pigsBuilder
{
public:
  typedef qmcDriver_t::configurations_t configurations_t;
  typedef system_t::pos_t pos_t;
  typedef system_t::geometry_t geometry_t;
  typedef system_t::rand_t rand_t;
  typedef qmcDriver_t::propagator_t propagator_t;
  
  pigsBuilder(qmcDriver_t * qmc_){qmcO=qmc_;};
  system_t::geometry_t* getGometry(){return qmcO->geo;}
  
  double getLbox(){return qmcO->geo->l_box;};
  double getTimeStep(){return qmcO->delta_tau;};
  qmcDriver_t::wave_t* getWave(){return qmcO->getWave();}
  void getParameter(string name,double & value);
  void build(configurations_t & configurations);
  propagator_t * getPropagator(){return qmcO->getMover()->getPropagator();}
  int getNBeads();
  int getN();
  system_t::rand_t* getRandomGenerator(){return qmcO->rand;}
  geometry_t * getGeometry(){return qmcO->geo;}
  
  void build(levyReconstructor & constructor);
  void build(pairParticleApproximationChain & prop);
  void build(freePropagator1D & gf);
  void build(energyThermodynamic & m,xml_input * main_input);
  void build(pairApproximationPropagator1DUnitary & gf);
  void build(wiggleMove & mov);
  void build(translateMove & mov);
  void build(moveVariationalTails & mov);
  void build(pigsMover & mov);
  void build(towerSampling & tower);
  void build(moveSingleBead & move);
  void build(pairCorrelationPigs & pairCorrelationPigs,xml_input * xmlElement);
  void build(measures<qmcDriver_t> & m);
  
private:
  
  qmcDriver_t * qmcO;
};

#endif


