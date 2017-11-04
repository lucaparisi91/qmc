#include "vmc.h"
#include "geometry.h"
#include "system.h"
#include "wavefunction.h"
#include "jastrow_collection.h"
#include <cmath>
#include "measures.h"
//perofrms a set of measurements on the system
/*
void vmc::make_measurements()
{
  vector <double> drift_force(geo->n_particles);
  double e;
  double e_f;
  m->increment();
  wave->kinetic_energy(p,e,e_f);// calculate the kinetic energy
  m->ms[0]->add(e,current_step);// add the local energy
  m->ms[1]->add(e_f,current_step); // add the force energy estimator
  
  m->make_measurements(p,wave);
  m->record(current_step);
}
*/
