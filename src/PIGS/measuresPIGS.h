#ifndef MEASURESPIGS_H
#define MEASURESPIGS_H

#include "../measures.h"
#include "pigsTools.h"

class energyTail : public measurement_scalar<qmcDriver_t::walker_t,qmcDriver_t::wave_t>
{
public:
  typedef qmcDriver_t::wave_t wave_t;
  typedef qmcDriver_t::walker_t walker_t;
  
  typedef typename wave_t::grad_t grad_t;
  
  energyTail(measure_scalar* mScal,string filename) :  measurement_scalar<walker_t,wave_t>::measurement_scalar(mScal){init(filename);};
  
  void init(string filename)
  {
    gradParticles=*buildAllParticlesGradient1D(filename);
  }
  
  void make_measurement(walker_t* w,wave_t* wave)
  {
    double e, dummy;
    wave->laplacianGradient((*w)[0],e,dummy,gradParticles);
    this->ms->add(e,0);
  }
  
private:
  grad_t gradParticles;
};

class energyThermodynamic : public measurement_scalar<qmcDriver_t::walker_t,qmcDriver_t::wave_t>
{
public:
  typedef qmcDriver_t::wave_t wave_t;
  typedef qmcDriver_t::walker_t walker_t;
  
  typedef typename wave_t::grad_t grad_t;
  typedef typename qmcDriver_t::propagator_t propagator_t;
  
  energyThermodynamic() :  measurement_scalar<walker_t,wave_t>::measurement_scalar(){};
  
  void make_measurement(walker_t* w,wave_t* wave);
  
  void setPropagator(propagator_t * p){prop=p;}
  
private:
  propagator_t *prop;
  
};



class pairCorrelationPigs :  public measurement<qmcDriver_t::walker_t,qmcDriver_t::wave_t,space_measure<measure_vector> >
{
public:
  typedef qmcDriver_t::configurations_t walker_t;
  typedef qmcDriver_t::wave_t wave_t;
  
  pairCorrelationPigs() : measurement<walker_t,wave_t,space_measure<measure_vector> >() {setA=0,setB=0;};
  
  virtual void make_measurement(walker_t * w,wave_t * wave);
  
  void setTimeSliceBegin(int i){iTimeSliceBegin=i;};
  void setTimeSliceEnd(int i){iTimeSliceEnd=i;}
  
private:
  
  int iTimeSliceBegin;
  int iTimeSliceEnd;
  int setA;
  int setB;
};

#endif
