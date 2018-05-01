#include "pigs.h"
#include "measuresPIGS.h"
#include "PIGSpropagator.h"
void pairCorrelationPigs::make_measurement(pairCorrelationPigs::walker_t * w,pairCorrelationPigs::wave_t * wave)
{
  
  for(int i=iTimeSliceBegin;i<iTimeSliceEnd;i++)
    {
      wave->pair_correlation_symm(&((*w)[i]),this->ms,this->ms->grid->b,setA);
    }
    
  
}


void energyThermodynamic::make_measurement(walker_t* w,wave_t* wave)
{
  this->ms->add(- prop->logEvaluateTimeDerivative(*w)/w->size(),0);
}
