#ifndef PHASE_H
#define PHASE_H


#include<vector>
#include<cmath>

template<class comp>
class phase
{
public:
  typedef typename comp::geometry_t geometry_t;
  
  phase(geometry_t  * geo_) {geo=geo_;};
  typedef typename comp::all_particles_t all_particles_t;
  virtual double evaluate(all_particles_t* p){return 0;};
  
protected:
  
  geometry_t * geo;
  
};

  
#endif



