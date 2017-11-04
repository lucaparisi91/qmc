#ifndef POTENTIAL_H
#define POTENTIAL_H

#include<vector>
#include<cmath>
#include<string>
#include<iostream>
#include "traits.h"
#include "system.h"
#include "xml-input.h"

using namespace std;
template<class comp>

class potential
{
public:
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename comp::geometry_t geometry_t;
  
  potential(geometry_t  * geo_)  {geo=geo_;};
  virtual double evaluate(all_particles_t* p){return 0;};
protected:
  geometry_t * geo;
  
};

template<class comp>
class empty_potential : public potential<comp>
{
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename comp::geometry_t geometry_t;
  
public:
  empty_potential(geometry_t* geo_) : potential<comp>(geo_){};
};
  
//the class with a certain degree of Rabi Coupling
class rabiCoupling
{
public:
  rabiCoupling(double omega_){omega=omega_;}
  inline double get_omega(){return omega;};
private:
  double omega;
  
};

class harmonicPotential
{
public:
  harmonicPotential(const double &omega_,const double &x0_,int set_=0){omega=omega_;x0=x0_;set=set_;};
  
  template<class all_particles_t,class geometry_t>
  double evaluate(geometry_t* geo,all_particles_t * p)
  {
    typedef typename all_particles_t::particles_t particles_t;
    particles_t* p1;
    int i;
    double value;
    value=0;
    p1=p->particle_sets[set];
   
    for(i=0;i<p1->n;i++)
      {
	value+=pow(geo->distance_pbc(p1->position[i],x0),2);
      }
    
    return value*1./2*omega*omega;
  }
  
private:
  double omega;
  double x0;
  int set;
};

template<class comp>
class speciesHarmonicPotential : public potential<comp>
{
public:
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename comp::geometry_t geometry_t;
  speciesHarmonicPotential(geometry_t * geo,vector<double> &omegas,vector<double> &x0,vector<int> &sets);
  double evaluate(all_particles_t* p);
  
  vector<harmonicPotential *> traps;
};


#include "potential.hpp"

#endif
