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
  typedef typename comp::wave_t wave_t;
  
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
template<class comp>
class rabiCouplingPotential : public potential<comp>
{
public:
  typedef typename potential<comp>::geometry_t geometry_t;
  typedef typename potential<comp>::all_particles_t all_particles_t;
 typedef typename comp::wave_t wave_t;
  
  rabiCouplingPotential(double omega_,geometry_t &geo_,wave_t& wave_) : potential<comp>(&geo_){setOmega(omega_);wave=&wave_;}
  
  void setOmega(double omega_){omega=omega_;}
  double getOmega() const{return omega;};
  
  double evaluate(all_particles_t *state)
  {
    double ev=0;
    
    for(int i=0;i<state[0].size();i++)
      {
	ev+=wave->spinFlip(i,*state);
      }
    
    return ev*omega;
  }
  
private:
  double omega;
  wave_t* wave;
  
};

class harmonicPotential
{
public:
  harmonicPotential(const double &omega_,const double &x0_,int set_=0){omega=omega_;x0=x0_;set=set_;};
  
  template<class all_particles_t,class geometry_t>
  double evaluate(geometry_t* geo,all_particles_t * p)
  {
    typedef typename all_particles_t::particles_t particles_t;
    
    int i;
    double value;
    value=0;

    particles_t& p1=(*p)[set];
   
   
    for(i=0;i<p1.size();i++)
      {
	value+=pow(geo->distance_pbc(p1[i].position(),x0),2);
      }
    
    return value*1./2*omega*omega;
  }
  
private:
  double omega;
  double x0;
  int set;
};


class boxPotentialManyBodySymmetric
{
  
public:
  boxPotentialManyBodySymmetric(){set=0;};

  void setParameter(string parameter,double value)
  {
    if (parameter=="V0")
      {
	V0=value;
      }
    else if (parameter=="L")
      {
	halfLength=value/2.;
      }
    else if (parameter=="set")
      {
	set=(int)value;
      }
    else
      {
	printf("Unkown parameter.");
	exit(2);
      }
  }

  double getParameter(string parameter)
  {
   
    if (parameter=="V0")
      {
	return V0;
      }
    else if (parameter=="L")
      {
	return halfLength*2.;
	
      }
    else if (parameter=="set")
      {
	return set;
      }
    else
      {
	printf("Unkown parameter.");
	exit(2);
      }
  }
  
  template<class all_particles_t,class geometry_t>
  double evaluate(geometry_t* geo,all_particles_t * p)
  {
    typedef typename all_particles_t::particles_t particles_t;
    
    int i;
    double value,x,d;
    value=0;

    particles_t& p1=(*p)[set];

     
    for(int i=0;i<p1.size();i++)
      {
	for(int j=0;j<i;j++)
	  {
	    x=geo->distance_pbc(p1[i].position(),p1[j].position());
	    d=abs(x);
	    if (d<halfLength)
	      {
		value+=V0;
	      }
	  }
      }
  
  return value;
  }
  
private:
  
  double halfLength;
  double V0;
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


template<class comp,class potential_t>
class totalPotential : public potential<comp>
{
public:
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename comp::geometry_t geometry_t;
  
  void setParameter(int i,string label,double value)
  {
    potentials[i].setParameter(label,value);
  }
  
  double getParameter(int i,string label)
  {
    return potentials[i].getParameter(label);
  }

  void addPotential()
  {
    potential_t toAdd;
    potentials.push_back(toAdd);
  }
  
private:
  
 vector<potential_t> potentials;
};



#include "potential.hpp"

#endif
