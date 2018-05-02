#ifndef PIGS_CONFIGURATIONS_H
#define PIGS_CONFIGURATIONS_H
#include "pigsTraits.h"
#include <iostream>

class configurationsPIGS
{
public:
  
  typedef system_t::all_particles_t configuration_t;
  typedef system_t::pos_t pos_t;
  inline configuration_t & operator[](int i){return configurations[i];}
  inline const configuration_t &  operator[](int i) const {return configurations[i];}
  
  void push_back(configuration_t & config)
  {
    configurations.push_back(config);
  }
  
  inline size_t size() {return configurations.size();}
  
  void resize(int n){configurations.resize(n);}
  
  double maxTime(){return deltaTau*configurations.size();}
  
  void setTimeStep(double deltaTau_){deltaTau=deltaTau_;}
  
  double getTimeStep(){return deltaTau;}
  
  friend ostream& operator<<(ostream& out,configurationsPIGS & configurations);
  
private:
  vector<configuration_t> configurations;
  double deltaTau;  
};

class distancesTable
{  
public:
  double & distance(unsiged int i,unsigned int j,unsigned int k);
  double & distance(unsigned int i,unsigned int k);
private:
  vector<double> singleParticleDistances;
  vector<vector<double> > pairParticleDistances;
  
}
#endif




