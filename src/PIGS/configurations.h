#ifndef PIGS_CONFIGURATIONS_H
#define PIGS_CONFIGURATIONS_H

#include "../dmc.h"

template<class comp>
class pigsWalker_t : public walker<comp>
{
public:

  void print();

  template<class qmc_t>
  pigsWalker_t(qmc_t* qmc_obj) : walker<comp>(qmc_obj) {};
  
};

template<class T>
class configurationsPIGS_t
{
  
  
public:
  typedef T configuration_t;
  
  inline configuration_t & operator[](int i){return configurations[i];}
  inline const configuration_t &  operator[](int i) const {return configurations[i];}
  
  void push_back(configuration_t & config)
  {
    configurations.push_back(config);
  }

  inline size_t size() {return configurations.size();}

  void resize(int n){configurations.resize(n);}
  
  void maxTime(){return deltaTau*configurations.size();}
  
  void setTimeStep(double deltaTau_){deltaTau=deltaTau_;}
  
  double getTimeStep(){return deltaTau;}
  
private:
  vector<configuration_t> configurations;
  double deltaTau;
  
};


template<class comp>
class buildConfigurationsPIGS
{
public:
  typedef comp qmc_t;
  typedef typename comp::pos_t pos_t;
  typedef typename comp::configurations_t configurations_t;
  typedef typename configurations_t::configuration_t configuration_t;
  
  buildConfigurationsPIGS(){};
  
  void initConfigurations(qmc_t * qmcO)
  {
    int nBeads=0;
    
    string kind="random";
    pos_t boxLength;
    
    boxLength=qmcO->geo->l_box;
    
    nBeads=qmcO->main_input->reset()->get_child("method")->get_child("nBeads")->get_value()->get_int();

    nBeads=qmcO->main_input->reset()->get_child("method")->get_child("nBeads")->get_value()->get_int();
    
    configurations_t & configurations=qmcO->getConfigurations();
    configurations.resize(nBeads);
    
    for(int i=0;i<nBeads;i++)
      {
	
	buildAllOrbitals("input.xml",configurations[i]);
	generateRandom(configurations[i],-boxLength/2.,boxLength/2.,qmcO->rand);
	qmcO->geo->all_particles_pbc(configurations[i]);
      }
    
  }

};


#endif
