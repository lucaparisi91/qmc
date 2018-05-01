#ifndef PIGSTOOLS_H
#define PIGSTOOLS_H

#include "../xml-input.h"
#include <string>

bool logMetropolis(double log_ratio,system_t::rand_t* rand);

string toString(double value);

string getName(xml_input * main_input);

string getAttribute(xml_input * main_input,string name,string defaultValue);

class towerSampling
{
public:
  towerSampling();
  int sample();
  void add(double p);
  void setRandomGenerator(system_t::rand_t * rand2){randO=rand2;}
  
private:
  system_t::rand_t * randO;
  vector<double> cumulative;
  double pTot;
};


#endif
