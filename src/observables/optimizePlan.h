#ifndef OPTIMIZEPLAN_H
#define OPTIMIZEPLAN_H

#include <vector>
#include <map>
#include <utility>
#include <string>
using namespace std;

class optimizePlan
{
public:
  typedef vector<pair<int,int> > paramRegisterType;
  typedef map<string,paramRegisterType> registerType;
  typedef registerType::iterator iterator;
  typedef registerType::const_iterator const_iterator;
  
  iterator begin(){return parameterRegister.begin();}
  const_iterator begin() const {return parameterRegister.begin();}
  iterator end(){return parameterRegister.end();}
  const_iterator end() const {return parameterRegister.end();}
  
  
  
  void add(string key,int waveKey,int paramKey)
  {
    parameterRegister[key].push_back(pair<int,int>(waveKey,paramKey));
  }
  
  paramRegisterType& operator[](string key){return parameterRegister[key];}
  
  void print();
  void setDelta(double delta_){delta=delta_;}
  
  double getDelta() const {return delta;};
  
  
private:
  double delta;
  registerType parameterRegister;
};

optimizePlan buildOptimizePlan(string filename);

#endif
