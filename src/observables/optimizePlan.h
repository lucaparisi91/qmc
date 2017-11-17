#ifndef OPTIMIZEPLAN_H
#define OPTIMIZEPLAN_H

#include <vector>
#include <map>
#include <utility>
#include <string>
using namespace std;

class optimizationParameterType
{
public:
  
  optimizationParameterType(){boundedBelow=false;boundedAbove=false;};
  
  int size() const {return waveParameters.size();}
  pair<int,int>& operator[](size_t i){return waveParameters[i];}
  const pair<int,int>& operator[](size_t i) const {return waveParameters[i];}
  
  void push_back(pair<int,int> value){waveParameters.push_back(value);}
  
  void setMinParameter(double minParameter_){minParameter=minParameter_;boundedBelow=true;}
  
  void setMaxParameter(double maxParameter_){maxParameter=maxParameter_;boundedAbove=true;}
  
  double getMinParameter() const {return minParameter;}
  double getMaxParameter() const {return maxParameter;}
  
  bool isBoundedBelow() const {return boundedBelow;}
  bool isBoundedAbove() const {return boundedAbove;}
  
private:
  vector<pair<int,int> > waveParameters;
  double minParameter;
  double maxParameter;
  bool boundedBelow;
  bool boundedAbove;
};

class optimizePlan
{
public:
  typedef optimizationParameterType paramRegisterType;
  typedef map<string,paramRegisterType> registerType;
  typedef registerType::iterator iterator;
  typedef registerType::const_iterator const_iterator;
  
  iterator begin(){return parameterRegister.begin();}
  const_iterator begin() const {return parameterRegister.begin();}
  iterator end(){return parameterRegister.end();}
  const_iterator end() const {return parameterRegister.end();}
  
  bool checkBoundsParameters(const vector<double> &parameters) const;
  
  void add(string key,int waveKey,int paramKey)
  {
    parameterRegister[key].push_back(pair<int,int>(waveKey,paramKey));
  }
  
  paramRegisterType& operator[](string key){return parameterRegister[key];}
  
  void print();
  void setDelta(double delta_){delta=delta_;}
  
  double getDelta() const {return delta;};
  
  void setMaxParameter(string key,double maxValue){parameterRegister[key].setMaxParameter(maxValue);}

  void setMinParameter(string key,double minValue){parameterRegister[key].setMinParameter(minValue);}
  
private:
  double delta;
  registerType parameterRegister;
};

optimizePlan buildOptimizePlan(string filename);

#endif
