#ifndef ACCUMULATOR_H
#define ACCUMULATOR_H
#include "../ptools.h"
#include <vector>

using namespace std;
  
template<class T >
class vectorAccumulator
{
public:
  vectorAccumulator(int n){resize(n);}
  void resize(int n){vectorAccumulatorMean.resize(n,0);weights.resize(n,0);}

  
  vectorAccumulator(){nMeasurements=0;}
  
  void accumulateMean(const vector<T> & vectorName);
  void accumulateMean(const vector<T> & vectorName,const vector<T> & weightsIn);

  void accumulateMean(int i,const T& value){vectorAccumulatorMean[i]+=value;weights[i]+=1;};
  void accumulateMean(int i,const T& value, const T& weight){vectorAccumulatorMean[i]+=value;weights[i]+=weight;};

  virtual void accumulate(int i,const T& value, const T& weight){accumulateMean(i,value,weight);};

  virtual void accumulate(int i,const T& value){accumulateMean(i,value);};

  
  int size() const {return vectorAccumulatorMean.size();};
  virtual void accumulate(const vector<T> & vectorName){accumulateMean(vectorName);}
  
  virtual void accumulate(const vector<T> & vectorName,const vector<T> & weightsIn){accumulateMean(vectorName,weightsIn);}
  
  void getMean(vector<T> &vecOut) const;

  
  virtual void reset();
  
  void reset(int i){weights[i]=0;vectorAccumulatorMean[i]=0;}
  void getMean(int i,T& value) const {value=vectorAccumulatorMean[i]/weights[i];} 
  bool isFilled(int i){return weights[i]==0 ? false : true; }
  virtual void transfer(int root); // sum over all tasks and store result in root. Other tasks are reset.

  T getTotWeight() const;
  int getNmeasurements() const {return nMeasurements;};
  void getMean(T& value,int i){value=vectorAccumulatorMean[i]/weights[i];}
  
  int get_pack_size()
  {
    int s;
    s= pTools::get_pack_size(nMeasurements)+ pTools::get_pack_size(weights) +pTools::get_pack_size(vectorAccumulatorMean);
    return s;
  }
  
  void pack(packed_data* packO)
  {
    packO->pack(nMeasurements);
    packO->pack(weights);
    packO->pack(vectorAccumulatorMean);
    
  }

  void unpack(packed_data* packO)
  {
    packO->unpack(nMeasurements);
    packO->unpack(weights);
    packO->unpack(vectorAccumulatorMean);
    
  }
  
private:
  int nMeasurements;
  vector<T> weights;
  vector<T> vectorAccumulatorMean;  
};

template<class T>
class vectorAccumulatorVariance : public vectorAccumulator<T>
{
public:
  
  typedef  vectorAccumulator<T> vectorAccumulatorType;
  
  void getMeanSquares(vector<T> &vecOut) const {return accumulatorSquares.getMean(vecOut);};
  
  void accumulateMeanSquares(const vector<T> &vecIn);
  void accumulateMeanSquares(const vector<T> &vecIn,const vector<T> &weightsIn);

  void accumulateMeanSquares(int i,const T& value){accumulatorSquares.accumulateMean(i,value);};
  void accumulateMeanSquares(int i,const T& value, const T& weight){accumulatorSquares.accumulateMean(i,value,weight);};
  
  virtual void accumulate(int i,const T& value, const T& weight){this->accumulateMean(i,value,weight);accumulateMeanSquares(i,value,weight);};

  virtual void accumulate(int i,const T& value){this->accumulateMean(i,value);accumulateMeanSquares(i,value);};

  virtual void accumulate(const vector<T> &vecIn){this->accumulateMean(vecIn);accumulateMeanSquares(vecIn);}
  
  virtual void accumulate(const vector<T> &vecIn,const vector<T> &weightsIn){this->accumulateMean(vecIn,weightsIn);accumulateMeanSquares(vecIn,weightsIn);}
  
  vectorAccumulatorVariance(int n):vectorAccumulatorType(n),accumulatorSquares(n){resize(n);}
  
  vectorAccumulatorVariance():vectorAccumulatorType(),accumulatorSquares(){}
  
  void resize(int n){vectorAccumulatorType::resize(n);accumulatorSquares.resize(n);scratch.resize(n);}
  
  void getVariances(vector<T> &out) const;
  
  virtual void reset(){vectorAccumulatorType::reset();accumulatorSquares.reset();};

  virtual void reset(int i){vectorAccumulatorType::reset(i);accumulatorSquares.reset(i);};
  
  virtual void transfer(int root){vectorAccumulatorType::transfer(root);accumulatorSquares.transfer(root);};
  void getMeanError(vector<T> &mean,vector<T> &out) const;
private:
  vectorAccumulatorType accumulatorSquares;
  vector<T> scratch;
};

#endif
