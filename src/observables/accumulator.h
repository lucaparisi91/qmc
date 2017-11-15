#ifndef ACCUMULATOR_H
#define ACCUMULATOR_H

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
  
  int getSize() const {return vectorAccumulatorMean.size();};
  virtual void accumulate(const vector<T> & vectorName){accumulateMean(vectorName);}
  
  virtual void accumulate(const vector<T> & vectorName,const vector<T> & weightsIn){accumulateMean(vectorName,weightsIn);}
  
  void getMean(vector<T> &vecOut) const;
  virtual void reset();
  
  virtual void transfer(int root); // sum over all tasks and store result in root. Other tasks are reset.

  T getTotWeight() const;
  int getNmeasurements() const {return nMeasurements;};
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
  
  virtual void accumulate(const vector<T> &vecIn){this->accumulateMean(vecIn);accumulateMeanSquares(vecIn);}
  
  virtual void accumulate(const vector<T> &vecIn,const vector<T> &weightsIn){this->accumulateMean(vecIn,weightsIn);accumulateMeanSquares(vecIn,weightsIn);}
  
  vectorAccumulatorVariance(int n):vectorAccumulatorType(n),accumulatorSquares(n){resize(n);}
  
  vectorAccumulatorVariance():vectorAccumulatorType(),accumulatorSquares(){}
  
  void resize(int n){vectorAccumulatorType::resize(n);accumulatorSquares.resize(n);scratch.resize(n);}
  
  void getVariances(vector<T> &out) const;
  
  virtual void reset(){vectorAccumulatorType::reset();accumulatorSquares.reset();};
  
  virtual void transfer(int root){vectorAccumulatorType::transfer(root);accumulatorSquares.transfer(root);};
  void getMeanError(vector<T> &mean,vector<T> &out) const;
private:
  vectorAccumulatorType accumulatorSquares;
  vector<T> scratch;
};

#endif
