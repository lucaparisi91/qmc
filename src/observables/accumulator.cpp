#include "accumulator.h"
#include "cassert"
#include "../ptools.h"

template<class T>
void vectorAccumulator<T>::accumulateMean(const vector<T> & vectorName)
{
  assert(vectorName.size()==size());
  for(int i=0;i<vectorName.size();i++)
    {
      vectorAccumulatorMean[i]+=vectorName[i];
      weights[i]+=1;
    }
  nMeasurements+=1;
}

template<class T>
void vectorAccumulator<T>::accumulateMean(const vector<T> & vectorName,const vector<T> &weightsIn)
{
  
  assert(vectorName.size()==size());
  assert(weightsIn.size()==size());
  assert(weights.size()==size());
  
  for(int i=0;i<vectorName.size();i++)
    {
      vectorAccumulatorMean[i]+=vectorName[i]*weightsIn[i];
      weights[i]+=weightsIn[i];
    }
  nMeasurements+=1;
}

template<class T>
void vectorAccumulator<T>::getMean(vector<T> & vectorOut) const
{
  vectorOut.resize(this->size());
  
  for(int i=0;i<vectorOut.size();i++)
    {
      assert(weights[i]>0);
      vectorOut[i]=vectorAccumulatorMean[i]/weights[i];
    }
}

template<class T>
void vectorAccumulator<T>::transfer(int root)
{
  pTools::transferSum(vectorAccumulatorMean,root);
  pTools::transferSum(weights,root);
  pTools::transferSum(nMeasurements,root); 
}

template<class T>
void vectorAccumulator<T>::reset()
{
  for(int i=0;i<vectorAccumulatorMean.size();i++)
    { 
      vectorAccumulatorMean[i]=0;
      weights[i]=0;
    }
  nMeasurements=0;
}

template<class T>
T vectorAccumulator<T>::getTotWeight() const
{
  T w;
  w=0;
  for(int i=0;i<weights.size();i++)
    { 
      w+=weights[i];
    }
  return w;
}


template<class T>
void  vectorAccumulatorVariance<T>::getVariances(vector<T> & vectorOut) const
{
  vector<T> tmp;
  tmp.resize(this->size());
  vectorOut.resize(this->size());
  this->getMean(vectorOut);
  getMeanSquares(tmp);
  for(int i=0;i<vectorOut.size();i++)
    {
      vectorOut[i]=tmp[i]-pow(vectorOut[i],2);
    }
}

template<class T>
void  vectorAccumulatorVariance<T>::accumulateMeanSquares(const vector<T> & vectorName)
{
  assert(vectorName.size()==this->size() );
  scratch.resize(this->size());
  
  for(int i=0;i<this->size();i++)
    {
      scratch[i]=vectorName[i]*vectorName[i];
    }
  
  accumulatorSquares.accumulate(scratch);
  
}

template<class T>
void  vectorAccumulatorVariance<T>::accumulateMeanSquares(const vector<T> & vectorName,const vector<T> & weightsIn)
{
  assert(vectorName.size()==this->size() );
  scratch.resize(this->size());
  
  for(int i=0;i<this->size();i++)
    {
      scratch[i]=vectorName[i]*vectorName[i];
    }
  
  accumulatorSquares.accumulate(scratch,weightsIn);
  
}

template<class T>
void  vectorAccumulatorVariance<T>::getMeanError(vector<T> & meanOut,vector<T> & errorOut) const
{
  assert(this->getNmeasurements()>0);
  this->getMean(meanOut);
  getMeanSquares(errorOut);
  for(int i=0;i<this->size();i++)
    {
      errorOut[i]=sqrt(abs(errorOut[i]-meanOut[i]))/this->getNmeasurements();
    }
  
}

template class vectorAccumulator<double>;
template class vectorAccumulatorVariance<double>;
