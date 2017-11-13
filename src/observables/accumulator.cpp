#include "accumulator.h"
#include "cassert"
#include "../ptools.h"

template<class T>
void vectorAccumulator<T>::accumulateMean(const vector<T> & vectorName)
{
  
  assert(vectorName.size()==getSize());
  for(int i=0;i<vectorName.size();i++)
    {
      vectorAccumulatorMean[i]+=vectorName[i];
      weights[i]+=1;
    }
}

template<class T>
void vectorAccumulator<T>::accumulateMean(const vector<T> & vectorName,const vector<T> &weightsIn)
{
  
  assert(vectorName.size()==getSize());
  assert(weightsIn.size()==getSize());
  assert(weights.size()==getSize());
  
  for(int i=0;i<vectorName.size();i++)
    {
      vectorAccumulatorMean[i]+=vectorName[i]*weightsIn[i];
      weights[i]+=weightsIn[i];
    }
}

template<class T>
void vectorAccumulator<T>::getMean(vector<T> & vectorOut) const
{
  T w;
  w=getTotWeight();
  assert(w>0);
  vectorOut.resize(this->getSize());
  
  for(int i=0;i<vectorOut.size();i++)
    {
      vectorOut[i]=vectorAccumulatorMean[i]/weights[i];
    }
}

template<class T>
void vectorAccumulator<T>::transfer(int root)
{
  pTools::transferSum(vectorAccumulatorMean,root);
  pTools::transferSum(weights,root); 
}

template<class T>
void vectorAccumulator<T>::reset()
{
  for(int i=0;i<vectorAccumulatorMean.size();i++)
    { 
      vectorAccumulatorMean[i]=0;
      weights[i]=0;
    }
}

template<class T>
T vectorAccumulator<T>::getTotWeight() const
{
  T w;
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
  tmp.resize(this->getSize());
  vectorOut.resize(this->getSize());
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
  assert(vectorName.size()==this->getSize() );
  scratch.resize(this->getSize());
  
  for(int i=0;i<this->getSize();i++)
    {
      scratch[i]=vectorName[i]*vectorName[i];
    }
  
  accumulatorSquares.accumulate(scratch);
  
}

template<class T>
void  vectorAccumulatorVariance<T>::accumulateMeanSquares(const vector<T> & vectorName,const vector<T> & weightsIn)
{
  assert(vectorName.size()==this->getSize() );
  scratch.resize(this->getSize());
  
  for(int i=0;i<this->getSize();i++)
    {
      scratch[i]=vectorName[i]*vectorName[i];
    }
  
  accumulatorSquares.accumulate(scratch,weightsIn);
  
}


template class vectorAccumulator<double>;
template class vectorAccumulatorVariance<double>;
