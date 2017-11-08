#include "accumulator.h"
#include "cassert"
#include "../ptools.h"


template<class T>
void vectorAccumulator<T>::accumulateMean(const vector<T> & vectorName)
{
  for(int i=0;i<vectorName.size();i++)
    {
      vectorAccumulatorMean[i]+=vectorName[i];
    }
  
  incrementCounter();
}

template<class T>
void vectorAccumulator<T>::getMean(vector<T> & vectorOut) const
{
  assert(getCounter()>0);
  
  for(int i=0;i<vectorOut.size();i++)
    {
      vectorOut[i]=vectorAccumulatorMean[i]/getCounter();
    }
}

template<class T>
void vectorAccumulator<T>::transfer(int root)
{
  pTools::transferSum(vectorAccumulatorMean,root);
  pTools::transferSum(nMeasurements,root); 
}

template<class T>
void vectorAccumulator<T>::reset()
{
  nMeasurements=0;
  for(int i=0;i<vectorAccumulatorMean.size();i++)
    {
      
      vectorAccumulatorMean[i]=0;
    }
}

template<class T>
 void vectorAccumulatorVariance<T>::reset()
{
  vectorAccumulator<T>::reset();
  nMeasuresSquares=0;
  for(int i=0;i<vectorAccumulatorSquares.size();i++)
    {
      vectorAccumulatorSquares[i]=0;
    }
}

template<class T>
void vectorAccumulatorVariance<T>::transfer(int root)
{
  
  vectorAccumulator<T>::transfer(root);
  
  pTools::transferSum(vectorAccumulatorSquares,root);
  pTools::transferSum(nMeasuresSquares,root);
  
}

template<class T>
void  vectorAccumulatorVariance<T>::getMeanSquares(vector<T> & vectorOut) const
{
  assert(nMeasuresSquares>0);
  
  for(int i=0;i<vectorOut.size();i++)
    {
      vectorOut[i]=vectorAccumulatorSquares[i]/nMeasuresSquares;
    }
}

template<class T>
void  vectorAccumulatorVariance<T>::accumulateMeanSquares(const vector<T> & vectorName)
{
  assert(vectorName.size()==vectorAccumulatorSquares.size() );
    
  for(int i=0;i<vectorName.size();i++)
    {
      vectorAccumulatorSquares[i]+=vectorName[i]*vectorName[i];
    }
  
  nMeasuresSquares+=1;
  
}



template class vectorAccumulator<double>;
template class vectorAccumulatorVariance<double>;
