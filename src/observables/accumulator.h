#include <vector>

using namespace std;

template<class T >
class vectorAccumulator
{
public:
  
  vectorAccumulator(int n){vectorAccumulatorMean.resize(n,0);nMeasurements=0;}

  
  void accumulateMean(const vector<T> & vectorName);
  virtual void accumulate(const vector<T> & vectorName){accumulateMean(vectorName);}
  void getMean(vector<T> &vecOut) const;
  void incrementCounter(){nMeasurements+=1;};
  virtual void reset();
  void decreaseCounter(){nMeasurements-=1;};
  int getCounter() const {return nMeasurements;}
  virtual void transfer(int root); // sum over all tasks and store result in root. Other tasks are reset.
private:
  int nMeasurements;
  vector<T> vectorAccumulatorMean;
  
};

template<class T>
class vectorAccumulatorVariance : public vectorAccumulator<T>
{
public:
  void getMeanSquares(vector<T> &vecOut) const;
  void accumulateMeanSquares(const vector<T> &vecIn);
  virtual void accumulate(const vector<T> &vecIn){this->accumulateMean(vecIn);accumulateMeanSquares(vecIn);}
  
  vectorAccumulatorVariance(int n):vectorAccumulator<T>(n){nMeasuresSquares=0;vectorAccumulatorSquares.resize(n,0.);}
  
  void getVariances(vector<T> &out) const;
  virtual void reset();
  virtual void transfer(int root);
  
private:
  int nMeasuresSquares;
  vector<T> vectorAccumulatorSquares;
  
};
