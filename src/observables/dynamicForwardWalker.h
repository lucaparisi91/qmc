#ifndef DYNAMICFORWARDWALKER_H
#define DYNAMICFORWARDWALKER_H
#include "../measures.h"
#include "accumulator.h"

template<class qt>
class scalarForwardWalking : public measure_dynamic<qt>
{
public:
  typedef typename qt::walker_t walker_t;
  typedef typename qt::wave_t wave_t;
  
  scalarForwardWalking(){currentIndex=0;}
  scalarForwardWalking(int bins){currentIndex=0;accumulator.resize(bins);}
  void resize(int bins){accumulator.resize(bins);}
  
  virtual void add(double m)
  {
    
    accumulator.reset(currentIndex);
    accumulator.accumulate(currentIndex,m);
    currentIndex=(currentIndex+1)%accumulator.size();
   }
  
  virtual double average()
  {
    double x;
    accumulator.getMean(currentIndex,x);
    return x;
  }

  virtual bool isFilled()
  {
    return accumulator.isFilled(currentIndex);
  }
  
  virtual void reset(){};
  virtual void print(){};
  
  virtual void pack(packed_data* packedO){packedO->pack(currentIndex);accumulator.pack(packedO);}

  virtual void unpack(packed_data* packedO){packedO->unpack(currentIndex);accumulator.unpack(packedO);}
  
  virtual int get_pack_size(){return pTools::get_pack_size(currentIndex) + accumulator.get_pack_size(); }

  
  measure_dynamic<qt> & operator=(measure_dynamic<qt> &m)
  {
    scalarForwardWalking<qt> * m2;
    m2=static_cast<scalarForwardWalking<qt> *>(&m);
    currentIndex=m2->currentIndex;
    accumulator=m2->accumulator;
  }
  
private:
  
  int currentIndex;
  vectorAccumulator<double> accumulator;
};

#endif
