#ifndef DYNAMICFORWARDWALKER_H
#define DYNAMICFORWARDWALKER_H
#include "../measures.h"
#include "accumulator.h"
#include "../tools.h"

template<class qt,class storage_t>
class forwardWalking : public measure_dynamic<qt>
{
public:
  typedef typename qt::walker_t walker_t;
  typedef typename qt::wave_t wave_t;
  
  forwardWalking(){currentIndex=0;}
  
  forwardWalking(int bins,const storage_t & initValue)
  {
    currentIndex=0;
    data.resize(bins,initValue);
    filled=0;
  }
  
  void resize(int bins,const storage_t & initValue){data.resize(bins,initValue);}
  
  virtual void add(const storage_t &m)
  {
    data[currentIndex]=m;
    currentIndex=(currentIndex+1)%data.size();
    if (currentIndex==0)
      {
	filled=1;
      }
    
  }
  
  virtual void incrementIndex()
  {
    currentIndex=(currentIndex+1)%data.size();
    if (currentIndex==0){filled=1;}
  }
  
  storage_t & current()
  {
    return data[currentIndex];
  }
  
  const storage_t & current() const
  {
    return data[currentIndex];
  }
  
  virtual int isFilled()
  {
    return filled;
  }
  
  virtual void print(){};
  
  virtual void reset()
  {
    currentIndex=0;
    filled=0;
  };
  
  
  virtual void pack(packed_data* packedO){packedO->pack(currentIndex);packedO->pack(data);packedO->pack(filled);}

  virtual void unpack(packed_data* packedO){packedO->unpack(currentIndex);packedO->unpack(data);packedO->unpack(filled);}
  
  virtual int get_pack_size(){return pTools::get_pack_size(currentIndex) + pTools::get_pack_size(data) +  pTools::get_pack_size(filled); }

  
  measure_dynamic<qt> & operator=(measure_dynamic<qt> &m)
  {
    forwardWalking<qt,storage_t> * m2;
    m2=static_cast< forwardWalking<qt,storage_t> *>(&m);
    currentIndex=m2->currentIndex;
    filled=m2->filled;
    data=m2->data;
  }

  vector<storage_t> & getData(){return data;}
  const vector<storage_t> & getData() const {return data;} ;
private:
  
  int filled;
  int currentIndex;
  vector<storage_t> data;
};

template<class qt>
class scalarForwardWalking : public forwardWalking<qt,double>
{
public:
  scalarForwardWalking() : forwardWalking<qt,double>() {};
  scalarForwardWalking(int bins) : forwardWalking<qt,double>(bins,0) {};
  
  scalarForwardWalking(int bins,double & m) : forwardWalking<qt,double>(bins,m) {};
  
  double average(){return this->current();}
  
  virtual vector<double> & currentVector(){return this->getData();}
  
private:
  
};

template<class qt>
class vectorForwardWalking : public forwardWalking<qt,vector<double> >
{
public:
  vectorForwardWalking() : forwardWalking<qt,vector<double> >() {};
  
  vectorForwardWalking(int bins,int binsVector) : forwardWalking<qt,vector<double> >(bins,vector<double>(binsVector,0)) {};
  
  virtual vector<double> & currentVector(){return this->current();}
private:  
  
};

#endif
