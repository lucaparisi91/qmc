#include "ptools.h"
template<class T>
class storageObject
{
public:
  
  typedef T scalar_t;
  virtual void pack(packed_data* md_pack)
  {
    md_pack->pack(n);
    md_pack->pack(m);
    md_pack->pack(n2);
    md_pack->pack(m2);
  }

  storageObject()
  {
    n=0;
    n2=0;
    filled=false;
  }
  
  virtual void unpack(packed_data* md_pack)
  {
    md_pack->unpack(n);
    md_pack->unpack(m);
    md_pack->unpack(n2);
    md_pack->unpack(m2);
  }
  
  virtual void clone(storageObject<T> & s)
  {
    n=s.n;
    m=s.m;
    n2=s.n2;
    m2=s.m2;
    
  }
  int get_pack_size()
  {
    return pTools::get_pack_size(n) + pTools::get_pack_size(m) + pTools::get_pack_size(n2) + pTools::get_pack_size(m2);
  }
  
  void clear()
  {
    n=0;
    
    n2=0;
    
  }
  
  void reset()
  {
    
    filled=true;
    
    n2=n;
    
    n=0;
  }
  
  inline scalar_t & get_sum(){return m;};
  inline int & get_n(){return n;};
  
  
  bool isFilled()
  {
    if (n2>0)
      {
	return true;
      }
    else
      {
	return false;
      }
  }
  
  void increment_index()
  {
    
    n=n+1;
  }
  
protected:
  scalar_t m;
  int n;
  scalar_t m2;
  int n2;
  bool filled;
};

template<class T>
class scalarStorage : public storageObject<T>
{
public:
  typedef T scalar_t;
  
  scalarStorage() : storageObject<scalar_t >()
  {
    this->m=0;
    this->m2=0;
  }
  
  scalarStorage<scalar_t> & operator+=(scalar_t m2)
  {
    
    this->m+=m2;
    this->n+=1;
  }
  
  void reset()
  {
    storageObject<scalar_t >::reset();
    this->m2=this->m;
    this->m=0;
  }
  
  void clear()
  {
    storageObject<scalar_t>::clear();
    this->m=0;
    this->m2=0;
  }
  
  scalar_t average()
  {
    return this->m2/this->n2;
  }
};

template<class T>
class vectorStorage : public storageObject<vector<T> >
{
public:
  typedef vector<T> scalar_t;
  vectorStorage(int ns) : storageObject<scalar_t >()
  {
    
    this->m.resize(ns);
    this->m2.resize(ns);
  }
  vectorStorage<scalar_t> & operator+=(scalar_t & m2)
  {
    
    assert(this->m.size()==m2.size() );
    for(int i=0;i<this->m.size();i++)
      {
	this->m[i]+=m2[i];
      }
    this->n+=1;
  }
  
  void increment_value(T x,int j,double step)
  {
    
    assert(j<this->m.size());
    assert(j>=0);
    this->m[j]+=x;
  }

  // reset the value
  void reset()
  {
    storageObject<scalar_t >::reset();
    this->m2=this->m;
    
    for(int i=0;i<this->m2.size();i++)
      {
	this->m[i]=0;
      }
    
  }

  void clear()
  {
    storageObject<scalar_t >::clear();
    for(int i=0;i<this->m2.size();i++)
      {
	this->m2[i]=0;
      }
    for(int i=0;i<this->m.size();i++)
      {
	this->m[i]=0;
      }
  }

  void average(vector<T> & v)
  {
    
    assert(v.size()==(this->m2).size());
    for(int i=0;i<v.size();i++)
      {
	
	v[i]=this->m2[i]/this->n2;
      }
    
  }
  
};

// storage when a space vector is needed
template<class T,class mesh_t>
class vectorSpaceStorage : public vectorStorage<T>
{
public:
  typedef vector<T> scalar_t;
  vectorSpaceStorage(int ns,mesh_t m_) : vectorStorage<T>::vectorStorage(ns),m(m_){}

  void increment_value(T value,T x,double step)
  {
    vectorStorage<T>::increment_value(value,m.index(x),step);
  }

  double get_step()
  {
    return m.get_step();
  }
  
private:
  mesh_t m;
  
};
