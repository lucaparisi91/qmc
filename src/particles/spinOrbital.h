#ifndef SPINORBITAL_H
#define SPINORBITAL_H

#include <iostream>
#include <vector>
#include "../ptools.h"
#include "../tools.h"

using namespace std;
class spinOrbital
{
public:
  spinOrbital(){spinValue=0;x=0;xBC=0;}
  inline double& position(){return x;};
  inline double& positionBC(){return xBC;};
  inline const double& positionBC() const{return xBC;};
  const double& position() const {return x;};
  int& spin(){return spinValue;};
  const int& spin() const {return spinValue;};
  int getPackSize() const;
  void pack(packed_data* packedData);
  void unpack(packed_data* packedData);
  int nParticles() const {return 1;}
  spinOrbital& operator=( const spinOrbital &s);
  friend ostream& operator<<(ostream& output ,const spinOrbital &orbital2);
  friend istream& operator>>(istream& input ,spinOrbital &orbital2);
  
private:
  int spinValue;
  double x;
  double xBC;
};

template<class orbital_t>
class orbitals
{
public:
  
  orbitals(){setLabel("Empty Orbital");}
  orbitals(int n){orbitals();resize(n);};
  orbital_t& operator[](int i){return orbitalsStorage[i];};
  const orbital_t& operator[](int i) const {return orbitalsStorage[i];};
  void resize(int n){orbitalsStorage.resize(n);};

  int getPackSize()
  {
    orbital_t dummy;
    return( dummy.getPackSize()*size() );
  }
  
  void pack(packed_data* packO)
  {
    
    for(int i=0;i<size();i++)
      {
	orbitalsStorage[i].pack(packO);
      }
    
  }
  void unpack(packed_data* packO)
  {
    for(int i=0;i<size();i++)
      {
	orbitalsStorage[i].unpack(packO);
      }
  
  }
  
  void setLabel(const string &label_){label=label_;};
  string getLabel() const {return label;};
  
  size_t size() const {return orbitalsStorage.size();}
  
  template<class T>
  friend ostream& operator<<(ostream& output,const orbitals<T> & orbitals2);

  template<class T>
  friend istream& operator>>(istream& input,orbitals<T> & orbitals2);
  
  int nParticles() const
  {
    int n;
    n=0;
    for(int i=0;i<size();i++)
      {
	n+=orbitalsStorage[i].nParticles();
      }
    return n;
  }
private:
  vector<orbital_t> orbitalsStorage;
  string label;
  
};

#endif
