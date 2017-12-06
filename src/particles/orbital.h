#ifndef ORBITAL_H
#define ORBITAL_H

#include <iostream>
#include <vector>
#include "../ptools.h"
#include "../tools.h"
#include "../xml-input.h"

using namespace std;
class spinOrbital
{
public:
  
  spinOrbital(){spinValue=0;x=0;xBC=0;}

  inline double& position(){return positionBC();};
  inline double& positionNoBC(){return x;};
  inline double& positionBC(){return x;};
  
  inline const double& position() const{return positionBC();};
  inline const double& positionNoBC() const{return x;};
  inline const double& positionBC() const{return x;};
  
  
  int& spin(){return spinValue;};
  const int& spin() const {return spinValue;};
  int get_pack_size() const;
  int size() const {return 1;}
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

using namespace std;
class orbital1D
{
  typedef double pos_t;
  
public:
  orbital1D(){x=0;xBC=0;}
  inline double& position(){return positionBC();};
  inline double& positionNoBC(){return x;};
  inline double& positionBC(){return xBC;};
  
  inline const double& position() const{return positionBC();};
  inline const double& positionNoBC() const{return x;};
  inline const double& positionBC() const{return xBC;};
  
  int size() const {return 1;}
  
  int get_pack_size() const;
  void pack(packed_data* packedData);
  void unpack(packed_data* packedData);
  int nParticles() const {return 1;}
  orbital1D & operator=( const orbital1D &s);
  friend ostream& operator<<(ostream& output ,const orbital1D &orbital2);
  friend istream& operator>>(istream& input ,orbital1D &orbital2);
  
private:
  
  double x;
  double xBC;
};



template<class orbital_t>
class orbitals
{
public: 
  typedef orbital_t particles_t;
  
  orbitals(){setLabel("EmptyOrbital");}
  orbitals(int n){orbitals();resize(n);};
  orbital_t& operator[](int i){return orbitalsStorage[i];};
  const orbital_t& operator[](int i) const {return orbitalsStorage[i];};
  void resize(int n){orbitalsStorage.resize(n);};

  int get_pack_size()
  {
    orbital_t dummy;
    return( dummy.get_pack_size()*size() );
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
  
  int sizes(vector<int> &ns) const
  {
    ns.resize(size());
    for(int i=0;i<size();i++)
      {
	ns[i]=orbitalsStorage[i].size();
      }
    
  }
  
  
  
  void push_back(const orbital_t &orbitalToAdd)
  {
    orbitalsStorage.push_back(orbitalToAdd); 
  }
  
private:
  vector<orbital_t> orbitalsStorage;
  string label;
  
};

void buildAllOrbitals(string filename,orbitals<orbitals<orbital1D> > & allOrbitals);




void generateUniform(orbitals<orbitals<orbital1D> > & orbitals,double xMin,double xMax);

void generateUniform(orbitals<orbital1D > & orbitals,double xMin,double xMax);



template<class randomGenerator_t>
void generateRandom(orbitals<orbital1D> & orbitals,double xMin,double xMax,randomGenerator_t* randO);

#endif
