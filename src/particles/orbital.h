#ifndef ORBITAL_H
#define ORBITAL_H

#include <iostream>
#include <vector>
#include "../ptools.h"
#include "../tools.h"
#include "../xml-input.h"

using namespace std;
class orbital1D 
{
  typedef double pos_t;
  
public:
  orbital1D(){x=0;xBC=0;}
  inline double& position(){return positionBC();};
  inline double& positionNoBC(){return x;};
  inline double& positionBC(){return xBC;};
  
  inline const double& position() const {return positionBC();};
  inline const double& positionNoBC() const{return x;};
  inline const double& positionBC() const{return xBC;};
  
  int size() const {return 1;}
  
  virtual int get_pack_size() const;
  virtual void pack(packed_data* packedData);
  virtual void unpack(packed_data* packedData);
  int nParticles() const {return 1;}
  
  virtual ostream& output(ostream &out) const
  {
    out<< x <<" "<< xBC;
    return out;
  }
  
  virtual istream& input(istream &in)
  {
    in >> x;
    in >> xBC;
    return in;
  }
  
private:
  double x;
  double xBC;
};

ostream& operator<<(ostream& out    ,const orbital1D &orbital2);
istream& operator>>(istream& input ,orbital1D &orbital2);

using namespace std;
class spinOrbital1D : public orbital1D
{
public:
  typedef double pos_t;

  
  spinOrbital1D() : orbital1D(){spinValue=1;}
  
  inline int& spin(){return spinValue;}
  inline const int& spin() const {return spinValue;}
  
  virtual int get_pack_size() const;
  virtual void pack(packed_data* packedData);
  virtual void unpack(packed_data* packedData);
  
  virtual ostream& output(ostream& out) const
  {
    orbital1D::output(out);
    
    out<<" "<<spinValue;
    return out;
  }
  
  virtual istream& input(istream& in)
  {
    orbital1D::input(in);
    in>>spinValue;
    return in;
  }
  
private:
  int spinValue;
  
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
    int size=0;
    
    for(int i=0;i<orbitalsStorage.size();i++)
      {
	size+=orbitalsStorage[i].get_pack_size();
      }
    
    return size;
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
  
  void sizes(vector<int> &ns) const
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

  
int getMagnetization(orbitals<spinOrbital1D> &p);
int getMagnetization(orbitals<orbital1D> &p);



void setNDown(orbitals<spinOrbital1D> &orbitals,int n);

void setNDown(orbitals<orbitals<spinOrbital1D> > &orbitals,int n);


template<class orbital_t>
void buildAllOrbitals(string filename,orbitals<orbitals<orbital_t> > & allOrbitals)
{
  xml_input xml_main_input;
  orbitals<orbital_t > orbitals;
  int n;
  
  xml_main_input.open(filename)->reset();

  xml_main_input.reset()->get_child("system")->get_child("particles");
  
  do
    {
      xml_main_input.get_attribute("n");
      n=xml_main_input.get_int();
      orbitals.resize(n);
      allOrbitals.push_back(orbitals);
      
      xml_main_input.get_next("particles");
    }
  
  while(xml_main_input.check());
  
  xml_main_input.reset();
  
}

template<class orbital_t>
void generateUniform(orbitals<orbital_t > & orbitals,double xMin,double xMax)
{
  double step;
  step=(xMax-xMin)/orbitals.size();
  for(int i=0;i<orbitals.size();i++)
    {
      orbitals[i].positionNoBC()=xMin+i*step;
      
    }
}

template<class orbital_t>
void generateUniform(orbitals<orbitals<orbital_t> > & orbitals,double xMin,double xMax)
{
  
  for(int i=0;i<orbitals.size();i++)
    {
      generateUniform(orbitals[i],xMin,xMax);
    }
  
}

template<class randomGenerator_t,class orbital_t>
void generateRandom(orbitals<orbital_t> & orbitals,double xMin,double xMax,randomGenerator_t* randO)
{
  
  for(int i=0;i<orbitals.size();i++)
    {
      orbitals[i].positionNoBC()=randO->uniform()*(xMax-xMin) + xMin;
    }
}

template<class randomGenerator_t,class orbital_t>
void generateRandom(orbitals<orbitals<orbital_t> > & orbitals,double xMin,double xMax,randomGenerator_t* randO)
{
  
  for(int i=0;i<orbitals.size();i++)
    {
      generateRandom(orbitals[i],xMin,xMax,randO);
    }
}

template<>
void buildAllOrbitals(string filename,orbitals<orbitals<spinOrbital1D> > & allOrbitals);
#endif
