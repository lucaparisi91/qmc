#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include "xml-input.h"

using namespace std;
template<class T>
class pbcTools
{
public:
  template<class obj> void bc(vector<obj> &vec)
  {
    int i;
    
    for (i=0;i<vec.size();i++)
    {
      vec[i]=static_cast<T*>(this)->bc(vec[i]);
    }
  }
  
  template<class particles_t>
  void particles_bc(particles_t& p)
  {
    int i;
    for(i=0;i<p.size();i++)
      {
	p[i].positionBC()=static_cast<T*>(this)->bc(p[i].positionNoBC());
      }
    
  }
  
  // pbc conditions on all particles
  template<class all_particles_t>
  void all_particles_bc(all_particles_t& p)
  {
    int i=0;
    for(i=0;i<p.size();i++)
      {
	particles_bc(p[i]);
      }
  }
  
};

class pbc1d : public pbcTools<pbc1d> 
{
private:
  double l_box;
public:
  
  pbc1d(double l_box_){l_box=l_box_;};
  pbc1d& operator= ( pbc1d &pbc2){l_box=pbc2.l_box;};
  double bc(double x);
  
  inline double distance_bc(double x1,double x2)
  {
    double x;
    x=x1-x2;
    //cout << x1<<" "<<x2<< " "<<x<<" "<<l_box<<endl;
    
    return x>l_box/2. ? x-l_box : ((x<-l_box/2.) ? x+l_box : x) ;
  }
  
};

string get_bc(xml_input* xml_input);
// uses no periodic boundary conditions

class noPbcD1 : public pbcTools<noPbcD1>
{
public:
  
  inline double bc(double x){return x;};
  inline double distance_bc(double x1,double x2){return (x1-x2);};
  noPbcD1(double l_box){};
  
};

template<class bc_class>
class geometry
{
 public:
  typedef bc_class bc_t;
  
  int dimensions; // dimensionality of the system
  double l_box;//length of the box
  geometry(double l_box_){dimensions=1;l_box=l_box_;bco=new bc_t(l_box_);};
  geometry<bc_class> & operator= ( geometry<bc_class> &geo2);
    
  geometry(string filename);
  
  inline double pbc(double x1)
  {
    return bco->bc(x1);
  }
  // returns the value applied period boundary conditions
  //inline double pbc2(double x)
  //{
  //  if (x<-l_box/2) x=x+l_box;
  //  if (x>= l_box/2) x=x-l_box;
    
  //  return x;
  //}
  
  inline void pbc(vector<double> &vec)
  {
    bco->bc(vec);
  }
  inline double distance_pbc(double x1,double x2)
  {
    return bco->distance_bc(x1,x2);
  }

  template<class particles_t>
  void particles_pbc(particles_t & p)
  {
    bco->particles_bc(p);
  }

  template<class all_particles_t>
  void all_particles_pbc(all_particles_t& p)
  {
    bco->all_particles_bc(p);
  }

  double getBoxDimensions()
  {
    return l_box;
  }
  void setBoxDimension(double lBox2)
  {
    l_box=lBox2;
  }
  
private:
  bc_t *bco;
};
template<class bc_class>
geometry<bc_class>::geometry(string filename)
{
  int n_particles;
  xml_input* main_input;
  main_input=new xml_input;
  main_input->open(filename);
  n_particles=0;
  n_particles=main_input->reset()->get_child("system")->get_child("particles")->get_attribute("n")->get_int();
  if (
      main_input->reset()->get_child("system")->get_child("lBox") != NULL
       )
    {
      l_box=main_input->reset()->get_child("system")->get_child("lBox")->get_value()->get_real();
    }
  else
    {
      l_box=n_particles;
    }
  //cout<<l_box;
  dimensions=1;
  bco=new bc_t(l_box);
  delete main_input;
}


#endif
