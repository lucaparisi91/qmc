#include "qmc.h"
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <iostream>
#include "random.h"
#include "input.h"
#include "jastrow.h"
#include "wavefunction.h"
#include "geometry.h"
#include "mesh.h"
#include "measures.h"
#include "xml-input.h"
#include "system.h"
#include "tools.h"
#include "factory.h"
#include "exception_classes.h"

using namespace std;
/*
wavefunction
- compute the value of the wavefunctions
- defines the measurements on the wavefunction(may not depened on the particular kind of wavefunction)
*/

template<class jastrow_t,class comp>
bill_jastrow_wavefunction<jastrow_t,comp>::bill_jastrow_wavefunction(bill_jastrow_wavefunction<jastrow_t,comp>::qmc_t* qmc_obj_, xml_input* xml_wave,string filename) : wavefunction<comp>(qmc_obj_), jastrowc(filename)
{
  string kind;
  string file;
  xmlNodePtr cur;
  
  cur=xml_wave->cur;
  this->src_particle_set=xml_wave->get_attribute("setA")->get_int();
  this->target_particle_set=xml_wave->get_attribute("setB")->get_int();
  
  xml_wave->get_child("jastrow");
  kind=xml_wave->get_attribute("kind")->get_string();
  file=xml_wave->get_attribute("file")->get_string();
  xml_wave->cur=cur;
  
}

template<class jastrow_t,class comp>
double bill_jastrow_wavefunction_two_body_asymmetric<jastrow_t,comp>::one_particle_log_evaluate(bill_jastrow_wavefunction_two_body_asymmetric<jastrow_t,comp>::all_particles_t * p,double r,int iP, int set_a)
{
  int set_to_cycle;
  double ro,x,value;
  
  int j;
  
  value=0;
  
  if (set_a==this->src_particle_set)
    {
      set_to_cycle=this->target_particle_set;
    }
  else
    {
      if (set_a==this->target_particle_set)
	{
	  set_to_cycle=this->src_particle_set;
	  
	}
      else
	{
	  return value;
	}
    }

  particles_t& p1=(*p)[set_to_cycle];
  
  for (j=0;j<p1.size();j++)
    {
      x=abs(this->qmc_obj->geo->distance_pbc(r,p1[j].position()));

      value+=log(this->jastrowc.d0(x));
      
    }
  return value;
}

template<class jastrow_t,class comp>
// returns the logarithm of the wavefunction
double bill_jastrow_wavefunction_two_body_asymmetric<jastrow_t,comp>::log_evaluate(bill_jastrow_wavefunction_two_body_asymmetric<jastrow_t,comp>::all_particles_t * p)
{
  int i,j;
  double x,tmp,value;
  particles_t& p1=(*p)[this->src_particle_set];
  particles_t& p2=(*p)[this->target_particle_set];
  
  value=0;    
  
  // assert that the required jastrows are defined
  for (i=0;i<p1.size();i++)
    {
      
      for (j=0;j<p2.size();j++)
	{
	  x=abs(this->qmc_obj->geo->distance_pbc(p1[i].position(),p2[j].position()));
	  value=value+log(this->jastrowc.d0(x));
	}	  
    }
  
  this->log_wavefunction_value=value;
  return value;
}

template<class jastrow_t,class comp>
double bill_jastrow_wavefunction_two_body_symmetric<jastrow_t,comp>::one_particle_log_evaluate(bill_jastrow_wavefunction_two_body_symmetric<jastrow_t,comp>::all_particles_t * p,double r,int iP, int set_a)
{
  double ro,x,value;
  
  particles_t& p1=(*p)[this->src_particle_set];
  int j;
  
  value=0;
  
  if (set_a != this->src_particle_set)
    {
      return 0;
    }
  
  
  for (j=0;j<iP;j++)
    {
      
      x=abs(this->qmc_obj->geo->distance_pbc(r,p1[j].position()));
      

      value+=log(this->jastrowc.d0(x));
     
    }
  
  for (j=iP+1;j<p1.size();j++)
    {
      
      x=abs(this->qmc_obj->geo->distance_pbc(r,p1[j].position()));

      value+=log(this->jastrowc.d0(x));
     
    }
  return value;
}


template<class X,class comp>
double bill_jastrow_wavefunction_two_body_symmetric<X,comp>::log_evaluate(bill_jastrow_wavefunction_two_body_symmetric<X,comp>::all_particles_t * p)
{
  int i,j;
  double x,tmp,value;
  
  particles_t& p1=(*p)[this->src_particle_set];
  value=0;
  
  // assert that the required jastrows are defined
  for (i=0;i<p1.size();i++)
    {
      
      for (j=0;j<i;j++)
      {
	x=abs(this->qmc_obj->geo->distance_pbc(p1[i].position(),p1[j].position()));
	
	value=value+log(this->jastrowc.d0(x));
      }
    
    }
  this->log_wavefunction_value=value;
  return value;
  
}

template<class X,class comp>
double bill_jastrow_wavefunction_spinor_two_body_symmetric<X,comp>::log_evaluate(bill_jastrow_wavefunction_spinor_two_body_symmetric<X,comp>::all_particles_t * p)
{
  int i,j;
  double x,tmp,value;
  particles_t * p1;
  
  value=0;
  
  p1=p->particle_sets[this->src_particle_set];  
  
  // assert that the required jastrows are defined
  for (i=0;i<p1->n;i++)
    {
      for (j=0;j<i;j++)
      {
	x=abs(this->qmc_obj->geo->distance_pbc(p1->position[i],p1->position[j]));
	value=value+log(abs(this->jastrowc.d0(x,p1->spinComp[i],p1->spinComp[j])));
      }
    
    }
  
  this->log_wavefunction_value=value;
  return value;
  
}


template<class X,class comp>
double bill_jastrow_wavefunction_one_body<X,comp>::one_particle_log_evaluate(bill_jastrow_wavefunction_one_body<X,comp>::all_particles_t * p,double r,int iP, int set_a)
{
  double x;
  
  
  if (set_a != this->src_particle_set)
    {
      return 0;
    }
  
 
 
  
  x=abs(this->qmc_obj->geo->distance_pbc(r,this->jastrowc.center));
  return (log(this->jastrowc.d0(x))); 
}

template<class X,class comp>
double bill_jastrow_wavefunction_one_body<X,comp>::log_evaluate(bill_jastrow_wavefunction_one_body<X,comp>::all_particles_t * p)
{
  int i,j;
  double x,tmp,value;
  particles_t& p1=(*p)[this->src_particle_set];
 
  value=0;
  // assert that the required jastrows are defined
  for (i=0;i<p1.size();i++)
    {  
      //x=(qmc_obj->geo->pbc(p1->position[i]));
      x=abs(this->qmc_obj->geo->distance_pbc(p1[i].position(),this->jastrowc.center));
      value=value+log(this->jastrowc.d0(x));
    }
  
  this->log_wavefunction_value=value;
  return value;
  
}
template<class X,class comp>
double bill_jastrow_wavefunction_one_body<X,comp>::evaluate_derivative(bill_jastrow_wavefunction_one_body<X,comp>::all_particles_t * p)
{
  int i,j;
  double x,tmp,value;
  particles_t& p1=(*p)[this->src_particle_set];
  
  
  value=0;
  // assert that the required jastrows are defined
  for (i=0;i<p1.size();i++)
    {  
      //x=(qmc_obj->geo->pbc(p1->position[i]));
      x=this->qmc_obj->geo->distance_pbc(p1[i].position(),this->jastrowc.center);
      
      value=value+this->jastrowc.dP0(x);
    }
  
  this->wavefunction_derivative=value;
  //cout << value<<endl;
  return value;
  
}

template<class X,class comp>
double bill_jastrow_wavefunction_two_body_asymmetric<X,comp>::evaluate_derivative(bill_jastrow_wavefunction_two_body_asymmetric<X,comp>::all_particles_t * p)
{
  int i,j;
  double x,tmp,value;
  
  particles_t& p1=(*p)[this->src_particle_set];
  particles_t& p2=(*p)[this->target_particle_set];
  
  value=0;
  // assert that the required jastrows are defined
  for (i=0;i<p1.size();i++)
    {
      for(j=0;j<p2.size();j++)
	{
	  //x=(qmc_obj->geo->pbc(p1->position[i]));
	  x=this->qmc_obj->geo->distance_pbc(p1[i].position(),p2[j].position());
	  value=value+this->jastrowc.dP0(x);
	}
    }

  
  
  this->wavefunction_derivative=value;
  //cout << value<<endl;
  return value;
  
}


// evaluates the second derivative
template<class X,class comp>
double bill_jastrow_wavefunction_one_body<X,comp>::evaluate_derivative_second(bill_jastrow_wavefunction_one_body<X,comp>::all_particles_t * p,const double phi)
{
  int i,j;
  double x,tmp,value;
  particles_t& p1=(*p)[this->src_particle_set];
  
  
  
  value=0;
  // assert that the required jastrows are defined
  for (i=0;i<p1.size();i++)
    {  
      //x=(qmc_obj->geo->pbc(p1->position[i]));
      x=this->qmc_obj->geo->distance_pbc(p1[1].position(),this->jastrowc.center);
      
      value=value+this->jastrowc.pD0(x);
    }
  
  //this->wavefunction_derivative=value;
  //cout << value<<endl;
  return (value +phi*phi) ;
  
}

template<class X,class comp>
double bill_jastrow_wavefunction_two_body_asymmetric<X,comp>::evaluate_derivative_second(bill_jastrow_wavefunction_two_body_asymmetric<X,comp>::all_particles_t * p,const double phi)
{
  int i,j;
  double x,tmp,value;
  
  
  particles_t& p1=(*p)[this->src_particle_set];
  particles_t& p2=(*p)[this->target_particle_set];
    
  value=0;
  
  // assert that the required jastrows are well defined
  for (i=0;i<p1.size();i++)
    {
      for(j=0;j<p2.size();j++)
	{
	  //x=(qmc_obj->geo->pbc(p1->position[i]));
	  x=this->qmc_obj->geo->distance_pbc(p1[i].position(),p2[j].position());
	  value=value+this->jastrowc.pD0(x);
	}
    }
  
  //this->wavefunction_derivative=value;
  //cout << value<<endl;
  return (value +phi*phi) ;
  
}


// evaluates the logarithm of the wavefunction
template <class comp>
double total_wavefunction<comp>::log_evaluate(total_wavefunction<comp>::all_particles_t* p)
{
  int i;
  double w;
  w=0;
  
  for(i=0;i<waves.size();i++)
    {
      w+=waves[i]->log_evaluate(p);
    }
  return w;
}
// returns the sign
// ! only after the wave function has already been evaluated

template<class comp>
int total_wavefunction<comp>::get_sign()
{
  int i;
  int sign=1;
  for(i=0;i<waves.size();i++)
    {
      sign*=waves[i]->sign;
    }
  return sign;
}
// returns the phase
// ! only after the wave function has already been evaluated
template<class comp>
double total_wavefunction<comp>::get_phase()
{
  int i;
  double phase=0;
  for(i=0;i<waves.size();i++)
    {
      phase+=waves[i]->phase;
    }
  return phase;
  
}
 
// loads from the input file the wavefunction with label 'label'
template<class comp>
total_wavefunction<comp>::total_wavefunction(typename total_wavefunction<comp>::qmc_t * qmc_obj_ )
{
  string kind;
  string jastrow_type;
  int setA;
  int setB;
  qmc_obj=qmc_obj_;
  
}

template<class comp>
void total_wavefunction<comp>::link_wavefunctions(xml_input* xml_wave,vector<wavefunction<comp>*> wave_table,string label)
{
  int i;
  bool found;
  xmlNodePtr cur;
  bool optimization;
  
  cur=xml_wave->cur;
  xml_wave->reset()->get_child("total_wavefunction","label",label);
  xml_wave->get_child("wavefunction");
  do
    {
      found=false;
      
      label=xml_wave->get_attribute("label")->get_string();
      
      for(i=0;i<wave_table.size();i++)
	{
	  if (label==wave_table[i]->label)
	    {
	      waves.push_back(wave_table[i]);
	      found=true;
	      break;
	    }
	}
      // if the jastrow is not found at all
      
      if (!found)
	{
	  throw jastrowNotFound(label);
	}
      xml_wave->get_next("wavefunction");
    }
  while(xml_wave->check());
  xml_wave->cur=cur;
  
}

template<class comp>
wavefunction<comp>::wavefunction(wavefunction<comp>::qmc_t * qmc_obj_)
{
  qmc_obj=qmc_obj_;
  sign=1;
  phase=0;
  toOptimize=false;
  energy_first_term=0;
  wavefunction_derivative=0;
  log_wavefunction_value=0;
}

template<class comp>
// prints all availible kind of jastrows
void total_wavefunction<comp>::print_jastrows()
{
  int i;
  
  for (i=0;i<waves.size();i++)
    {
      waves[i]->print(i);
    }
}

// clean the state of the system

// loads wavefunctions from an external file

template<class jastrow_t,class comp>
void bill_jastrow_wavefunction<jastrow_t,comp>::print(int i)
{
  
  jastrowc.print(string("jastrow")+ int_to_string(i)+ string(".dat"),string("jastrow_1d")+ int_to_string(i) + string(".dat"),string("jastrow_2d")+ int_to_string(i) + ".dat");
}

template<class comp>
void load_wavefunctions(xml_input * xml_wave, vector< typename comp::swave_t* > &waves,comp* qmc_obj)
{
  
  int setA,setB;
  string jastrow_type;
  string kind,label,jastrow_kind;
  xmlNodePtr cur;
  string filename;
  string wavefunctionId;
  string jastrowId;
  int optParameter;
  waveFactory< comp > fac;
  
  // register the correct type to be saved
  fac.registerType("bill_jastrowsymm2bdelta",& ( createBillJastrowTwoBodySymmetric<jastrowOptimized<jastrow_delta>,comp> ) );
  fac.registerType("bill_jastrowsymm2bspline",& ( createBillJastrowTwoBodySymmetric<jastrowOptimized<jastrow_spline>,comp> ) );
  fac.registerType("bill_jastrowasymm2bspline",& ( createBillJastrowTwoBodyAsymmetric<jastrowOptimized<jastrow_spline>,comp> ) );
  fac.registerType("bill_jastrowasymm2bdelta",& ( createBillJastrowTwoBodyAsymmetric<jastrowOptimized<jastrow_delta>,comp> ) );
  fac.registerType("bill_jastrowasymm2bdelta_phonons",& ( createBillJastrowTwoBodyAsymmetric<jastrowOptimized<jastrow_delta_phonons>,comp> ) );
  fac.registerType("bill_jastrowsymm2bdelta_phonons",& ( createBillJastrowTwoBodySymmetric<jastrowOptimized<jastrow_delta_phonons>,comp> ) );
  fac.registerType("bill_jastrowsymm2bdelta_bound_state",& ( createBillJastrowTwoBodySymmetric<jastrowOptimized<jastrow_delta_bound_state>,comp> ) );
  fac.registerType("bill_jastrowasymm2bdelta_bound_state",& ( createBillJastrowTwoBodyAsymmetric< jastrowOptimized<jastrow_delta_bound_state>,comp> ) );
  
  fac.registerType("bill_jastrowasymm2bdelta_bound_state_no_pbc",& ( createBillJastrowTwoBodyAsymmetric<jastrow_delta_bound_state_no_pbc,comp> ) );

  fac.registerType("bill_jastrowasymm2bdelta_bound_state_no_pbc2",& ( createBillJastrowTwoBodyAsymmetric<jastrow_delta_bound_state_no_pbc2,comp> ) );

  fac.registerType("bill_jastrowasymm2bdelta_bound_state_no_pbc3",& ( createBillJastrowTwoBodyAsymmetric<jastrow_delta_bound_state_no_pbc3,comp> ) );
  
  fac.registerType("bill_jastrowsymm2bdelta_bound_state_no_pbc",& ( createBillJastrowTwoBodySymmetric<jastrow_delta_bound_state_no_pbc,comp> ) );
  
  fac.registerType("bill_jastrowasymm2bdelta_in_trap_exponential",& ( createBillJastrowTwoBodyAsymmetric< jastrow_delta_in_trap_exponential,comp> ) );
  
  fac.registerType("bill_jastrowsymm2bdelta_in_trap",& ( createBillJastrowTwoBodySymmetric<jastrowOptimized<jastrow_delta_in_trap> ,comp> ) );

  fac.registerType("bill_jastrowsymm2bdelta_in_trap_exponential",& ( createBillJastrowTwoBodySymmetric<jastrow_delta_in_trap_exponential,comp> ) );

  fac.registerType("bill_jastrow_spin_orbitalsymm2bdelta_in_trap",& ( buildBillJastrowSpinTwoBody<jastrowSpinOrbitalTwoBody<jastrow_delta_in_trap,jastrow_delta_in_trap> ,comp> ) );

  fac.registerType("bill_jastrow_spin_orbitalsymm2bdelta_in_trapdelta_in_trap",& ( buildBillJastrowSpinTwoBody<jastrowSpinOrbitalTwoBody<jastrow_delta_in_trap,jastrow_delta_in_trap> ,comp> ) );

  fac.registerType("bill_jastrow_spin_orbitalsymm2bdelta_phonons",& ( buildBillJastrowSpinTwoBody<jastrowSpinOrbitalTwoBody<jastrow_delta_phonons,jastrow_delta_phonons> ,comp> ) );
  

  fac.registerType("bill_jastrow_spin_orbitalsymm2bdelta",& ( buildBillJastrowSpinTwoBody<jastrowSpinOrbitalTwoBody<jastrow_delta,jastrow_delta> ,comp> ) );
  
  fac.registerType("bill_jastrow_spin_orbitalsymm2bdeltadelta_bound_state",& ( buildBillJastrowSpinTwoBody<jastrowSpinOrbitalTwoBody<jastrow_delta,jastrow_delta_bound_state> ,comp> ) );
  
  fac.registerType("bill_jastrow_spin_orbital1bgauss",& ( buildBillJastrowSpinOneBody<jastrowSpinOrbitalOneBody<jastrow_gaussian,jastrow_gaussian> ,comp> ) );

  fac.registerType("bill_jastrow_spin_orbital1bgaussgauss",& ( buildBillJastrowSpinOneBody<jastrowSpinOrbitalOneBody<jastrow_gaussian,jastrow_gaussian> ,comp> ) );

    fac.registerType("bill_jastrow_spin_orbital1bspinOrbital",& ( buildBillJastrowSpinOneBody<jastrowSpinOrbitalOneBody<jastrowSpinOrbital,jastrowSpinOrbital> ,comp> ) );
  
  fac.registerType("bill_jastrowasymm2bdelta_in_trap",& ( createBillJastrowTwoBodyAsymmetric<jastrowOptimized<jastrow_delta_in_trap>,comp> ) );
  fac.registerType("bill_jastrow1bdelta_in_trap",& ( createBillJastrowOneBody<jastrowOptimized<jastrow_delta_in_trap>,comp> ) );
  fac.registerType("bill_jastrow1bdelta",& ( createBillJastrowOneBody<jastrowOptimized<jastrow_delta>,comp> ) );
  fac.registerType("bill_jastrow1bdelta_bound_state",& ( createBillJastrowOneBody<jastrowOptimized<jastrow_delta_bound_state>,comp> ) );
  fac.registerType("bill_jastrow1bgauss",& ( createBillJastrowOneBody<jastrowOptimized<jastrow_gaussian>,comp> ) );
  fac.registerType("bill_jastrow1bOpt0gauss",& ( createBillJastrowOneBody<jastrowOptimized<jastrow_gaussian,jastrowDerivativeAlphaGaussian>,comp> ) );
  fac.registerType("bill_jastrow1bOpt1gauss",& ( createBillJastrowOneBody<jastrowOptimized<jastrow_gaussian,jastrowDerivativeBetaGaussian>,comp> ) );
  fac.registerType("bill_jastrow1bOpt2gauss",& ( createBillJastrowOneBody<jastrowOptimized<jastrow_gaussian,jastrowDerivativeBetaGaussianReverse>,comp> ) );
  
  xml_wave->reset()->get_child("wavefunctions")->get_child("wavefunction");
  
  do
    {
      
      label=xml_wave->get_attribute("label")->get_string();
      
      cout << label<<endl;
      // adds a new wavefunction
      kind=xml_wave->get_attribute("kind")->get_string();
      jastrow_type=xml_wave->get_attribute("jastrow_type")->get_string();
      setA=xml_wave->get_attribute("setA")->get_int();
      setB=xml_wave->get_attribute("setB")->get_int();
      
      if (xml_wave->get_attribute("optParameter")!=NULL)
	{
	  optParameter=xml_wave->get_int();
	}
      else
	{
	  optParameter=-1;
	}
      
      wavefunctionId=createWavefunctionId(xml_wave);
      
      cur=xml_wave->cur;
      xml_wave->get_child("jastrow");
      jastrowId=createJastrowId(xml_wave);
      
      filename=xml_wave->get_attribute("file")->get_string();
      jastrow_kind=xml_wave->get_attribute("kind")->get_string();
      
      xml_wave->cur=cur;
      
      /*waves.push_back(
	new bill_jastrow_wavefunction_two_body_symmetric<
	jastrow_delta,comp>(qmc_obj,xml_wave,filename));
      */
      
      waves.push_back(fac.create(wavefunctionId + jastrowId,qmc_obj,xml_wave,filename));
      
      waves.back()->src_particle_set=setA;
      waves.back()->target_particle_set=setB;
      waves.back()->label=label;
      
      // flag the wavefunction if an optimization parameter is set
      if (optParameter>=0)
	{
	  (waves.back())->setOptimized(true);
	  waves.back()->setOptParameter(optParameter);
	}
      
      xml_wave->get_next("wavefunction");
      
    }
  
  while( xml_wave->check() );
  
}
