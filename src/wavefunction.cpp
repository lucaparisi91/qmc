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
#include "dmc.h"
#include "vmc.h"

/*
template<>
// kinetic energy using the fixed phase approximati
void total_wavefunction<vmc<spinor1D> >::kinetic_energy(total_wavefunction<vmc<spinor1D> >::all_particles_t* p, double &e,double &e_f)
{
  int i,j,k;
  double T,e_tmp;
  e=0;
  T=0;
  e_f=0;
  
  p->reset_drift_force();
  // computes some of the kinetic terms separately
  for(i=0;i<waves.size();i++)
    {
      waves[i]->kinetic_energy_first_term(p,T);
      e=e+T;
    }
  
  e_f=0;
  
  // compute the drift force term
  
  for(i=0;i<p->particle_sets.size();i++)
    {
      for(j=0;j<p->particle_sets[i]->n;j++)
	{
	  e_f+=real(pow(p->particle_sets[i]->drift_force[j],2));
	}
    }
  // sum the two contributions
  e=e + e_f;
  // multiply for the diffusion coefficient
  e=-e/2.;
  e_f=e_f/2;
}

template<>
// kinetic energy using the fixed phase approximati
void total_wavefunction<dmc<spinor1D> >::kinetic_energy(total_wavefunction<dmc<spinor1D> >::all_particles_t* p, double &e,double &e_f)
{
  int i,j,k;
  double T,e_tmp;
  e=0;
  T=0;
  e_f=0;
  
  p->reset_drift_force();
  // computes some of the kinetic terms separately
  for(i=0;i<waves.size();i++)
    {
      waves[i]->kinetic_energy_first_term(p,T);
      e=e+T;
    }
  
  e_f=0;
  
  // compute the drift force term
  
  for(i=0;i<p->particle_sets.size();i++)
    {
      for(j=0;j<p->particle_sets[i]->n;j++)
	{
	  e_f+=real(pow(p->particle_sets[i]->drift_force[j],2));
	}
    }
  // sum the two contributions
  e=e + e_f;
  // multiply for the diffusion coefficient
  e=-e/2.;
  e_f=e_f/2;
}

*/
/*
template<>
void load_wavefunctions(xml_input * xml_wave, vector<  dmc<spinor1D>::swave_t* > &waves,dmc<spinor1D>* qmc_obj)
{
  string label,kind;
  int setA;
  int setB;
  xmlNodePtr cur;
  
  xml_wave->reset()->get_child("wavefunctions")->get_child("wavefunction");
  do
    {
  label=xml_wave->get_attribute("label")->get_string();
  setA=xml_wave->get_attribute("setA")->get_int();
  setB=xml_wave->get_attribute("setB")->get_int();
  
  // adds a new wavefunction
  kind=xml_wave->get_attribute("kind")->get_string();
  
  if (kind=="bill_jastrow_spinor")
    {
      string fileParams;
      fileParams=xml_wave->get_attribute("file")->get_string();
      if (setA == setB)
		{
		  // load the required jastrows
		  cur=xml_wave->cur;
		  string fileUp,fileDown,kindDown,kindUp,kindUpDown,fileUpDown;
		  fileUp=xml_wave->get_child("jastrowUp")->get_attribute("file")->get_string();
		  kindUp=xml_wave->get_attribute("kind")->get_string();
		  xml_wave->cur=cur;
		  
		  fileDown=xml_wave->get_child("jastrowDown")->get_attribute("file")->get_string();
		  kindDown=xml_wave->get_attribute("kind")->get_string();
		  xml_wave->cur=cur;
		  
		  fileUpDown=xml_wave->get_child("jastrowUpDown")->get_attribute("file")->get_string();

		  kindUpDown=xml_wave->get_attribute("kind")->get_string();
		  xml_wave->cur=cur;
		  
		  if (kindUp == "delta")
		    {
		      waves.push_back(new bill_jastrow_wavefunction_spinor_two_body_symmetric< jastrow_spinor_free_sampling< jastrow_delta,jastrow_delta,jastrow_delta>,dmc<spinor1D> >(qmc_obj,fileUp,fileDown,fileUpDown,fileParams));
		    }
		  else
		    {
		      cout << "Not supported yet."<<endl;
		      exit(1);
		    }
		  
		}
	      else
		{
		  cout << "Not supported yet";
		  exit(1);
		}
    }
  
  waves.back()->src_particle_set=setA;
  waves.back()->target_particle_set=setB;
  waves.back()->label=label;
  xml_wave->get_next("wavefunction");
    }
  while( xml_wave->check() );
}

template<>
void load_wavefunctions(xml_input * xml_wave, vector<  vmc<spinor1D>::swave_t* > &waves,vmc<spinor1D>* qmc_obj)
{
  string label,kind;
  int setA;
  int setB;
  xmlNodePtr cur;
  
  xml_wave->reset()->get_child("wavefunctions")->get_child("wavefunction");
  do
    {
  label=xml_wave->get_attribute("label")->get_string();
  setA=xml_wave->get_attribute("setA")->get_int();
  setB=xml_wave->get_attribute("setB")->get_int();
  
  // adds a new wavefunction
  kind=xml_wave->get_attribute("kind")->get_string();
  
  if (kind=="bill_jastrow_spinor")
    {
      string fileParams;
      fileParams=xml_wave->get_attribute("file")->get_string();
      if (setA == setB)
		{
		  // load the required jastrows
		  cur=xml_wave->cur;
		  string fileUp,fileDown,kindDown,kindUp,kindUpDown,fileUpDown;
		  
		  fileUp=xml_wave->get_child("jastrowUp")->get_attribute("file")->get_string();
		  kindUp=xml_wave->get_attribute("kind")->get_string();
		  xml_wave->cur=cur;
		  
		  fileDown=xml_wave->get_child("jastrowDown")->get_attribute("file")->get_string();
		  kindDown=xml_wave->get_attribute("kind")->get_string();
		  xml_wave->cur=cur;
		  
		  fileUpDown=xml_wave->get_child("jastrowUpDown")->get_attribute("file")->get_string();

		  kindUpDown=xml_wave->get_attribute("kind")->get_string();
		  xml_wave->cur=cur;
		  
		  if (kindUp == "delta")
		    {
		      waves.push_back(new bill_jastrow_wavefunction_spinor_two_body_symmetric< jastrow_spinor_free_sampling< jastrow_delta,jastrow_delta,jastrow_delta>,vmc<spinor1D> >(qmc_obj,fileUp,fileDown,fileUpDown,fileParams));
		    }
		  else
		    {
		      cout << "Not supported yet."<<endl;
		      exit(1);
		    }
		  
		}
	      else
		{
		  cout << "Not supported yet";
		  exit(1);
		}
    }
  
  waves.back()->src_particle_set=setA;
  waves.back()->target_particle_set=setB;
  waves.back()->label=label;
  xml_wave->get_next("wavefunction");
    }
  while( xml_wave->check() );
}

*/

// template class total_wavefunction<dmc<D1_t> >;
// template class total_wavefunction< dmc<spinor1D> >;
// template class bill_jastrow_wavefunction<jastrow_delta,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_two_body_symmetric<jastrow_delta,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_two_body_symmetric<jastrow_delta_bound_state,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_two_body_symmetric<jastrow_delta_in_trap,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_two_body_symmetric<jastrow_delta_phonons,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_two_body_asymmetric<jastrow_delta,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_two_body_asymmetric<jastrow_delta_bound_state,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_two_body_asymmetric<jastrow_delta_phonons,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_two_body_asymmetric<jastrow_delta_in_trap,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_one_body<jastrow_delta,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_one_body<jastrow_delta_bound_state,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_one_body<jastrow_gaussian,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_one_body<jastrow_delta_in_trap,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_one_body<jastrow_delta_phonons,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_spinor_two_body_symmetric<jastrow_spinor_free_sampling<jastrow_delta,jastrow_delta,jastrow_delta> ,dmc<spinor1D> >;

string createWavefunctionId(xml_input* xml_wave)
{
  string kind;
  int setA,setB;
  string symmetricStatus;
  string nBody;
  int optParameter;
  string stringId;
  stringId="";
  optParameter=-1;
  kind=xml_wave->get_attribute("kind")->get_string();
  setA=xml_wave->get_attribute("setA")->get_int();
  setB=xml_wave->get_attribute("setB")->get_int();
  nBody=xml_wave->get_attribute("jastrow_type")->get_string();

  if (xml_wave->get_attribute("optParameter") != NULL )
    {
      
      optParameter=xml_wave->get_int();
    }
  else
    {
      optParameter=-1;
    }
  
  if (nBody=="1b")
    {
      stringId= kind + nBody;
    }
  else
    {
      if (setA==setB)
	{
	  symmetricStatus="symm";
	}
      else
	{
	  symmetricStatus="asymm";
	};
      
      // builds the unique name of the wavefunction
      stringId=kind + symmetricStatus + nBody;
      
    }
  
  if (optParameter >= 0)
    {
      stringId=stringId+string("Opt")+int_to_string(optParameter);
    }
  
  return stringId;
}

  // create the jastrow type id
string createJastrowId(xml_input* xml_wave)
{
  string jastrow_kind;
  jastrow_kind=xml_wave->get_attribute("kind")->get_string();
  
  
  return jastrow_kind;
  
};

template<>
void load_wavefunctions(xml_input * xml_wave, vector<  dmc<spinor1D>::swave_t* > &waves,dmc<spinor1D>* qmc_obj);

template<>
void load_wavefunctions(xml_input * xml_wave, vector<  vmc<spinor1D>::swave_t* > &waves, vmc<spinor1D>* qmc_obj);



