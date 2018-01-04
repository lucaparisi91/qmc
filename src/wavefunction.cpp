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
  string kind1="";
  jastrow_kind=xml_wave->get_attribute("kind")->get_string();
  
  if ( xml_wave->get_attribute("kind1") != NULL )
    {
      kind1=xml_wave->get_string();
    }
  jastrow_kind+=kind1;
  
  return jastrow_kind;
  
};

template<>
void load_wavefunctions(xml_input * xml_wave, vector<  dmc<spinor1D>::swave_t* > &waves,dmc<spinor1D>* qmc_obj);

template<>
void load_wavefunctions(xml_input * xml_wave, vector<  vmc<spinor1D>::swave_t* > &waves, vmc<spinor1D>* qmc_obj);



