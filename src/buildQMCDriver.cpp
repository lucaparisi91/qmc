
#include <iostream>
#include <string>
#include "xml-input.h"
#include "omp.h"
#include "traits.h"
#include "input.h"
#include "geometry.h"
using namespace std;

string createQMCDriverId()
{
  string bcInput;
  string particlesType;
  string calculation;
  bool isSpinor;
  string particleType;
  xml_input* main_input;
  main_input=new xml_input;
  main_input->open("input.xml");
  calculation=main_input->get_child("method")->get_attribute("kind")->get_string();
  if (main_input->get_attribute("spinor") )
    {
      isSpinor=main_input->get_bool();
    }
  else
    {
      isSpinor=false;
    }
  
  if (isSpinor)
    {
      particleType="spinor";
    }
  else
    {
      particleType=="normal";
    }
  
  bcInput=get_bc(main_input);
  
  return
    (
     string("bcInput=")+bcInput + string("particles=")+particlesType
     );
}
