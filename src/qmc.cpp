#include "qmc.h"
#include "geometry.h"
#include "random.h"
#include "wavefunction.h"
#include "input.h"
#include "measures.h"
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include "xml-input.h"
#include "tools.h"
#include "mpi.h"
#include "timer.h"
#include "potential.h"

using namespace std;



empty_t* build_potential(xml_input* input,const empty_t*   pot)
{
  
  return new empty_t();
}

// build a potential from scratch
rabiCoupling* build_potential(xml_input* input,const rabiCoupling* pot)
{
  double omega;
  
  if(! input->reset()->get_child("system")->get_child("potential") )
    {
      cout << "No potential found";
      exit(1);
    }
  else
    {
      omega=input->get_attribute("omega")->get_real();
      return new rabiCoupling(omega);  
    }
}



// builds an harmonic potential from the ground

void print(rabiCoupling* om)
{
  cout << "Omega: "<< om->get_omega()<<endl;
};
void print(empty_t* p)
{
  
};


