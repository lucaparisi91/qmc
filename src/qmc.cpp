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



// builds an harmonic potential from the ground


