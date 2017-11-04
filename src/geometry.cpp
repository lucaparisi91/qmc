#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include "input.h"
#include "geometry.h"
#include <sstream>
#include "xml-input.h"
#include "qmc.h"

double pbc1d::bc(double x)
{
  // returns the periodic boundary conditions
  int n_box;
  double y;
  
  n_box=floor((abs(x) + l_box*0.5)/l_box);
  y=abs(x)+l_box*0.5-l_box*n_box;
  
  if ( x >= 0)
    {
      x=-l_box*0.5 + y;
    }
  else
    {
      x=l_box*0.5-y;
    }
  return x;
  assert(abs(x) <= l_box*1./2);
}
template<class bc_class>
geometry<bc_class> & geometry<bc_class>::operator= ( geometry<bc_class> & geo2)
{
  l_box=geo2.l_box;
  (*bco)=(*geo2->bco);
}

string get_bc(xml_input* main_input)
{
  xmlNodePtr cur,cur2;
  string bc;
  cur=main_input->cur;
  main_input->reset()->get_child("system")->get_child("bc");
  if ( main_input->cur != NULL )
    {
      bc=main_input->get_value()->get_string();
    }
  else
    {
      bc="periodic";
    }
  
  main_input->cur=cur;
  return bc;
}
