#include "functions.h"
#include "xml-input.h"
#include <cmath>

delta_wavefunction::delta_wavefunction(string filename)
  {
    parameters.resize(4);
    load_parameters(filename);
    
  }

void delta_wavefunction::load_parameters(string filename)
  {
    xml_input *xml_jastrow=new xml_input;
    xml_jastrow->open(filename);
    parameters[0]=xml_jastrow->reset()->get_child("k")->get_value()->get_real();
    parameters[1]=xml_jastrow->reset()->get_child("delta")->get_value()->get_real();
    parameters[2]=xml_jastrow->reset()->get_child("l_box")->get_value()->get_real();
    parameters[3]=xml_jastrow->reset()->get_child("cut_off")->get_value()->get_real();
    
    a=0;
    b=parameters[2]/2;
  }
  
double delta_wavefunction::d0(double x)
  {
    
    if (x < parameters[3])
      {
      return sin(parameters[0]*x + parameters[1]);
      }
    else
      {
      return sin(parameters[0]*parameters[3] + parameters[1]);
      }
  }
  
double delta_wavefunction::d1(double x)
  {
    if (x < parameters[3])
      {
      return parameters[0]*cos(parameters[0]*x + parameters[1]);
      }
    else
      {
	return 0;
	}
  }
  
double delta_wavefunction::d2(double x)
  {
    if (x < parameters[3])
      {
	return -pow(parameters[0],2)*sin(parameters[0]*x + parameters[1]);
      }
    else
      {
      return 0;
      }
  }

// defines a function from the product of two functions. Take care of the derivatives
template <class T>
double double_func<T>::d0(double x)
{
  return f1->d0(x)*f2->d0(x);
}
template <class T>
double double_func<T>::d1(double x)
{
  return f1->d1(x)*f2->d0(x) + f1->d0(x)*f2->d1(x);
}
template <class T>
double double_func<T>::d2(double x)
{
  return f1->d2(x)*f2->d0(x) + f1->d1(x)*(f2->d2(x) + f2->d1(x)) + f1->d1(x)*f2->d2(x);
}
