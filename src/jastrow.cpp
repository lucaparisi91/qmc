#include <iostream>
#include <vector>
#include <vector>
#include <cassert>
#include <fstream>
#include "jastrow.h"
#include "xml-input.h"
#include <cmath>
#include "typelist.h"

using namespace std;
//-------------------------------delta jastrow definition--------

jastrow_delta::jastrow_delta(string filename)
  {
    parameters.resize(5);
    load_parameters(filename);
  }
jastrow_delta_in_trap::jastrow_delta_in_trap(string filename)
  {
    parameters.resize(3);
    load_parameters(filename);
  }
jastrow_delta_phonons::jastrow_delta_phonons(string filename)
  {
    parameters.resize(6);
    load_parameters(filename);
    //process();
    
    cout  << "k: "<< parameters[0]<<endl;
    cout  << "beta: "<< parameters[2]<<endl;
    cout  << "delta: "<< parameters[1]<<endl;
    
  }

// load parameter
void jastrow_delta_in_trap::load_parameters(string filename)
  {
    xml_input *xml_jastrow=new xml_input;
    xml_jastrow->open(filename);
    parameters[1]=xml_jastrow->reset()->get_child("l_box")->get_value()->get_real();
    parameters[0]=xml_jastrow->reset()->get_child("c")->get_value()->get_real();
    parameters[2]=xml_jastrow->reset()->get_child("position")->get_value()->get_real();
    
    center=parameters[2];
    a=0;
    b=parameters[1]/2.;
  }

void jastrow_delta_phonons::load_parameters(string filename)
{
  xml_input *xml_jastrow=new xml_input;
  xml_jastrow->open(filename);
  parameters[0]=xml_jastrow->reset()->get_child("k")->get_value()->get_real();
  parameters[1]=xml_jastrow->reset()->get_child("delta")->get_value()->get_real();
  parameters[2]=xml_jastrow->reset()->get_child("beta")->get_value()->get_real();
  parameters[3]=xml_jastrow->reset()->get_child("z")->get_value()->get_real();
  parameters[4]=xml_jastrow->reset()->get_child("l_box")->get_value()->get_real();
  parameters[5]=xml_jastrow->reset()->get_child("g")->get_value()->get_real();
  
  k2=M_PI/parameters[4];
  center=0;
  a=0;
  b=parameters[4]/2.;
  
}

void jastrow_delta::load_parameters(string filename)
  {
    xml_input *xml_jastrow=new xml_input;
    xml_jastrow->open(filename);
    parameters[0]=xml_jastrow->reset()->get_child("k")->get_value()->get_real();
    parameters[1]=xml_jastrow->reset()->get_child("delta")->get_value()->get_real();
    parameters[2]=xml_jastrow->reset()->get_child("l_box")->get_value()->get_real();
    parameters[3]=xml_jastrow->reset()->get_child("cut_off")->get_value()->get_real();
    parameters[4]=xml_jastrow->reset()->get_child("position")->get_value()->get_real();
    center=parameters[4];
    a=0;
    b=parameters[2]/2;
    c=sin(parameters[0]*parameters[3] + parameters[1]);
    d=-pow(parameters[0],2);
  }


double jastrow_delta::d1d0(const double &x)
  {
    
    if (x < parameters[3])
      {
      return parameters[0]/tan(parameters[0]*x + parameters[1]);
      }
    else
      {
	return 0;
	}
  }

double jastrow_delta::d2d0(const double &x)
  {
    
    if (x < parameters[3])
      {
	return d;
      }
    else
      {
	return 0;
	}
  }


//----- delta bound state jastrow
jastrow_delta_bound_state::jastrow_delta_bound_state(string filename)
  {
    
    parameters.resize(5);
    load_parameters(filename);    
  }

void jastrow_delta_bound_state::load_parameters(string filename)
{
  xml_input *xml_jastrow=new xml_input;
  xml_jastrow->open(filename);
  parameters[0]=0; // k
  parameters[1]=0; // B
  parameters[4]=xml_jastrow->reset()->get_child("g")->get_value()->get_real();
  parameters[2]=xml_jastrow->reset()->get_child("l_box")->get_value()->get_real();
  parameters[3]=xml_jastrow->reset()->get_child("cut_off")->get_value()->get_real();
  
  a=0;
  b=parameters[2]/2;
  
  process(parameters[3]);
  
}

void jastrow_delta_bound_state::process(double xI)
{
  double aP,bP,k,g,tmp;
  f.cut_off=parameters[3];
  f.g=parameters[4];
  aP=0.01;
  bP=10;
  
  while(f(aP)>0)
    {
      aP=aP/2;
    }
  
  while( f(bP)<0 )
    {
      bP=bP*2;
    }
  
  parameters[0]=findRootBrente(f,aP,bP,0.001);
  parameters[1]=exp(-2*parameters[0]*parameters[3]);
  
  // set up paraemters for the derivative with respect to the parameter
  k=parameters[0];
  g=parameters[4];
  xI=parameters[3];
  // first derivative
  dkdXI=(g*pow(cosh(k*xI),2)-k*k)/(sinh(2*k*xI)/2 + k*xI);
  // second derivative
  dkdXI2= (g*sinh(2*k*xI)*(dkdXI*xI + k)-2*k*dkdXI-dkdXI*(k+xI*dkdXI)*(1+cosh(2*k*xI)))/(sinh(2*k*xI)*0.5 + xI*k);
  
  
}



double jastrow_delta_bound_state::d1d0(const double &x)
  {
    double y;
    y=x-parameters[3];
    if (y<0)
      {
	return parameters[0]*tanh(parameters[0]*y);
      }
    else
      {
	return 0; 
      }
  }

double jastrow_delta_bound_state::d2d0(const double &x)
  {
    double y;
    y=x-parameters[3];
    if (y<0)
      {
	return parameters[0]*parameters[0];
      }
    else
      {
	return 0; 
      }
  }

// gaussian jastrow wavefunction

jastrow_gaussian::jastrow_gaussian(string filename)
{
  
  parameters.resize(2);
  load_parameters(filename);
}

void jastrow_gaussian::load_parameters(string filename)
{
  double lBox;
  xml_input *xml_jastrow=new xml_input;
  xml_jastrow->open(filename);
  
  parameters[0]=xml_jastrow->reset()->get_child("alpha")->get_value()->get_real();
  
  this->center=xml_jastrow->reset()->get_child("position")->get_value()->get_real();
  delete xml_jastrow;
  xml_jastrow=new xml_input;
  xml_jastrow->open("input.xml");
  lBox=xml_jastrow->reset()->get_child("system")->get_child("lBox")->get_value()->get_real();
  a=0;
  b=lBox/2;
  delete xml_jastrow;

  
}

double jastrow_gaussian::d0(const double &x)
{
  
  return exp(-x*x*parameters[0]);
}
double jastrow_gaussian::d1(const double &x)
{
  return -2*x*parameters[0]*exp(-x*x*parameters[0]);
}
double jastrow_gaussian::d1d0(const double &x)
{
  return -2*x*parameters[0];
}

double jastrow_gaussian::d2d0(const double &x)
{
  return -2*parameters[0]*(1-2*parameters[0]*x*x);
}

double jastrow_gaussian::d2(const double &x)
{
  return -2*parameters[0]*(1-2*parameters[0]*x*x)*exp(-x*x*parameters[0]);
}


