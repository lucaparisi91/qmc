#ifndef JASTROW_H
#define JASTROW_H

#include<vector>
#include<string>
#include <iostream>
#include <fstream>
#include "traits.h"
#include <complex>
#include "tools.h"
#include "xml-input.h"
#include "exception_classes.h"

template<class jastrowT>
class jastrow
{
  
  typedef typename traits<jastrowT>::position_t position_t;
  typedef typename traits<jastrowT>::value_t value_t;
  
 public:
  position_t a; // lower range
  position_t b; // higher range
  vector<real_t> parameters;
  real_t center;
  
  void print(string name_0d,string name_1d,string name_2d);//print to 3 files the jastrow value, first and second derivative
  inline value_t d0(const position_t &x){return static_cast<jastrowT*>(this)->d0(x);}; // value of the jastrow
  inline value_t d1(const position_t &x){return static_cast<jastrowT*>(this)->d1(x);}; // value of the jastrow
  inline value_t d2(const position_t &x){return static_cast<jastrowT*>(this)->d2(x);}; // value of the jastrow
  
  void setParameter(double x,int i)
  {
    throw notYetSsupported("setParameter");
  }
  
  double getParameter(int i) const
  {
    throw notYetSsupported("getParameter");
  }
  
};

class jastrow_gaussian : public jastrow<jastrow_gaussian>
{
public:
  typedef double position_t;
  typedef double value_t;
  jastrow_gaussian(string filename);
  void load_parameters(string filename);
  inline double d0(const double &x){return exp(-x*x*parameters[0]);};
  double d1d0(const double &x);
  double d2d0(const double &x);
  inline double d1(const double &x){return -2*x*parameters[0]*exp(-x*x*parameters[0]);}
  inline double d2(const double &x){return -2*parameters[0]*(1-2*parameters[0]*x*x)*exp(-x*x*parameters[0]);};
  double dP1(const double &x){return 0;};
  double dP2(const double &x){return 0;};
  void setParameter(double x,int i)
  {

    if(i==0)
      {
	parameters[0]=x;
      }
    else
      {
	if(i==1)
	  {
         	center=x;
	      
	  }
	else
	  {
	    if(i==2)
	      {
		center=-x;
	      }
	    else
	      {
		throw unkown_parameter(int_to_string(i));
	      }
	  }
      }
  
   
  }
  
 double getParameter(int i) const
  {
    if(i==0)
      {
	return parameters[0];
      }
    else
      {
	if(i==1)
	  {
	    return center;
	  }
	
	else
	  {
	    if (i==2)
	      {
		return -center;
	      }
	    else
	      {
		throw unkown_parameter(int_to_string(i));
	      }
	  }
      }
  }

};

class jastrowSpinOrbital : public jastrow<jastrowSpinOrbital>
{
public:
  typedef double position_t;
  typedef double value_t;
  
  jastrowSpinOrbital(string filename)
  {
    parameters.resize(3);
    load_parameters(filename);
  }
  
  void load_parameters(string filename)
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
  
  inline double d0(const double &x){return parameters[0];};
  double d1d0(const double &x){return 0.;};
  double d2d0(const double &x){return 0.;};
  inline double d1(const double &x){return 0.;};
  inline double d2(const double &x){return 0;};
  double dP1(const double &x){return 0;};
  double dP2(const double &x){return 0;};
  
  void setParameter(double x,int i)
  {
    assert(i==0);
    parameters[i]=x;
  }
  
  double getParameter(int i) const
  {
    assert(i==0);
    return parameters[0];
  }
  
};


class jastrow_delta : public jastrow<jastrow_delta>
{
  typedef double position_t;
  typedef real_t value_t;
  
 public:
  
  jastrow_delta(string filename);
  void load_parameters(string filename);
  inline double d0(const double &x){return (x < parameters[3]) ? sin(parameters[0]*x + parameters[1]): c;};
  
  double d1d0(const double &x);

  inline double d1(const double &x){return (x<parameters[3]) ? cos(parameters[0]*x + parameters[1])*parameters[0]: 0;};

  inline double d2(const double &x){return (x < parameters[3]) ? -sin(parameters[0]*x + parameters[1])*parameters[0]*parameters[0]: 0;};
  
  
  double d2d0(const double &x);
  double dP1(const double &x){return 0;};
  double dP2(const double &x){return 0;};
  void setParameter(double x,int i){}
  double getParameter(int i) const {return 0;};
private:
  double d;
  double c;
};

class jastrow_delta_in_trap : public jastrow<jastrow_delta_in_trap>
{
  typedef double position_t;
  typedef real_t value_t;
public:
  
  jastrow_delta_in_trap(string filename);
  void load_parameters(string filename);
 
  inline double d0(const double &x){return (x + this->parameters[0]);}
  inline double d1(const double &x){return 1;};
  inline double d1d0(const double &x){return 1/(x + this->parameters[0]);};
  inline double d2d0(const double &x){return 0;};
  inline double d2(const double &x){return 0;};
  inline double dP1(const double &x){return 0;};
  inline double dP2(const double &x){return 0;};
  void setParameter(double x,int i){}
  double getParameter(int i) const {return 0;};
  
};

class jastrow_delta_phonons : public jastrow<jastrow_delta_phonons>
{
  typedef double position_t;
  typedef real_t value_t;
public:
  jastrow_delta_phonons(string filename);
  void load_parameters(string filename);
  double dP1(const double &x){return 0;};
  double dP2(const double &x){return 0;};
  
  void setParameter(double x,int i)
  {
    
    if (i==0)
      {
	cout << "set to "<< x << endl;
	parameters[3]=x ;
	process();
      }
    else
      {
	cout << "Unkown parameter to set. Exit !";
      }
    
  }
  
  double getParameter(int i) const
  {
    if (i==0)
      {
	return parameters[3];
      }
    else
      {
	cout << "Unkown parameter. Exit !";
      }
  };
  
  __attribute__((always_inline)) double d0(const double &x)
  {
    return (x < parameters[3]) ? sin(parameters[0]*x + parameters[1]) :  pow(sin(k2*x),parameters[2]);
  }
  
  __attribute__((always_inline)) double d1(const double &x)
  {
    return (x < parameters[3]) ? parameters[0]*cos(parameters[0]*x + parameters[1]) : pow(sin(k2*x),parameters[2]-1)*parameters[2]*cos(k2*x)*k2;
  }
  
  __attribute__((always_inline)) double d2(const double &x)
  {
    return (x < parameters[3]) ? -parameters[0]*parameters[0]*sin(parameters[0]*x + parameters[1]) : pow(sin(k2*x),parameters[2]-2)*parameters[2]*(parameters[2]-1)*pow(cos(k2*x)*k2,2) - pow(sin(k2*x),parameters[2])*parameters[2]*k2*k2;
  }
  
  void process()
  {
    double aP=1e-5;
    double bP=3*aP;
    
    fR.z=parameters[3];
    fR.g=parameters[5];
    fR.l_box=parameters[4];
    
    while( fR(bP)*fR(aP) >=0 )
    {
      bP*=2;
    }
    
    parameters[0]=findRootBrente(fR,aP,bP,1e-4);
    parameters[1]=atan(parameters[0]/parameters[5]);

    parameters[2]=(parameters[0]*tan((M_PI/parameters[4])*parameters[3]))/(M_PI/parameters[4]  * tan(parameters[0]*parameters[3]+parameters[1]));
    
    k2=M_PI/parameters[4];
    center=0;
    a=0;
    b=parameters[4]/2.;
  }

  
  class fRoot
  {
  public:
    double operator()(double x)
    {
      double delta,beta;
      delta=atan(x/g);
      beta=(x*tan((M_PI/l_box)*z))/( (M_PI/l_box)  * tan(x*z+delta));
      return (sin(x*z+delta) - pow(sin((M_PI/l_box)*z),beta));
    }
    
    double z;
    double l_box;
    double g;  
  };
  
private:
  
  double k2;
  fRoot fR;
};

class jastrow_delta_bound_state : public jastrow<jastrow_delta_bound_state>
{
public:
  typedef double position_t;
  typedef real_t value_t;
  jastrow_delta_bound_state(string filename);
  void load_parameters(string filename);
  void process(double);
  inline double d0(const double x){double y=x-parameters[3]; return (y<0) ? (exp(-parameters[0]*y) + exp(parameters[0]*y))/2 : 1;}
  
  inline double d1(const double x)  { double y=x-parameters[3]; return (y<0) ? parameters[0]*(exp(parameters[0]*y) - exp(-parameters[0]*y)   )/2 : 0;
  }
  
  inline double d2(const double x)
  { double y=x-parameters[3]; return (y<0) ? parameters[0]*parameters[0]*(exp(parameters[0]*y) + exp(-parameters[0]*y))/2 : 0;}
  
  // set the the parameter
  void setParameter(double x,int i)
  {
    parameters[3]=x;
    process(x);
  }
  
  double getParameter(int i) const {return parameters[3];};
  double dP1(const double &x);
  double dP2(const double &x);
  
  double d2d0(const double &x);
  double d1d0(const double &x);
private:
  
  class fRoot
  {
  public:
    double operator()(double x)
    {
      return x-g/(2*tanh(x*cut_off));
    }
    
    double cut_off;
    double g;
  };
  fRoot f;
protected:
  double dkdXI; // dependance on k from the free parameter cut_off
  double dkdXI2;
};

template<class jastrowSame,class jastrowDifferent>
class jastrowRealSpinorSym
{
public:
  void setParameter(double x,int i){}
  double getParameter(int i) const {return 0;};
  typedef tinyVector<double,2> spinComp_t;
  
  jastrowRealSpinorSym(string filenameSame,string filenameDiff) : jS(filenameSame),jD(filenameDiff){};
  
    inline double d0(const double &x,const spinComp_t & s1,const spinComp_t & s2) const
    {
      return jS.d0(x)*(s1[0]*s2[0] + s1[1]*s2[1] ) + jD.d0(x)*(s1[0]*s2[1] + s1[1]*s2[0]);
    }

    inline double d1(const double &x,const spinComp_t & s1 , const spinComp_t & s2) const
    {
      return jS.d1(x)*(s1[0]*s2[0] + s1[1]*s2[1] ) + jD.d1(x)*(s1[0]*s2[1] + s1[1]*s2[0]);
    }
  
    inline double d2(const double &x,const spinComp_t & s1 , const spinComp_t & s2) const
    {
      return jS.d2(x)*(s1[0]*s2[0] + s1[1]*s2[1] ) + jD.d2(x)*(s1[0]*s2[1] + s1[1]*s2[0]);
    }
  
private:
  jastrowDifferent jD;
  jastrowSame jS;
  // to be used as variational parameters for the system
  double a;
  double b;
  
};

template<class jastrowUp,class jastrowDown,class jastrowUpDown>
class jastrow_spinor
{
 public:
  void setParameter(double x,int i){}
  double getParameter(int i) const {return 0;};
  
  void load_parameters(string filename);
  // creates a spinor of jastrow wavefunctions
  jastrow_spinor(string filenameUp,string filenameDown,string filenameUpDown,string fileParams ) : jastrow1(filenameUp),jastrow2(filenameDown),jastrow3(filenameUpDown){load_parameters(fileParams);a=1;b=-1;}
  inline complex<double> d0(const double &x,const vector< complex<double> > &s1,const vector<complex<double> > &s2){return a*a*jastrow1.d0(x)*s1[0]*s2[0] +b*b*jastrow2.d0(x)*s1[1]*s2[1] + a*b*jastrow3.d0(x)*(s1[0]*s2[1] + s1[1]*s2[0]);
  }
  inline complex<double> d1(const double &x,const vector< complex<double> > &s1,const vector<complex<double> > &s2){return a*a*jastrow1.d1(x)*s1[0]*s2[0] +b*b*jastrow2.d1(x)*s1[1]*s2[1] + a*b*jastrow3.d1(x)*(s1[0]*s2[1] + s1[1]*s2[0])  ;}
  
  inline complex<double> d2(const double &x,const vector< complex<double> > &s1,const vector< complex<double> > &s2){return a*a*jastrow1.d2(x)*s1[0]*s2[0] +b*b*jastrow2.d2(x)*s1[1]*s2[1] + a*b*jastrow3.d2(x)*(s1[0]*s2[1] + s1[1]*s2[0])  ;}
  
  void print(string name_0d,string name_1d,string name_2d)
  {
    jastrow1.print(name_0d + "Up.dat",name_1d + "Up.dat",name_2d + "Up.dat");
    jastrow2.print(name_0d + "Down.dat",name_1d + "Down.dat",name_2d + "Down.dat");
    jastrow3.print(name_0d + "UpDown.dat",name_1d + "UpDown.dat",name_2d + "UpDown.dat");
    
  }


  inline double dP1(const double &x){return 0;};
  inline double dP2(const double &x){return 0;};
  
protected:
  jastrowUp jastrow1;
  jastrowDown jastrow2;
  jastrowUpDown jastrow3;
  
  complex<double> a;
  complex<double> b;
  
};



// spinor jastrow with the possibility to compute the spin gradient for a certain wavefunction
template<class jastrowUp,class jastrowDown,class jastrowUpDown>
class jastrow_spinor_free_sampling : public jastrow_spinor<jastrowUp,jastrowDown,jastrowUpDown>
{
public:
  void setParameter(double x,int i){}
  double getParameter(int i) const {return 0;};
  inline double dP1(const double &x){return 0;};
  inline double dP2(const double &x){return 0;};
  jastrow_spinor_free_sampling(string filenameUp,string filenameDown,string filenameUpDown,string fileParams ) : jastrow_spinor<jastrowUp,jastrowDown,jastrowUpDown>(filenameUp,filenameDown,filenameUpDown,fileParams){};
  // first spin derivative respect to the first spin component
  inline complex<double> d11(const double &x,const vector< complex<double> > &s1,const vector<complex<double> > &s2){complex<double> i(0,1);return i*this->a*this->a*this->jastrow1.d0(x)*s1[0]*s2[0] -i* this->b*this->b*this->jastrow2.d0(x)*s1[1]*s2[1] + i*this->a*this->b*this->jastrow3.d0(x)*(s1[0]*s2[1] - s1[1]*s2[0]);}
  // first spin derivative with respect to the second spin component
  inline complex<double> d12(const double &x,const vector< complex<double> > &s1,const vector<complex<double> > &s2){complex<double> i(0,1);return i*this->a*this->a*this->jastrow1.d0(x)*s1[0]*s2[0] -i* this->b*this->b*this->jastrow2.d0(x)*s1[1]*s2[1] + i*this->a*this->b*this->jastrow3.d0(x)*(-s1[0]*s2[1] + s1[1]*s2[0]);}
  
};


template<class jastrowT>
void jastrow<jastrowT>::print(string name_0d,string name_1d,string name_2d)
{
  // print the jastrows on the files
  int bins=1000;
  int i;
  double step,x;
  ofstream file_0d,file_1d,file_2d;
  
  file_0d.open(name_0d.c_str());
  file_1d.open(name_1d.c_str());
  file_2d.open(name_2d.c_str());
  
  step=(b-a)*1./bins;
  
  for (i=0;i<=bins;i++)
    {
      // x position on the array
      x=i*step+a;
      // print out to a file
      file_0d << x <<" "<< d0(x)<<endl;
      file_1d << x <<" "<< d1(x)<<endl;
      file_2d << x << " "<< d2(x)<<endl;
      
    }
  
}


template<class jastrowUp,class jastrowDown,class jastrowUpDown>
void  jastrow_spinor<jastrowUp,jastrowDown,jastrowUpDown>::load_parameters(string filename)
{
  xml_input *xml_jastrow=new xml_input;
  xml_jastrow->open(filename);
  // get the 'a' and 'b' coefficients for the system
  a=xml_jastrow->reset()->get_child("a")->get_value()->get_complex();
  b=xml_jastrow->reset()->get_child("b")->get_value()->get_complex();
  
}

class emptyJastrow
{
public:
  typedef double value_t;
  emptyJastrow(const vector<double> &parameters){}
  double d0(double x){throw noOptimization();return 0;};
  double d1(double x){throw noOptimization();return 0;};
  double d2(double x){throw noOptimization();return 0;};
  double pD0(double x){throw noOptimization();return 0;};
  void setParameter(double x,int i){throw noOptimization();};
  
};

// optimization of the jastro parameter
template<class jastrow_t,class jastrowDerivative=emptyJastrow>
class jastrowOptimized : public jastrow_t
{
  
public:
jastrowOptimized(string filename) : jastrow_t(filename)
  {
    jD=new jastrowDerivative(this->parameters);
  };
  
typedef typename jastrowDerivative::value_t value_t;
  value_t dP0(double x){return jD->d0(x);};
  value_t dP1(double x){return jD->d1(x);};
  value_t dP2(double x){return jD->d2(x);};
  value_t pD0(double x){return jD->pD0(x);};
  
  void setParameter(double x,int i){
    jastrow_t::setParameter(x,i);
      }
  ;
  
  double getParameter(int i) const {return jastrow_t::getParameter(i);};
  
public :
jastrowDerivative *jD;
};

// optimized jastrow with a bound state
template<>
class jastrowOptimized<jastrow_delta_bound_state,emptyJastrow> : public jastrow_delta_bound_state
{
public:
  jastrowOptimized(string filename) : jastrow_delta_bound_state(filename){};
  
  value_t dP0(double x)
  {
    double k;
    double y;
    y=x-parameters[3];
    k=parameters[0];
    return tanh(k*y)*(dkdXI*y-k);
  };
  
  value_t dP1(double x)
  {
    double k;
    double y;
    y=x-parameters[3];
    k=parameters[0];
    return (dkdXI*tanh(k*y)+(dkdXI*y-k)*k/pow(cosh(k*y),2));
  }
  
  value_t dP2(double x)
  {
    double k;
    double y;
    y=x-parameters[3];
    k=parameters[0]; 
    return k*(2*dkdXI-k*tanh(k*y)*(dkdXI*y-k) )/pow(cosh(k*y),2);
  };
  
  value_t pD0(double x)
  {
    double k;
    double y;
    y=x-parameters[3];
    k=parameters[0]; 
    return (dkdXI2*y-2*dkdXI)*tanh(k*y)+pow((dkdXI-k),2)/(pow(cosh(k*y),2));
  }
  
};


// derivate with respect to zero
class jastrowDerivativeAlphaGaussian
{
public:
  jastrowDerivativeAlphaGaussian(const vector<double> & parameters_)
  {
    parameters=parameters_;
  }
  
  typedef double value_t;
  value_t d0(double x){return -x*x;};
  value_t d1(double x){return -2*x;};
  value_t d2(double x){return -2;};
  // returns the value of the derivative in a certain point
  value_t pD0(double x){return 0;}
  void setParameter(double x,int i)
  {
    if (i==0)
      {
	parameters[0]=x;
      }
    else
      {
	
      }
  }
  
private:
  // parameters to optimize on
  vector<double> parameters;
};
// derive weith respect to the center of the gaussian
class jastrowDerivativeBetaGaussian
{
public:
  jastrowDerivativeBetaGaussian(const vector<double> & parameters_)
  {
    parameters=parameters_;
  }
  typedef double value_t;
  
  value_t d0(double x){return 2*parameters[0]*x;}
  value_t d1(double x){return 2*parameters[0];};
  value_t d2(double x){return 0;};
  value_t pD0(double x){return -2*parameters[0];}
  void setParameter(double x,int i)
  {
    if (i==0)
      {
	parameters[0]=x;
      }
    else
      {
	
      }
  }
  
private:
  vector<double> parameters;
};

// derivative of -beta instead of beta
class jastrowDerivativeBetaGaussianReverse
{
public:
  jastrowDerivativeBetaGaussianReverse(const vector<double> & parameters_)
  {
    parameters=parameters_;
  }
  typedef double value_t;
  
  value_t d0(double x){return -2*parameters[0]*x;}
  value_t d1(double x){return -2*parameters[0];};
  value_t d2(double x){return 0;};
  value_t pD0(double x){return 2*parameters[0];}
  void setParameter(double x,int i)
  {
    if (i==0)
      {
	parameters[0]=x;
      }
    else
      {
	
      }
  }
  
private:
  vector<double> parameters;
};

#include "jastrow_spline.h"
#include "jastrowSpinOrbital.h"
#endif

