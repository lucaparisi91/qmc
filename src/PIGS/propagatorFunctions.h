#ifndef PROPAGATORFUNCTIONS_H
#define PROPAGATORFUNCTIONS_H

#include <cmath>
#include <cassert>
#include <vector>

using namespace std;

class complementaryErrorFunctionExponential
{
public:
  inline double operator()(double u)
  {
    return 2./(  sqrt(M_PI)*(u+sqrt(u*u+4/M_PI)));
  }
  
};

class pairApproximationPropagator1D
{
public:
  
  pairApproximationPropagator1D(double g_,double tau_)
  {
    g=g_;
    tau=tau_;
    c=g*sqrt(M_PI*tau/2.);
    d=sqrt(1/(2*tau));
    e=g*tau;
  }
  
  inline double operator()(double x,double y)
  {
    double res=1-c*exp(-(x*y+abs(x*y))/tau)*cerrfexp(d*(abs(x)+abs(y)+e));
    assert(res>0);
    return res;
  }
  
  inline double logEvaluate(double x,double y)
  {
    return log(operator()(x,y));
  }
  
  /*
  // returns the logarithm of the derivative of the propagator
  double logDerivative(double x,double y)
  {
    double u=d*(abs(x)+abs(y)+e);
    double h=(abs(x*y)+x*y)/(tau);
    return c*( h/tau*cerrfexp(u)+f(u)/(2*tau) - d*(abs(x)+abs(y)-e)/(2*tau)*(2*u*f(u)-2./M_PI))/(c*f(u)-exp(h));  
  }
  */  
private:
  
  double g;
  double tau;
  complementaryErrorFunctionExponential cerrfexp;
  
  double c;
  double d;
  double e;
  
};

class pairApproximationPropagator1DUnitary
{
public:
  typedef system_t::particles_t::particles_t orbital_t;
  typedef system_t::geometry_t geometry_t;
  
  pairApproximationPropagator1DUnitary(){};
  
  inline double operator()(const orbital_t & x1,const orbital_t & x2,const orbital_t & x3,const orbital_t & x4)
  {
    double x,y;
    double xnpbc,ynpbc;
    double res=0;
    x=abs(geo->distance_pbc(x1.position(),x2.position()));
    y=abs(geo->distance_pbc(x3.position(),x4.position()));
    
    xnpbc=x1.positionNoBC()-x2.positionNoBC();
    ynpbc=x3.positionNoBC()-x4.positionNoBC();
    
    if (xnpbc*ynpbc<0) return 0;
    //if (x*y>30) {return 1;};
    
    res=1-exp(-(2*x*y)/tau);
    
    
    assert(res>=0);
    return res;
  }
  
  inline double logEvaluate(const orbital_t & x1,const orbital_t & x2,const orbital_t & x3,const orbital_t & x4)
  {
    return log(operator()(x1,x2,x3,x4));
  }

  inline double logEvaluateTimeDerivative(const orbital_t & x1,const orbital_t & x2,const orbital_t & x3,const orbital_t & x4)
  {
    double x,y;
    double res=0;
    x=abs(geo->distance_pbc(x1.position(),x2.position()));
    y=abs(geo->distance_pbc(x3.position(),x4.position()));
    
    //x=x1.positionNoBC()-x2.positionNoBC();
    //y=x3.positionNoBC()-x4.positionNoBC();
    
    res=exp(-2*x*y/tau);
    return -2.*x*y*res/((1-res)*tau*tau);
    
  }
  /*
  // returns the logarithm of the derivative of the propagator
  double logDerivative(double x,double y)
  {
    double u=d*(abs(x)+abs(y)+e);
    double h=(abs(x*y)+x*y)/(tau);
    return c*( h/tau*cerrfexp(u)+f(u)/(2*tau) - d*(abs(x)+abs(y)-e)/(2*tau)*(2*u*f(u)-2./M_PI))/(c*f(u)-exp(h));  
  }
  */
  
  void setTimeStep(double tau_){tau=tau_;}
  double getTimeStep() const {return tau;}

  void setGeometry(geometry_t* geo_){geo=geo_;}
  geometry_t* getGeometry(){return geo;}
private:
  geometry_t* geo;
  double tau;
};



class freePropagator1D
{
public:
  typedef system_t::geometry_t geometry_t;
  typedef system_t::particles_t::particles_t orbital_t;
  
  freePropagator1D(){A=0;tau=1;C=0;}
  void setTimeStep(double tau_){tau=tau_;A=1./sqrt(2*M_PI*tau);C=log(A);D=0.5/tau;}
  
  inline double operator()(orbital_t& orbital1,orbital_t& orbital2)
  {
    double z=abs(geo->distance_pbc(orbital1.position(),orbital2.position()));
    return exp(-pow(z,2)*D)*A;
  }
  
  inline double logEvaluate(orbital_t& orbital1,orbital_t& orbital2)
  {
    double z=abs(geo->distance_pbc(orbital1.position(),orbital2.position()));
    return -z*z*D + C;
  }
  
  inline double logEvaluateTimeDerivative(double dis)
  {
    return -D+dis*dis/(2*tau*tau);
  }
  
  void setGeometry(geometry_t* geo_){geo=geo_;}
  geometry_t* getGeometry(){return geo;}
  
  private:
  geometry_t* geo;
  double A;
  double D;
  double tau;
  double C;
  
  };

#endif
