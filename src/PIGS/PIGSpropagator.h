class complementaryErrorFunctionExponential
{
  void operator()()
  {
    return 2./(  sqrt(M_PI)*(u+sqrt(u*u+4/M_PI)));
  }
  
}
  
template<class configurations_t>
class pairApproximationPropagator1D
{
public:
  
  void pairApproximationPropagator1D(double g_,double tau_)
  {
    g=g_;
    tau=tau_;
    c=g*sqrt(pi*tau/2.);
    d=sqrt(1/(2*tau));
    e=g*tau;
  }
  
  void operator(int x,int y)
  {
    return 1-c*exp(-(x*y+abs(x*y))/tau)*cerrfexp(d*(abs(x)+abs(y)+e));
  }
  
  // returns the logarithm of the derivative of the propagator
  double logDerivative(double x,double y)
  {
    double u=d*(abs(x)+abs(y)+e);
    double h=(abs(x*y)+x*y)/(tau);
    return c*( h/tau*cerrfexp(u)+f(u)/(2*tau) - d*(abs(x)+abs(y)-e)/(2*tau)*(2*u*f(u)-2./M_PI))/(c*f(u)-exp(h));  
  }
      
private:
  
  double g;
  double tau;
  complementaryErrorFunctionExponential cerrfexp;
  
  double c;
  double d;
  double e;
  
};


