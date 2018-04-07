#ifndef PIGSPROPAGATOR_H
#define PIGSPROPAGATOR_H

#include <cmath>
#include <cassert>
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
  
  void pairApproximationPropagator1D(double tau_)
  {
    tau=tau_;
  }
  
  inline double operator()(double x,double y)
  {
    double res=0;
    if (x*y>0)
      {
	res=1-exp(-(2*x*y)/tau);
      }
    assert(res>=0);
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
  void setTimeStep(double tau_){tau=tau_;}
  double getTimeStep() const {return tau;}
  
private:
 
  double tau;
  
};


class freePropagator1D
{
public:
  freePropagator1D(){A=0;tau=1;C=0;}
  void setTimeStep(double tau_){tau=tau_;A=1./sqrt(2*M_PI*tau);C=log(A);}
  
  inline double operator()(double z)
  {
    return exp(-pow(z,2))*A;
  }
  
  inline double logEvaluate(double z)
  {
    return -pow(z,2) + C;
  }
    
  void setSet(int set_){set=set_;}
  int getSet() const{return set;}
  
  private:
  
  double A;
  double tau;
  double C;
  int set;
  };


  
template<class comp>
class pairParticleApproximationChain
{
public:
  typedef typename comp::freeParticlePropagator_t freeParticlePropagator_t;
  typedef typename comp::pairParticlePropagator_t pairParticlePropagator_t;
  typedef typename comp::geometry_t geometry_t;
  typedef typename comp::wave_t wave_t;
  typedef typename comp::pos_t pos_t;
  typedef typename comp::configurations_t configurations_t;
  
  void evaluateLog(const configurations_t &configurations)
  {
    
    pos_t delta;
    // single particle propagator
    prob=0;
    for(int j=0;j< freeParticlePropagators.size();j++)
      {
	int i=freeParticlePropagators[j].getSet();
	// loop on imaginary time
	for (int l=0;l<configurations.size()-1;l++)
	  {
	    // loop on particles
	    for(int k=0;k<configurations[l][i].size();k++)
	      {
		delta=geo->distance_pbc(configurations[l][i][k].position(),configurations[l+1][i][k].position());
		prob+=freeParticlePropagators[j]->logEvaluate(delta);
	      }
	  }
	
      }
    
    pos_t delta2;
    // pair particle propagator
    for(int j=0;j<pairParticlePropagators.size();j++)
      {
	
	int i=pairParticlePropagators[j].getSet();
	// loop on imaginary time
	for (int l=0;l<configurations.size()-1;l++)
	  {
	    // loop on particles
	    for(int k1=0;k1<configurations[l][i].size();k1++)
	      {
		for(int k2=0;k2<k1;k2++)
		  {
		    delta=geo->distance_pbc(configurations[l][i][k1].position(),configurations[l][i][k2].position());
		    delta2=geo->distance_pbc(configurations[l+1][i][k1].position(),configurations[l+1][i][k2].position());
		    prob+=pairParticlePropagators[j]->logEvaluate(delta,delta2);
		  }
	      }
	
	  }
      }
    
    // tails
    prob+=wave->log_evaluate(configurations[0]);
    prob+=wave->log_evaluate(configurations[configurations.size()-1]);
    
  }
  
  void setGeometry(geometry_t & geo_){geo=geo_;};
  void setWave(wave_t & wave_){wave=wave_;};
  
  double getCurrentProbability() const {return prob;};
  
  void attachFreePropagator(freeParticlePropagator_t & gr){freeParticlePropagators.push_back(&gr);};
  void attachPairPropagator(pairParticlePropagator_t & gr){pairParticlePropagators.push_back(&gr);};
  //evaluate the log(pi(x_new)) - log(pi(x_current)) by changing one bead position
  
  //evaluate the ratio when moving just one particle
  double evaluateLogRatioHead(int set,int P,const pos_t & newPos,const configurations_t & configurations);
  
  double evaluateLogFreeLink(int iTime,int set, int iP,const configurations_t & configurations)
  {
    pos_t delta;
    delta=geo->distance_pbc(configurations[iTime][set][iP].position(),configurations[iTime+1][set][iP].position());
    return freeParticlePropagators[0]->logEvaluate(delta);
  }
  
  double evaluateLogPairLink(int iTime,int set, int iP,const configurations_t & configurations)
  {
    pos_t delta;
    pos_t delta1;
    double link=0;
    
    for(int i=0;i<configurations[iTime][set].size();i++)
      {
	
       
	delta=geo->distance_pbc(configurations[iTime][set][i].position(),configurations[iTime][set][iP].position());
	delta1=geo->distance_pbc(configurations[iTime+1][set][i].position(),configurations[iTime+1][set][iP].position());
	    
	    
	link+=pairParticlePropagators[0]->logEvaluate(delta,delta1);
      }
    
    return link;
  }
  
  double evaluateLogRatio(int iTime,int set,int iP,const pos_t & newPos,configurations_t & configurations)
  {
    pos_t delta,delta1,delta2,delta3;
    double logRatio=0;
    pos_t tmp;
    
    assert(iTime>0);
    assert(iTime<configurations.size()-1);
    
    //old links
    logRatio-=evaluateLogFreeLink(iTime,set, iP,configurations);
    logRatio-=evaluateLogFreeLink(iTime-1,set, iP,configurations);
    
    logRatio-=evaluateLogPairLink(iTime,set, iP,configurations);
    logRatio-=evaluateLogPairLink(iTime-1,set, iP,configurations);

    // new links
    tmp=configurations[iTime][set][iP].position();
    configurations[iTime][set][iP].position()=newPos;
    
    logRatio+=evaluateLogFreeLink(iTime,set, iP,configurations);
    logRatio+=evaluateLogFreeLink(iTime-1,set, iP,configurations);
    logRatio+=evaluateLogPairLink(iTime,set, iP,configurations);
    logRatio+=evaluateLogPairLink(iTime-1,set, iP,configurations);
    
    // resetting old configurations
    configurations[iTime][set][iP].position()=tmp;
    
    return logRatio;
    
  }
  
  void incrementRatio(double deltaProb){prob+=deltaProb;}
  
  private:
  geometry_t * geo;
  wave_t * wave;
  vector<freeParticlePropagator_t*> freeParticlePropagators;
  vector<pairParticlePropagator_t*> pairParticlePropagators;
  double currentRato;
  double prob;
  
};

  
template<class comp>
class pigsMover
  {
    typedef typename comp::rand_t rand_t;
    typedef typename comp::propagator_t propagator_t;
    typedef typename comp::configurations_t configurations_t;
    typedef typename comp::pos_t pos_t;
    
  public:

    pigsMover()
    {
      deltaRMax=0.1;
      nMoveOneParticle=0;
      nMoveOneParticleSuccess=0;
    }
    void move(configurations_t & configurations)
    {
      moveOneParticle(configurations);  
    }
    
    // moves particles one at a time
    void moveOneParticle(configurations_t & configurations)
    {
      bool accept;
      double logRatio=0;
      
      int l=int(randO->uniform()*(configurations.size()-2))+1;
      //int set=int(randO->uniform()*configurations[l].size());
      int set=0;
      int iP=int(randO->uniform()*configurations[l][set].size());
      pos_t Rnew=(randO->uniform()-0.5)*2*deltaRMax + configurations[l][set][iP].position();
      logRatio=propagator->evaluateLogRatio(l,set,iP,Rnew ,configurations);
      accept=metropolis(logRatio,randO);
      
      if (accept)
	{
	  configurations[l][set][iP].position()=Rnew;
	  propagator->incrementRatio(logRatio);
	  nMoveOneParticleSuccess+=1;
	}
      nMoveOneParticle+=1;
    }
    
  private:
    rand_t *randO;
    propagator_t *propagator;
    //---------- parameters for the list of events
    double deltaRMax; // amount to move a bead in PIGS
    double nMoveOneParticle,nMoveOneParticleSuccess;  
};

#endif
