#ifndef MOVES_H
#define MOVES_H

#include <cassert>
#include <vector>
using namespace std;

template<class comp>
class qmcMover
{
public:
  
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename all_particles_t::particles_t particles_t;
  typedef typename comp::rand_t rand_t;
  typedef typename comp::grad_t grad_t;
  typedef typename particles_t::particles_t orbital_t;
  
  qmcMover(int seed,double delta_tau_)  : rand_o(seed)
  {
    setTimeStep(delta_tau);
  }
  
  double getTimeStep() const {return delta_tau;};
  void setTimeStep(double delta_tau_)
  {
    delta_tau=delta_tau_;
  }
  qmcMover(rand_t & rand_,double delta_tau_) : rand_o(rand_)
  {
    setTimeStep(delta_tau_);
  }
  
  virtual void reserve(size_t n)=0;
  virtual void move( all_particles_t & p, const grad_t & grad)=0;
  virtual void computeEffectiveDriftForce(all_particles_t &p,const grad_t & grad){cout << "No effective drift force to compute"<<endl;exit(2);};

  virtual void moveSpinRabi(particles_t & p,vector<double> &spinFlipsRatios){cout << "No spin rabi."<<endl;exit(2);};

  virtual void moveSpinRabi(orbital_t & p,double spinFlipsRatio){cout << "No individual spin rabi."<<endl;exit(2);};
  
  virtual grad_t & getEffectiveDriftForce(){cout << "No effective drift force"<<endl;exit(2);};
  
  virtual int getOrder()=0;
  
  virtual void moveGaussian(all_particles_t & p){};
protected:
  rand_t& getRandomEngine(){return rand_o;}
  double delta_tau;
  rand_t rand_o;
  
};

template<class comp>
class qmcMover1Order : public qmcMover<comp>
{
public:
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename comp::rand_t rand_t;
  typedef typename comp::grad_t grad_t;
  
  qmcMover1Order(int seed,double delta_tau_)  : qmcMover<comp>::qmcMover(seed,delta_tau_)
  {
    
  };
  
  qmcMover1Order(rand_t & rand_,double delta_tau_) : qmcMover<comp>::qmcMover(rand_,delta_tau_)
  {
    
  };
  virtual int getOrder()
  {
    return 1;
  }
  void reserve(size_t n)
  {
    work.resize(n);
  }
  
  // gaussian step
  void moveGaussian(all_particles_t & p)
  {
    int i,j,ni;
    for(i=0;i<p.size();i++)
      {
	
	ni=p[i].size();
	//cout << "ni: "<<ni << " " << work.size()<<endl;
       
	assert(work.size() >= ni );
	// generates the random numbers
	this->rand_o.gaussian(work,ni);
	// diffusion step
	for(j=0;j<ni;j++)
	  {
	    //cout << work[j]<<endl;
	    p[i][j].positionNoBC()+=sqrt(this->delta_tau)*work[j];
	  }
      }
  }
  
  // drift step
  
  void drift(all_particles_t & p,const grad_t & grad)
  {
    
    for(int i=0;i<p.size();i++)
      {
	for(int j=0;j<p[i].size();j++)
	  {
	    p[i][j].positionNoBC()+=this->delta_tau*grad[i][j];
	  }
      }
  }
  
  // first order move
  void move( all_particles_t & p, const grad_t & grad)
  {
    drift(p,grad);
    moveGaussian(p);
    
  }
  
private:
  
  vector<double> work;
  
};

#include "moveSecondOrder.h"

#include "movesSpin.h"

template<class comp>
qmcMover<comp>* buildQMCMover(comp* qmcObj)
{
  // load the type of mover from the input file
  string moverType;
  string qmcType;
  
  qmcType=qmcObj->main_input->reset()->get_child("method")->get_attribute("kind")->get_string();
  
  qmcObj->main_input->get_child("propagator");
  
  if ( qmcObj->main_input->check())
    {
      moverType=qmcObj->main_input->get_value()->get_string();
    }
  else
    {
      moverType="default";
    }
  
  cout << "moverType: "<<moverType<<endl;
  if (qmcType=="dmc" or qmcType=="svmc")
    {
      // build the mover
      if (moverType=="1order" or moverType=="default")
	{
	  return new qmcMover1Order<comp>(*qmcObj->rand,qmcObj->delta_tau);
	  
	}
      
      if (moverType=="2order")
	{
	  qmcMover2Order<comp>* ret;
	  ret=new qmcMover2Order<comp>(*qmcObj->rand,qmcObj->delta_tau,qmcObj->wave,qmcObj->geo);
	  ret->init(qmcObj->getInputFileName());
	  return ret;
	}
      
      if (moverType=="spin")
	{
	  qmcMover1OrderSpin<comp>* ret;
	  double omegaRabi=qmcObj->main_input->reset()->get_child("system")->get_child("oneBodyPotential")->get_child("oneBodyPotential")->get_attribute("omega")->get_real();
	  
	  double pFlip=sinh(omegaRabi*qmcObj->delta_tau)/( sinh(omegaRabi*qmcObj->delta_tau) + cosh(omegaRabi*qmcObj->delta_tau) );
	  
	  ret=new qmcMover1OrderSpin<comp>(*qmcObj->rand,qmcObj->delta_tau);
	  ret->setSpinTimeStep(pFlip);
	  return ret;
	}

      cout << "Unokmn propagator";
      exit(1);
      
    }
  else if (qmcType=="vmc")
    {
      if (moverType=="default")
	{
	  return new qmcMover1Order<comp>(*qmcObj->rand,qmcObj->delta_tau);
	}
      
      else if(moverType=="spin")
	{
	  return new qmcMover1OrderSpin<comp>(*qmcObj->rand,qmcObj->delta_tau);
	}
      else
	{
	  cout << "Unkown propagator"<<endl;
	  exit(1);
	}
    }
  
}

#endif
