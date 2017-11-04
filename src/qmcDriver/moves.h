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
  typedef typename comp::rand_t rand_t;
  typedef typename comp::grad_t grad_t;
  
  qmcMover(int seed,double delta_tau_)  : rand_o(seed)
  {
    setTimeStep(delta_tau);
  }
  
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
  
  virtual grad_t & getEffectiveDriftForce(){cout << "No effective drift force"<<endl;exit(2);};
  
  virtual int getOrder()=0;
protected:
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
    for(i=0;i<p.getN();i++)
      {
	
	ni=p[i].getN();
	//cout << "ni: "<<ni << " " << work.size()<<endl;
       
	assert(work.size() >= ni );
	// generates the random numbers
	this->rand_o.gaussian(work,ni);
	// diffusion step
	for(j=0;j<ni;j++)
	  {
	    //cout << work[j]<<endl;
	    p[i].getNoBC(j)+=sqrt(this->delta_tau)*work[j];
	  }
      }
  }
  
  // drift step
  
  void drift(all_particles_t & p,const grad_t & grad)
  {
    
    for(int i=0;i<p.getN();i++)
      {
	for(int j=0;j<p[i].getN();j++)
	  {
	    p[i].getNoBC(j)+=this->delta_tau*grad[i][j];
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

template<class comp>
qmcMover<comp>* buildQMCMover(comp* qmcObj)
{
  typedef typename comp::wave_t wave_t;
  // load the type of mover from the input file
  string moverType;
  qmcObj->main_input->reset()->get_child("method")->get_child("propagator");
  
  if ( qmcObj->main_input->check())
    {
      moverType=qmcObj->main_input->get_value()->get_string();
    }
  else
    {
      moverType="1order";
    }
  
  cout << "moverType: "<<moverType<<endl;
  
  // build the mover
  if (moverType=="1order")
    {
      return new qmcMover1Order<comp>(*qmcObj->rand,qmcObj->delta_tau);
    }
  
  if (moverType=="2order")
    {
      qmcMover2Order<comp>* ret;
      ret=new qmcMover2Order<comp>(*qmcObj->rand,qmcObj->delta_tau,qmcObj->wave,qmcObj->geo);
      ret->init(qmcObj->main_input);
      return ret;
    }
  
}

#endif
