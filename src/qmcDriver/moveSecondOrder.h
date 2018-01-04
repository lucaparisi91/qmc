#ifndef MOVESECONDORDER_H
#define MOVESECONDORDER_H

#include <cassert>
#include <vector>
using namespace std;

template<class comp>
class qmcMover2Order : public qmcMover<comp>
{
public:
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename comp::rand_t rand_t;
  typedef typename comp::grad_t grad_t;
  typedef typename comp::wave_t wave_t;
  typedef typename comp::geometry_t geometry_t;
  qmcMover2Order(int seed,double delta_tau_,wave_t* wave_,geometry_t * geo_)  : qmcMover<comp>::qmcMover(seed,delta_tau_)
  {
    wave=wave_;
    geo=geo_;
  };
  
  qmcMover2Order(rand_t & rand_,double delta_tau_,wave_t * wave_, geometry_t * geo_) : qmcMover<comp>::qmcMover(rand_,delta_tau_)
  {
    wave=wave_;
    geo=geo_;
  };
  
  void reserve(size_t n)
  {
    work.resize(n);
  }
  
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

  
  void drift(all_particles_t & p,const grad_t & grad,double delta)
  {
    
    for(int i=0;i<p.size();++i)
      {
	for(int j=0;j<p[i].size();++j)
	  {
	    p[i][j].positionNoBC()+=delta*grad[i][j];
	  }
      }
  }
  
  virtual int getOrder()
  {
    return 2;
  }
  
  // initializes temporary gradient and particles object from inputfile
  
  void init(string filename)
  {
    buildAllOrbitals(filename,ptmp);
    ptmp.sizes(ns);
    gradTmp.resize(ns);
  }
  
  // compute the second order effective drift force
  void computeEffectiveDriftForce(all_particles_t & p,const grad_t & grad)
  {
    ptmp=p;
    //first drift
    drift(ptmp,grad,this->delta_tau/2.);
    geo->all_particles_pbc(ptmp);
    // second drift
    wave->gradient(ptmp,gradTmp);
    ptmp=p;
    gradTmp+=grad;
    drift(ptmp,gradTmp,this->delta_tau/4.);
    geo->all_particles_pbc(ptmp);
    
    wave->gradient(ptmp,gradTmp);
    
  }
  
  grad_t & getEffectiveDriftForce()
  {
    return gradTmp;
  }
  
  void move( all_particles_t & p, const grad_t & grad)
  {
    // third drift
    drift(p,gradTmp,this->delta_tau);
    moveGaussian(p);
  }
  
private:
  vector<int> ns;
  vector<double> work;
  grad_t gradTmp; // used to store intermediate drift forces
  wave_t * wave;
  geometry_t *geo;
  all_particles_t ptmp;
  
};

#endif
