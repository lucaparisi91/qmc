template<class comp>
class qmcMover1OrderSpin : public qmcMover<comp>
{
public:
  typedef typename qmcMover<comp>::all_particles_t all_particles_t;
  typedef typename qmcMover<comp>::grad_t grad_t;
  typedef typename all_particles_t::particles_t particles_t;
  
  
  
  qmcMover1OrderSpin(int seed,double delta_tau_)  : qmcMover<comp>(seed,delta_tau_), particleMover(seed,delta_tau_){tauSpin=0.005;omega=1;};
  
  qmcMover1OrderSpin(rand_t & rand_,double delta_tau_) : qmcMover<comp>(rand_,delta_tau_),particleMover(rand_,delta_tau_){tauSpin=0.005;omega=1;};
  
  void setSpinTimeStep(double tau_)
  {
    tauSpin=tau_;
  }

  double getSpinTimeStep() const {return tauSpin;}
  
  virtual int getOrder()
  {
    return 1;
  }
  
  void reserve(size_t n)
  {
    work.resize(n);
    particleMover.reserve(n);
  }
  
  void moveUniformSpin(all_particles_t & p)
  {
    this->getRandomEngine().uniform(work);
    int k=0;
    int newSpin;
    for(int i=0;i<p.size();i++)
      {
	for(int j=0;j<p[i].size();j++)
	  {
	    
	    assert(work[k]>=0 and work[k]<=1);
	    if (work[k]<tauSpin)
	      {
		p[i][j].spin()*=-1;
	      }
	    k++;
	  }
      }
    
  }
  
  void moveSpinRabi(particles_t & p,vector<double> &spinFlipsRatios)
  {
    double pSpinFlip,tmp;
    this->getRandomEngine().uniform(work);
    
    for(int i=0;i<p.size();i++)
      {
	tmp=omega*particleMover.getTimeStep()*spinFlipsRatios[i];
	pSpinFlip=sinh(tmp);
	pSpinFlip/=(pSpinFlip + cosh(tmp));
	
	if (work[i]<pSpinFlip)
	  {
	    p[i].spin()*=-1;
	  }	    
      }
      
  }
  
  void moveGaussian(all_particles_t & p)
  { 
    particleMover.moveGaussian(p);
    
    moveUniformSpin(p);
    
    
    //while(getMagnetization(p[0])<0);
    
  }
  
  virtual void move( all_particles_t & p, const grad_t & grad)
  {
    particleMover.move(p,grad);
    //moveUniformSpin(p);
  }
  
private:
  vector<double> work;
  qmcMover1Order<comp> particleMover;
  double tauSpin;
  double omega;
};
