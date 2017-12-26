
template<class jastrow_t,class tm>
class bill_jastrow_wavefunction_two_body_symmetricSpinOrbital : public bill_jastrow_wavefunction<jastrow_t,tm>
{
  
  typedef typename wavefunction<comp>::particles_t particles_t;
  typedef typename wavefunction<comp>::qmc_t qmc_t;
  typedef typename wavefunction<comp>::all_particles_t all_particles_t;
  typedef typename wavefunction<comp>::value_t value_t;
  typedef typename wavefunction<comp>::grad_t grad_t;

  
  bill_jastrow_wavefunction_two_body_asymmetric(qmc_t* qmc_obj_,jastrow_t& jastrowO) : bill_jastrow_wavefunction< jastrow_t,comp>(qmc_obj_,jastrowO) {};
  
  void laplacianMinusGradientSquared(all_particles_t & p,grad_t & grad,value_t & e)
{
  typedef typename grad_t::gradParticles_t gradParticles_t;
  
  value_t tmp,tmp1,tmp2;
  const particles_t &p1=p[this->src_particle_set];
  gradParticles_t & grad1=grad[this->src_particle_set];
  double x,sign,d;
  
  // reset to zero output
  e=0;
  
  for(int i=0;i<p1.size();i++)
    {
      for(int j=0;j<i;j++)
	{
	  x=this->qmc_obj->geo->distance_pbc(p1[i].position(),p1[j].position());
	  d=abs(x);
	  if (d!=0)
	    {
	      sign=x/d;
	      tmp=this->jastrowc.d0(d,p1[i].spin(),p1[j].spin());
	      tmp1=this->jastrowc.d1(d,(d,p1[i].spin(),p1[j].spin()));
	      tmp2=this->jastrowc.d2(d,p1[i].spin(),p1[j].spin());
	      
	      e+=2*(tmp2/tmp - pow((tmp1/tmp),2));
	      grad1[i]+=sign*tmp1/tmp;
	      grad1[j]-=sign*tmp1/tmp;
	    }
	  else
	    {
	      e=1000000;
	    }
	  
	}
    }
}
  
  void gradient(all_particles_t & p,grad_t & grad)
  {

    typedef typename grad_t::gradParticles_t gradParticles_t;
  
    value_t tmp,tmp1;
    const particles_t &p1=p[this->src_particle_set];
    gradParticles_t & grad1=grad[this->src_particle_set];
    double x,sign,d;
  
  
    for(int i=0;i<p1.size();i++)
      {
	for(int j=0;j<i;j++)
	  {
	    x=this->qmc_obj->geo->distance_pbc(p1[i].position(),p1[j].position());
	    d=abs(x);
	    if (d!=0)
	      {
		sign=x/d;
		tmp=this->jastrowc.d0(d,p1[i].spin(),p1[j].spin());
		tmp1=this->jastrowc.d1(d,p1[i].spin(),p1[j].spin());
		grad1[i]+=sign*tmp1/tmp;
		grad1[j]-=sign*tmp1/tmp;
	      }
	    
	  }
      }
  }

  double  log_evaluate(all_particles_t * p)
  {
  
    int i,j;
    double x,tmp,value;
  
    particles_t& p1=(*p)[this->src_particle_set];
    value=0;
  
    // assert that the required jastrows are defined
    for (i=0;i<p1.size();i++)
      {
	
	for (j=0;j<i;j++)
	  {
	    x=abs(this->qmc_obj->geo->distance_pbc(p1[i].position(),p1[j].position()));
	
	    value=value+log(this->jastrowc.d0(x,p1[i],p1[j]));
	  }
    
      }
    this->log_wavefunction_value=value;
    return value;
  }
  
};

template<class J,class comp>
wavefunction* createBillJastrowTwoBodySymmetricSpinOrbit(comp* qmc_obj,xml_input* xml_wave ,const string &filename)
{
  
}
