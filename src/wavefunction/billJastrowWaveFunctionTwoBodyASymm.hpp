
template<class X,class comp>
void bill_jastrow_wavefunction_two_body_asymmetric<X,comp>::laplacianMinusGradientSquared(const typename bill_jastrow_wavefunction_two_body_asymmetric<X,comp>::all_particles_t & p,typename bill_jastrow_wavefunction_two_body_asymmetric<X,comp>::grad_t & grad,typename bill_jastrow_wavefunction_two_body_asymmetric<X,comp>::value_t & e)
{
  // -----init
  typedef typename grad_t::gradParticles_t gradParticles_t;
  value_t tmp,tmp1,tmp2;
  const particles_t &p1=p[this->src_particle_set];
  const particles_t &p2=p[this->target_particle_set];
  gradParticles_t & grad1=grad[this->src_particle_set];
  gradParticles_t & grad2=grad[this->target_particle_set];
  double x,sign,d;
  int i,j;
  
  // computes the total average energy
  
  e=0;
  
  for(i=0;i<p1.getN();i++)
    {
      
      for(j=0;j<p2.getN();j++)
	{	
	  x=this->qmc_obj->geo->distance_pbc(p1[i],p2[j]);
	  d=abs(x);
	  if (d!=0)
	    {
	      
	      sign=x/d;
	      tmp=this->jastrowc.d0(d);
	      tmp1=this->jastrowc.d1(d);
	      tmp2=this->jastrowc.d2(d);
	  
	      grad1[i]+=sign*tmp1/tmp;
	      grad2[j]-=sign*tmp1/tmp;
	  
	      e+=2*(tmp2/tmp - pow((tmp1/tmp),2));
	    }
	  else
	    {
	      e=1000000;
	    }
	  
	}
    }
};

template<class X,class comp>
void bill_jastrow_wavefunction_two_body_asymmetric<X,comp>::gradient(const typename bill_jastrow_wavefunction_two_body_asymmetric<X,comp>::all_particles_t & p,typename bill_jastrow_wavefunction_two_body_asymmetric<X,comp>::grad_t & grad)
{
  // -----init
  typedef typename grad_t::gradParticles_t gradParticles_t;
  value_t tmp,tmp1,tmp2;
  const particles_t &p1=p[this->src_particle_set];
  const particles_t &p2=p[this->target_particle_set];
  gradParticles_t & grad1=grad[this->src_particle_set];
  gradParticles_t & grad2=grad[this->target_particle_set];
  double x,sign,d;
  int i,j;
  
  // computes the total average energy
  
  for(i=0;i<p1.getN();i++)
    {
      
      for(j=0;j<p2.getN();j++)
	{	
	  x=this->qmc_obj->geo->distance_pbc(p1[i],p2[j]);
	  d=abs(x);
	  
	  sign=x/d;
	  tmp=this->jastrowc.d0(d);
	  tmp1=this->jastrowc.d1(d);
	  
	  grad1[i]+=sign*tmp1/tmp;
	  grad2[j]-=sign*tmp1/tmp;
	  
	}
    }
};


