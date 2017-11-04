
template<class X,class comp>
void bill_jastrow_wavefunction_two_body_symmetric<X,comp>::laplacianMinusGradientSquared(const typename bill_jastrow_wavefunction_two_body_symmetric<X,comp>::all_particles_t & p,typename bill_jastrow_wavefunction_two_body_symmetric<X,comp>::grad_t & grad,typename bill_jastrow_wavefunction_two_body_symmetric<X,comp>::value_t & e)
{
  typedef typename grad_t::gradParticles_t gradParticles_t;
  
  value_t tmp,tmp1,tmp2;
  const particles_t &p1=p[this->src_particle_set];
  gradParticles_t & grad1=grad[this->src_particle_set];
  double x,sign,d;
  
  // reset to zero output
  e=0;
  
  for(int i=0;i<p1.getN();i++)
    {
      for(int j=0;j<i;j++)
	{
	  x=this->qmc_obj->geo->distance_pbc(p1[i],p1[j]);
	  d=abs(x);
	  if (d!=0)
	    {
	      sign=x/d;
	      tmp=this->jastrowc.d0(d);
	      tmp1=this->jastrowc.d1(d);
	      tmp2=this->jastrowc.d2(d);
	  
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
};

template<class X,class comp>
void bill_jastrow_wavefunction_two_body_symmetric<X,comp>::gradient(const typename bill_jastrow_wavefunction_two_body_symmetric<X,comp>::all_particles_t & p,typename bill_jastrow_wavefunction_two_body_symmetric<X,comp>::grad_t & grad)
{
  typedef typename grad_t::gradParticles_t gradParticles_t;
  
  value_t tmp,tmp1;
  const particles_t &p1=p[this->src_particle_set];
  gradParticles_t & grad1=grad[this->src_particle_set];
  double x,sign,d;
  
  
  for(int i=0;i<p1.getN();i++)
    {
      for(int j=0;j<i;j++)
	{
	  x=this->qmc_obj->geo->distance_pbc(p1[i],p1[j]);
	  d=abs(x);
	  if (d!=0)
	    {
	      sign=x/d;
	      tmp=this->jastrowc.d0(d);
	      tmp1=this->jastrowc.d1(d);
	      grad1[i]+=sign*tmp1/tmp;
	      grad1[j]-=sign*tmp1/tmp;
	    }
	  
       }
    }
};


