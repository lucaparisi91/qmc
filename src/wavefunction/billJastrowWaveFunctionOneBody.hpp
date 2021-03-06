template<class X,class comp>
void bill_jastrow_wavefunction_one_body<X,comp>::laplacianMinusGradientSquared(const typename bill_jastrow_wavefunction_one_body<X,comp>::all_particles_t & p,typename bill_jastrow_wavefunction_one_body<X,comp>::grad_t & grad,typename bill_jastrow_wavefunction_one_body<X,comp>::value_t & e)
{
  typedef typename grad_t::gradParticles_t gradParticles_t;
  value_t tmp,tmp1,tmp2;
  const particles_t &p1=p[this->src_particle_set];
  gradParticles_t & grad1=grad[this->src_particle_set];
  double x,sign,d;
  e=0;
  
  for(int i=0;i<p1.size();i++)
    {
      x=this->getGeometry()->distance_pbc(p1[i].position(),this->jastrowc.center);
      
      d=abs(x);
      //returns the direction of the wavefunction
      sign=x/d;
      
      tmp1=this->jastrowc.d1d0(d);
      tmp2=this->jastrowc.d2d0(d);
      
      e=e+(tmp2 - tmp1*tmp1);
      grad1[i]+=sign*tmp1;
      
    }
  
};


template<class X,class comp>
void bill_jastrow_wavefunction_one_body<X,comp>::laplacianMinusGradientSquaredLogWave(const typename bill_jastrow_wavefunction_one_body<X,comp>::all_particles_t & p,typename bill_jastrow_wavefunction_one_body<X,comp>::grad_t & grad,typename bill_jastrow_wavefunction_one_body<X,comp>::value_t & e,typename bill_jastrow_wavefunction_one_body<X,comp>::value_t & waveValue)
{
  typedef typename grad_t::gradParticles_t gradParticles_t;
  value_t tmp,tmp1,tmp2;
  const particles_t &p1=p[this->src_particle_set];
  gradParticles_t & grad1=grad[this->src_particle_set];
  double x,sign,d;
  e=0;
  waveValue=0;
  for(int i=0;i<p1.size();i++)
    {
      x=this->getGeometry()->distance_pbc(p1[i].position(),this->jastrowc.center);
      d=abs(x);
      //returns the direction of the wavefunction
      sign=x/d;
      
      tmp=this->jastrowc.d0(d);
      tmp1=this->jastrowc.d1(d);
      tmp2=this->jastrowc.d2(d);
      
      e=e+(tmp2/tmp - pow(tmp1/tmp,2));
      grad1[i]+=sign*tmp1/tmp;
      waveValue+=log(tmp);
      
    }
  
};


template<class X,class comp>
void bill_jastrow_wavefunction_one_body<X,comp>::gradient(const typename bill_jastrow_wavefunction_one_body<X,comp>::all_particles_t & p,typename bill_jastrow_wavefunction_one_body<X,comp>::grad_t & grad)
{
  typedef typename grad_t::gradParticles_t gradParticles_t;
  value_t tmp,tmp1,tmp2;
  const particles_t &p1=p[this->src_particle_set];
  gradParticles_t & grad1=grad[this->src_particle_set];
  double x,sign,d;  
  
  for(int i=0;i<p1.size();i++)
    {
      x=this->getGeometry()->distance_pbc(p1[i].position(),this->jastrowc.center);
      d=abs(x);
      //returns the direction of the wavefunction
      sign=x/d;
      tmp1=this->jastrowc.d1d0(d);
      grad1[i]+=sign*tmp1;
    }
  
};


