
template<class jastrow_t,class comp>
class bill_jastrow_wavefunction_two_body_symmetricSpinOrbital : public bill_jastrow_wavefunction<jastrow_t,comp>
{
public:
  typedef typename wavefunction<comp>::particles_t particles_t;
  typedef typename wavefunction<comp>::qmc_t qmc_t;
  typedef typename wavefunction<comp>::all_particles_t all_particles_t;
  typedef typename wavefunction<comp>::value_t value_t;
  typedef typename wavefunction<comp>::grad_t grad_t;
  
  bill_jastrow_wavefunction_two_body_symmetricSpinOrbital(qmc_t* qmc_obj_,jastrow_t& jastrowO) : bill_jastrow_wavefunction< jastrow_t,comp>(qmc_obj_,jastrowO) {};
  
  virtual void laplacianMinusGradientSquared(const all_particles_t & p,grad_t & grad,value_t & e)
{
  typedef typename grad_t::gradParticles_t gradParticles_t;
  
  value_t tmp,tmp1,tmp2;
  const particles_t &p1=p[this->src_particle_set];
  gradParticles_t & grad1=grad[this->src_particle_set];
  double x,sign,d;
  
  // reset to zero output
  e=0;
  
  for(int i=0;i<p1.size();++i)
    {
      for(int j=0;j<i;++j)
	{
	  x=this->qmc_obj->geo->distance_pbc(p1[i].position(),p1[j].position());
	  d=abs(x);
	      
	  sign=x/d;
	  tmp=this->jastrowc.d0(d,p1[i].spin(),p1[j].spin());
	  tmp1=this->jastrowc.d1(d,p1[i].spin(),p1[j].spin());
	  tmp2=this->jastrowc.d2(d,p1[i].spin(),p1[j].spin());
	      
	  e+=2*(tmp2/tmp - pow((tmp1/tmp),2));
	  grad1[i]+=sign*tmp1/tmp;
	  grad1[j]-=sign*tmp1/tmp;
	  
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
	
	    value=value+log(this->jastrowc.d0(x,p1[i].spin(),p1[j].spin()));
	  }
    
      }
   this->log_wavefunction_value=value;
   return value;
  }
  
  virtual double spinFlip(int nFlip,const all_particles_t &p)
  {
    const particles_t &p1=p[this->src_particle_set];
    assert(nFlip<p1.size());
    double value=0;
    double x;
    
    for (int i=0;i<p1.size();i++)
      {
	if (i!=nFlip)
	  {
	    x=abs(this->qmc_obj->geo->distance_pbc(p1[i].position(),p1[nFlip].position()));
	
	    value=value+log(this->jastrowc.d0(x,p1[i].spin(),-p1[nFlip].spin()))-log(this->jastrowc.d0(x,p1[i].spin(),p1[nFlip].spin()));
	  }
	
      }
    
    return exp(value);
  }
  
  virtual bill_jastrow_wavefunction_two_body_symmetricSpinOrbital<jastrow_t,comp> * clone()
  {
    bill_jastrow_wavefunction_two_body_symmetricSpinOrbital<jastrow_t,comp> * wave2;
    
    wave2=new bill_jastrow_wavefunction_two_body_symmetricSpinOrbital(this->qmc_obj,this->jastrowc);
    this->copyTo(wave2);
    return wave2;
  }
  
private:
    
};
    



template<class jastrow_t,class comp>
typename comp::swave_t* buildBillJastrowSpinTwoBody(comp* qmc_obj,xml_input* xml_wave ,const string &filename)
{
  int src_particle_set;
  int target_particle_set;
  string file1,file2,kind;
  bill_jastrow_wavefunction_two_body_symmetricSpinOrbital<jastrow_t,comp> * newWave;
  xmlNodePtr cur;
  cur=xml_wave->cur;

  src_particle_set=xml_wave->get_attribute("setA")->get_int();
  target_particle_set=xml_wave->get_attribute("setB")->get_int();
  xml_wave->get_child("jastrow");
  kind=xml_wave->get_attribute("kind")->get_string();
  file1=xml_wave->get_attribute("file")->get_string();
  file2=xml_wave->get_attribute("file1")->get_string();

  jastrow_t jastrow(file1,file2);

  newWave=new bill_jastrow_wavefunction_two_body_symmetricSpinOrbital<jastrow_t,comp>(qmc_obj,jastrow);

  newWave->setSrcSet(src_particle_set);
  newWave->setTargetSet(target_particle_set);

  xml_wave->cur=cur;

  return newWave;
  
}
