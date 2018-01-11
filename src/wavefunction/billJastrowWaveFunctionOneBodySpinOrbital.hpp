
template<class jastrow_t,class comp>
class bill_jastrow_wavefunction_one_bodySpinOrbital : public bill_jastrow_wavefunction<jastrow_t,comp>
{
public:
  typedef typename wavefunction<comp>::particles_t particles_t;
  typedef typename wavefunction<comp>::qmc_t qmc_t;
  typedef typename wavefunction<comp>::all_particles_t all_particles_t;
  typedef typename wavefunction<comp>::value_t value_t;
  typedef typename wavefunction<comp>::grad_t grad_t;

  
  bill_jastrow_wavefunction_one_bodySpinOrbital(qmc_t* qmc_obj_,jastrow_t& jastrowO) : bill_jastrow_wavefunction< jastrow_t,comp>(qmc_obj_,jastrowO) {};
  
  virtual void laplacianMinusGradientSquared(const all_particles_t & p,grad_t & grad,value_t & e)
{
  typedef typename grad_t::gradParticles_t gradParticles_t;

   typedef typename grad_t::gradParticles_t gradParticles_t;
  value_t tmp,tmp1,tmp2;
  const particles_t &p1=p[this->src_particle_set];
  gradParticles_t & grad1=grad[this->src_particle_set];
  double x,sign,d;
  e=0;
  
  for(int i=0;i<p1.size();i++)
    {
      x=this->qmc_obj->geo->distance_pbc(p1[i].position(),this->jastrowc.getCenter(p1[i].spin()));
      d=abs(x);
      //returns the direction of the wavefunction
      sign=x/d;
      
      tmp1=this->jastrowc.d1d0(d,p1[i].spin());
      
      tmp2=this->jastrowc.d2d0(d,p1[i].spin());
      
      e=e+(tmp2 - tmp1*tmp1);
      grad1[i]+=sign*tmp1;
      
    }
 
  
}


  virtual void laplacianMinusGradientSquaredLogWave(const all_particles_t & p,grad_t & grad,value_t & e,value_t & waveValue)
{
  typedef typename grad_t::gradParticles_t gradParticles_t;
  
  typedef typename grad_t::gradParticles_t gradParticles_t;
  value_t tmp,tmp1,tmp2;
  const particles_t &p1=p[this->src_particle_set];
  gradParticles_t & grad1=grad[this->src_particle_set];
  double x,sign,d;
  e=0;
  waveValue=0;
  for(int i=0;i<p1.size();i++)
    {
      x=this->qmc_obj->geo->distance_pbc(p1[i].position(),this->jastrowc.getCenter(p1[i].spin()));
      d=abs(x);
      //returns the direction of the wavefunction
      sign=x/d;
      
      tmp=this->jastrowc.d0(d,p1[i].spin());
      tmp1=this->jastrowc.d1(d,p1[i].spin())/tmp;
      tmp2=this->jastrowc.d2(d,p1[i].spin())/tmp;
      
      e+=(tmp2 - tmp1*tmp1);
      grad1[i]+=sign*tmp1;
      waveValue+=log(tmp);
    }
 
  
}
  void gradient(all_particles_t & p,grad_t & grad)
  {
    
    typedef typename grad_t::gradParticles_t gradParticles_t;
    value_t tmp,tmp1,tmp2;
    const particles_t &p1=p[this->src_particle_set];
    gradParticles_t & grad1=grad[this->src_particle_set];
    double x,sign,d;  
  
    for(int i=0;i<p1.size();i++)
      {
	x=this->qmc_obj->geo->distance_pbc(p1[i].position(),this->jastrowc.getCenter(p1[i].spin()));
	d=abs(x);
	//returns the direction of the wavefunction
	sign=x/d;
	tmp1=this->jastrowc.d1d0(d,p1[i].spin());
	grad1[i]+=sign*tmp1;
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
	//x=(qmc_obj->geo->pbc(p1->position[i]));
	x=abs(this->qmc_obj->geo->distance_pbc(p1[i].position(),this->jastrowc.getCenter(p1[i].spin())));
	value=value+log(this->jastrowc.d0(x,p1[i].spin()));
      }
  
    this->log_wavefunction_value=value;
    return value;
};
  
  virtual double spinFlip(int nFlip,const all_particles_t &p)
  {
    const particles_t &p1=p[this->src_particle_set];
    assert(nFlip<p1.size());
    double x,xFlip;
    
    assert(nFlip<p1.size());
    
    x=abs(this->qmc_obj->geo->distance_pbc(p1[nFlip].position(),this->jastrowc.getCenter(p1[nFlip].spin())));

    xFlip=abs(this->qmc_obj->geo->distance_pbc(p1[nFlip].position(),this->jastrowc.getCenter(-p1[nFlip].spin())));
    
    return this->jastrowc.d0(xFlip,-p1[nFlip].spin())/this->jastrowc.d0(x,p1[nFlip].spin());
	
  }
  
  virtual bill_jastrow_wavefunction_one_bodySpinOrbital<jastrow_t,comp> * clone()
  {
    bill_jastrow_wavefunction_one_bodySpinOrbital<jastrow_t,comp> * wave2;
    
    wave2=new bill_jastrow_wavefunction_one_bodySpinOrbital(this->qmc_obj,this->jastrowc);
    this->copyTo(wave2);
    return wave2;
  }

  
private:
  
};

template<class jastrow_t,class comp>
typename comp::swave_t* buildBillJastrowSpinOneBody(comp* qmc_obj,xml_input* xml_wave ,const string &filename)
{
  int src_particle_set;
  int target_particle_set;
  string file1,file2,kind;
  bill_jastrow_wavefunction_one_bodySpinOrbital<jastrow_t,comp> * newWave;
  xmlNodePtr cur;
  cur=xml_wave->cur;

  src_particle_set=xml_wave->get_attribute("setA")->get_int();
  target_particle_set=xml_wave->get_attribute("setB")->get_int();
  xml_wave->get_child("jastrow");
  kind=xml_wave->get_attribute("kind")->get_string();
  file1=xml_wave->get_attribute("file")->get_string();
  file2=xml_wave->get_attribute("file1")->get_string();

  jastrow_t jastrow(file1,file2);

  newWave=new bill_jastrow_wavefunction_one_bodySpinOrbital<jastrow_t,comp>(qmc_obj,jastrow);

  newWave->setSrcSet(src_particle_set);
  newWave->setTargetSet(src_particle_set);

  xml_wave->cur=cur;

  return newWave;
  
}
