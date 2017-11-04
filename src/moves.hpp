template<class all_particles_t,class rand_t>
void moveGaussian(all_particles_t* p,rand_t* rand_o,const double &delta_tau)
{
  typedef typename all_particles_t::particles_t particles_t;
  int i,j;
  particles_t* p1;
    // gaussian step
  for(i=0;i<p->particle_sets.size();i++)
    {
      p1=p->particle_sets[i];
      rand_o->gaussian(p1->work);
      for(j=0;j<p1->n;j++)
	{
	  p1->position_no_pbc[j]=p1->position_no_pbc[j] + sqrt(delta_tau)*(p1->work[j]);
	}
    }
}

template<class all_particles_t,class rand_t>
void GaussianMoveSpin(all_particles_t* p,rand_t * rand_o, const double & delta_tau )
{
  typedef typename all_particles_t::particles_t particles_t;
  int i,j;
  particles_t * p1;
  double s;
  for(i=0;i<p->particle_sets.size();i++)
    {
      p1=p->particle_sets[i];
      rand_o->gaussian(p1->work);
      
      for(j=0;j<p1->n;j++)
	{
	  // update the phase S
	  //---> gaussian
	  s=arg(p1->spinComp[j][0]);
	  
	  s+=sqrt(delta_tau)*(p1->work[j]);
	  //----->drift
	  
	  // update the components of the system
	  p1->spinComp[j][0]=std::polar(1.,s);
	  p1->spinComp[j][1]=polar(1.,-s);
	  
	}
    }  
}

