template<class comp>
template<class qmc_type>
walker<comp>::walker(qmc_type* qmc_obj)
{
  int i;
  vector<int> ns;
  state=new all_particles_t();
  state_tmp=new all_particles_t();
  state->init(qmc_obj->main_input);
  state_tmp->init(qmc_obj->main_input);
  ev=0;
  e=0;
  state->getNs(ns);
  
  particlesGradient.resize(ns);
  particlesGradientBackup.resize(ns);  
}

template<class comp>
template<class qmc_type>
void dmc_walker<comp>::set(qmc_type * qmc_obj)
{
  walker<comp>::set(qmc_obj);
  descendants=1;
}

template<class comp>
template<class qmc_t>
void walker<comp>::set(qmc_t* qmc_obj)
{
  
  wavefunction_value=qmc_obj->wave->log_evaluate(state);
  qmc_obj->wave->laplacianGradient(*state,e,e_f,particlesGradient);
  particlesGradientBackup.clone(particlesGradient);
  ev=qmc_obj->potential_obj->evaluate(state);
  e=e+ev;
  e_old=e;
  weight=0;
}


template<class comp>
void walker<comp>::set_particles(walker<comp>::all_particles_t * p)
{
  *state=*p;
}

template<class comp>
template<class qmc_type>
void walker<comp>::load(xml_input* xml_walkers,qmc_type* qmc_o)
{
  state->load(xml_walkers);
  qmc_o->geo->all_particles_pbc(state);
}

// export the correct number of measurements
template<class qt>
void walker_load_dynamic_measures(vector<measure_dynamic<qt> *> & md ,xml_input* xml_md)
{
  bool futureWalkers;
  int bins;
  double l_box;
  
  l_box=xml_md->reset()->get_child("system")->get_child("lBox")->get_value()->get_real();
  
  xml_md->reset()->get_child("measures")->get_first_child();
  
  while(xml_md->check())
    {
      // sets the number of future walkers
      if (xml_md->get_attribute("futureWalkers") != NULL)
	{
	  futureWalkers=xml_md->get_bool();
	}
      else
	{
	  futureWalkers=false;
	}

      if (xml_md->get_attribute("bins") != NULL)
	{
	  bins=xml_md->get_int();
	}
      else
	{
	  bins=0;
	}
      
      if(xml_md->get_name()=="center_of_mass_difference" and futureWalkers==true )
	{
	  md.push_back(new futureWalkerScalar<qt,scalarStorage<double> >() ); 
	}
      
      if(xml_md->get_name()=="center_of_mass_differenceSquared" and futureWalkers==true )
	{
	  md.push_back(new futureWalkerScalar<qt,scalarStorage<double> >() ); 
	}
      
      if(xml_md->get_name()=="static_structure_factor" and futureWalkers==true )
	{
	  
	  md.push_back(buildStructureFactorWalker<qt>(xml_md,l_box)); 
	}

       if(xml_md->get_name()=="pair_correlation" and futureWalkers==true )
	{
	  
	  md.push_back(buildPairCorrelationFuture<qt>(xml_md,l_box)); 
	}
      
      
      if (xml_md->get_name()=="winding_number")
	{
	  md.push_back(build_center_of_mass_w<qt>(xml_md));
 	}
      
      
     
      xml_md->get_next();
    }
}

// returns the weight for the system
template<class comp>
double dmc_walker<comp>::get_weight(const double e_ref,dmc<comp>*  dmc_obj)
{
  double barrier;
  
  barrier=(dmc_obj->delta_tau)*(
			      (this->e + this->e_old)/2.
			      - e_ref
				);
  return exp(-barrier);
}

template<class comp>
template<class qmc_type>
dmc_walker<comp>::dmc_walker(qmc_type* dmc_obj) : walker<comp>(dmc_obj)
{
  int i=0;
  descendants=1;
  md.resize(0);
  walker_load_dynamic_measures<dmc<comp> >(md,dmc_obj->main_input);
  //build_md(md);
  
}

template<class comp>
template<class walker2_t>
void walker<comp>::clone(walker2_t & w)
{
  int i;
  assert(state != NULL);
  assert(w.state != NULL);
  (*state)=(*(w.state));
  e=w.e;
  ev=w.ev;
  weight=w.weight;
  e_difference=w.e_difference;
  wavefunction_value=w.wavefunction_value;
  particlesGradient.clone(w.particlesGradient);
  
  
}
template<class comp>
dmc_walker_fixed_node<comp>& dmc_walker_fixed_node<comp>::operator=(dmc_walker_fixed_node<comp> & w)
{
  dmc_walker<comp>::operator=(w);
  this->sign=w.sign;
}
// dmc walker fixed phase
template<class comp>
dmc_walker_fixed_phase<comp>& dmc_walker_fixed_phase<comp>::operator=(dmc_walker_fixed_phase<comp> & w)
{
  dmc_walker<comp>::operator=(w);
  this->phase=w.phase; 
}

template<class comp>
dmc_walker<comp>& dmc_walker<comp>::operator=(dmc_walker<comp> & w)
{
  int i;
  
  this->clone(w);
  descendants=w.descendants;
  for(i=0;i<md.size();i++)
    {
      (*md[i])=(*w.md[i]);
    }
  
  return *this;
}

template<class comp>
template<class qmc_type>
void dmc_walker_fixed_node<comp>::update(qmc_type* dmc_obj)
{
  dmc_walker_fixed_node<comp>::update(dmc_obj);
  if(this->accept)
    {
      this->sign_current=dmc_obj->wave->get_sign();
    }
}

template<class comp>
template<class qmc_type>
void dmc_walker_fixed_phase<comp>::update(qmc_type* dmc_obj)
{
  dmc_walker_fixed_phase<comp>::update(dmc_obj);
  if(this->accept)
    {
      this->phase_current=dmc_obj->wave->get_phase();
    }
}

template<class comp>
template<class qmc_type>
void dmc_walker<comp>::update(qmc_type* dmc_obj)
{
  int i;
  double random_number;
  double barrier;
  double e_new;
  double old_wavefunction_value;

  
  (*this->state_tmp)=(*this->state);
  
  
  //state->print();
  //state_tmp->print();
  //cout << dmc_obj->geo->l_box<<endl;
  old_wavefunction_value=this->wavefunction_value;
  
  //state->move(dmc_obj->rand,dmc_obj->delta_tau);
  //dmc_obj->move2order(state,state_tmp);
  //move1order(this->state,dmc_obj);

  dmc_obj->timers[2]->start();
  
  if (dmc_obj->qmcMoverO->getOrder()==2)
    {
      // compute and store F(R2)
      dmc_obj->qmcMoverO->computeEffectiveDriftForce((*(this->state)),this->particlesGradient);
      this->particlesGradientBackup.clone(dmc_obj->qmcMoverO->getEffectiveDriftForce());
    }
  else if(dmc_obj->qmcMoverO->getOrder()==1)
    {
      // store F(R)
      this->particlesGradientBackup.clone(this->particlesGradient);
    }
  
  dmc_obj->timers[2]->stop();
  
  dmc_obj->timers[3]->start();
  dmc_obj->qmcMoverO->move( (*(this->state)),this->particlesGradient);
  dmc_obj->timers[3]->stop();
  
  //state->particle_sets[1]->position_no_pbc[0]=0;
  dmc_obj->geo->all_particles_pbc(this->state);
  
  // update the wavefunction
  dmc_obj->timers[5]->start();
  this->wavefunction_value=dmc_obj->wave->log_evaluate(this->state);
  dmc_obj->timers[5]->stop();

  dmc_obj->timers[4]->start();
  // update the kinetic energy
  dmc_obj->wave->laplacianGradient(*(this->state),e_new,this->e_f,this->particlesGradient);
  dmc_obj->timers[4]->stop();
  // set the quantum drift force for the particles
  //dmc_obj->wave->set_drift_force(state);
  
  // compute the non-symmetric correction to metropolis
 
   dmc_obj->timers[2]->start();
  ++(dmc_obj->n_metropolis);
  
  //cout << 2*(wavefunction_value - old_wavefunction_value) + Q_probability(state,state_tmp,dmc_obj->delta_tau) <<endl;

  
  if (dmc_obj->qmcMoverO->getOrder()==2)
    {
      // compute F(R2 prime)
      dmc_obj->qmcMoverO->computeEffectiveDriftForce((*(this->state)),this->particlesGradient);

      this->accept=metropolis(2*(this->wavefunction_value - old_wavefunction_value)+ Q_probability( dmc_obj->qmcMoverO->getEffectiveDriftForce()  ,this->particlesGradientBackup,*(this->state),*(this->state_tmp),dmc_obj->delta_tau),dmc_obj->rand);
    }
  else
  if (dmc_obj->qmcMoverO->getOrder()==1)
    {
      
      
      this->accept=metropolis(2*(this->wavefunction_value - old_wavefunction_value)+ Q_probability(this->particlesGradient,this->particlesGradientBackup,*(this->state),*(this->state_tmp),dmc_obj->delta_tau),dmc_obj->rand);

      
    }

  
  dmc_obj->timers[2]->stop();
  if ( this->accept )
    {
      this->e_old=this->e;
      // print out the energy of the system
      this->e=e_new;
      ++(dmc_obj->success_metropolis);
      
      // compute the potential energy
      this->ev=dmc_obj->potential_obj->evaluate(this->state);
      this->e=e_new+this->ev;
      
    }
  else
    {
      this->wavefunction_value=old_wavefunction_value;
      (*this->state)=(*this->state_tmp);
      this->particlesGradient.clone(this->particlesGradientBackup);
      this->e_old=this->e;
    }
  
}
//sets the number of descendants
template<class comp>
template<class qmc_type>
void dmc_walker<comp>::set_descendants(double &e_ref,qmc_type* dmc_obj)
{
  double barrier;
  double random_number;
  // computes the barrier
  barrier=(dmc_obj->delta_tau)*(
			      (this->e + this->e_old)/2.
			      - e_ref
			      );
  
  random_number=dmc_obj->rand->uniform();
      // computs the numbers of descendatns
  if (abs(barrier) > 30)
    {
      descendants=0;
    }
  else
    {
      descendants=(int)(exp(-barrier) + random_number);
    }
  //cout <<"a:"<< accept<<endl;
  if (!this->accept)
    {
      descendants=1;
    }
  
  
}


template<class comp>
template<class qmc_type>
void dmc_walker_fixed_node<comp>::set_descendants(double &e_ref,qmc_type* dmc_obj)
{
  if (sign_current==sign)
    {
      dmc_walker<comp>::set_descendants(e_ref,dmc_obj);
    }
  else
    {
      this->descendants=0;
    }
}

template<class comp>
template<class qmc_type>
void dmc_walker_fixed_phase<comp>::set_descendants(double &e_ref,qmc_type* dmc_obj)
{
  if (phase_current==phase)
    {
      dmc_walker<comp>::set_descendants(e_ref,dmc_obj);
    }
  else
    {
      this->descendants=0;
    }
}

// performs some kind of measurement
template<class comp>
template<class measure_kind,class wavefunction_type>
void dmc_walker<comp>::make_measurements(measure_kind* m,wavefunction_type* wave,double current_step)
{
  m->make_measurements(this,wave);
}

template<class comp>
void walker<comp>::printGeneral()
{
  cout <<"w--------------"<<endl;
  cout<<"pos: ";
  this->state->print();
  cout<<endl;
  cout<<"e: "<<this->e<<endl;
  
}

template<class comp>
void walkers<comp>::print()
{
  int i;
  cout<<"------------walkers("<<n<<")---------------"<<endl;
  for(i=0;i<n;i++)
    { 
      ws[i]->print();
    }
}
// saves the walkers to the main file system


//packs some walkers in a certain pack
template<class comp>
void walker<comp>::pack(packed_data* walker_pack)
{
  walker_pack->pack(wavefunction_value,1);
  walker_pack->pack(e,1);  
  state->pack(walker_pack);
  particlesGradient.pack(walker_pack);
  
}

template<class comp>
int walker<comp>::get_pack_size()
{ 
  return pTools::get_pack_size(wavefunction_value) + pTools::get_pack_size(e) + state->get_pack_size() + particlesGradient.getPackSize();
  
}

template<class comp>
void dmc_walker<comp>::pack(packed_data* walker_pack)
{
  walker_pack->pack(this->descendants,1);
  walker<comp>::pack(walker_pack);
  for(int i=0;i<md.size();i++)
    {
      md[i]->pack(walker_pack);
    }
}
template<class comp>
int dmc_walker<comp>::get_pack_size()
{
  int size=0;
  
  size= walker<comp>::get_pack_size() + pTools::get_pack_size(this->descendants);
  
  for(int i=0;i<md.size();i++)
    {
      size+=md[i]->get_pack_size();
    }
  return size;
  
}

// unpack some walkers in a certain pack
template<class comp>
void walker<comp>::unpack(packed_data* walker_pack)
{ 
  walker_pack->unpack(wavefunction_value,1);
  walker_pack->unpack(e,1);
  state->unpack(walker_pack);
  particlesGradient.unpack(walker_pack);
}

template<class comp>
void dmc_walker<comp>::unpack(packed_data* walker_pack)
{
  walker_pack->unpack(this->descendants,1);
  walker<comp>::unpack(walker_pack);
  for(int i=0;i<md.size();i++)
    {
      md[i]->unpack(walker_pack);
    }
}

// saves xml elements
template<class comp>
void walker<comp>::save(xml_input* xml_walkers)
{
  state->save(xml_walkers);
 
}
