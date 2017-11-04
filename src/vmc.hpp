template<class comp>
template<class qmc_type>
#include "observables/optimizationObservablesLinearMethod.h"
#include "observables/correlatedEnergyDifference.h"

void vmc_walker<comp>::update(qmc_type* vmc_obj)
{
  double old_wavefunction_value;
  // make a copy of the old configuration
  (*this->state_tmp)=(*this->state);
  
  //state->print();
  //state_tmp->print();
  //cout << dmc_obj->geo->l_box<<endl;
  
  old_wavefunction_value=this->wavefunction_value;
  
  vmc_obj->moveEngine.moveGaussian(*(this->state));
  
  vmc_obj->geo->all_particles_pbc(this->state);
  
  // update the wavefunction
  
  this->wavefunction_value=vmc_obj->wave->log_evaluate(this->state);
  
  ++(vmc_obj->n_metropolis);
  
  this->accept=metropolis(2*(this->wavefunction_value - old_wavefunction_value),vmc_obj->rand);
  
  if ( this->accept )
    {
      ++(vmc_obj->success_metropolis); 
    }
  else
    {
      this->wavefunction_value=old_wavefunction_value;
      (*this->state)=(*this->state_tmp);
      
    };
}

template<class comp>
template<class measure_kind,class wavefunction_type>
void vmc_walker<comp>::make_measurements(measure_kind* m,wavefunction_type* wave,double current_step)
{
  m->make_measurements(this,wave);
}

template<class comp>
void vmc<comp>::step()
{
  this->current_step++;
  //cout << current_step<< " "<<ws->n<<endl;
  m->increment();
  
  w->update(this);
  // actually compute the kinetic energy
  this->wave->laplacianGradient((*w->state),w->e,w->e_f,w->getParticlesGradient());
  
  w->ev=this->potential_obj->evaluate(w->state);
  w->e+=w->ev;
  // performs the measurement
  w->make_measurements(m,this->wave,this->current_step);
  m->record(this->current_step);
}

template<class comp>
void vmc<comp>::optimizationStep()
{
  this->current_step++;
  //cout << current_step<< " "<<ws->n<<endl;
  
  w->update(this);
  // actually compute the kinetic energy
  this->wave->laplacianGradient((*w->state),w->e,w->e_f,w->getParticlesGradient());
  
  w->ev=this->potential_obj->evaluate(w->state);
  w->e+=w->ev;
  // accumulate the necessary quantitities for optimization
  mO->make_measurement(w,this->wave);
  
}
template<class comp>
void vmc<comp>::correlatedEnergyOptimizationStep()
{
  this->current_step++;
  //cout << current_step<< " "<<ws->n<<endl;
  
  w->update(this);
  // actually compute the kinetic energy
  this->wave->laplacianGradient((*w->state),w->e,w->e_f,w->getParticlesGradient());
  
  w->ev=this->potential_obj->evaluate(w->state);
  w->e+=w->ev;
  // accumulate the necessary quantitities for optimization
  mEnergyCorrelated->make_measurement(w,this->wave);
  
}

// initial warming step
template<class comp>
void vmc<comp>::warmup_step()
{
  this->current_step++;
  //cout << current_step<< " "<<ws->n<<endl;
  m->increment();
  w->update(this);
}

template<class comp>
vmc<comp>::vmc() : qmc<comp>(), moveEngine(*this->rand,this->delta_tau)
{
  e_t=0;
  load_wavefunctions<vmc<comp> >(this->main_input,waves,this);
  wave=new total_wavefunction<vmc<comp> >(this);
  wave->link_wavefunctions(this->main_input,waves,"main");
  
  wave->print_jastrows();
  
  // generate the new measurment
  m=new measures_t("input.xml",this);
  w=new walker_t(this);
  // whatever I want to save different walkers as input to a DMC calculation
  if (this->main_input->reset()->get_child("method")->get_attribute("inputDmc") != NULL)
    {
      inputDmc=this->main_input->get_bool();
    }
  else
    {
      inputDmc=false;
    }
  // check whatever to optimize the wavefunction
  
  if (this->main_input->reset()->get_child("method")->get_attribute("optimize") != NULL)
    
    {
      optimize=this->main_input->get_bool();
    }
  
  if (optimize==true)
    {
      mO=buildOptimizationObservablesLinearMethod<walker_t,wave_t>(this->main_input,this->wave);

      
      mEnergyCorrelated=buildCorrelatedEnergy<walker_t,wave_t>(this->main_input,this->wave);
      
      params.resize(3);
      
      warmupOptimizeSteps=this->main_input->reset()->get_child("warmupOptimize")->get_value()->get_int();

      correlatedEnergySteps=this->main_input->reset()->get_child("correlatedEnergySteps")->get_value()->get_int();
      
    }
  
  

  
}



template<class comp>
void vmc<comp>::save()
{
  qmc<comp>::saveGeneralQmc();
  // save the walker
  xml_input* xml_save;
  xml_save= new xml_input;
  xml_save->new_doc("walkers");
  // saves some settings
  xml_save->reset()->add_child("walker","");
  //cout <<xml_save->reset()->get_name();
  w->save(xml_save);
  
  xml_save->save("walkers.xml");
  delete xml_save;
}
// add a walker to the configurations file

template<class comp>
void vmc<comp>::saveAddWalker()
{
  qmc<comp>::saveGeneralQmc();
  xml_input* xml_save;
  
  xml_save= new xml_input;
  if ( check_file_exists("walkers.xml"))
    {
      xml_save->open("walkers.xml");
    }
  else
    {
      xml_save->new_doc("walkers");
    }
  
  xml_save->reset()->add_child("walker","");
  //cout <<xml_save->reset()->get_name();
  w->save(xml_save);
  xml_save->save("walkers.xml");
  
}

template<class comp>
// load the last walker from the file
void vmc<comp>::load()
{
  qmc<comp>::loadGeneralQmc();
  
  xml_input* load_xml;
  load_xml=new xml_input;
  if (check_file_exists("walkers.xml"))
    {
      
      load_xml->open("walkers.xml");
      load_xml->get_child("walker");
      delete w;
      w=new walker_t(this);
      
      if (load_xml->check())
	{
	  w->load(load_xml,this);
	}
      else
	{
	  cout << "No walker to load."<<endl;
	  exit(1);
	}
    }
  else
    {
	  
      cout << "Generating a starting a position..."<<endl;
      // set a uniform distribution for the particles
      
      w->state->set_uniform(this->geo->l_box,this->rand);
      
    }
      
  // initialize properties of the walker
  w->set(this);
  // reserbe memory for the move engine
  moveEngine.reserve(w->state->getNTot());
  w->make_measurements(m,this->wave,this->current_step);
  
}

template<class comp>
void vmc<comp>::out()
{
  unsigned int i;
  int mpi_task;
  m->reduce();
  if (this->mpi_task == 0)
    {
      cout<<"----------------------------"<<endl;
      m->ms[0]->print();
      cout<<"step:"<<this->current_step<<endl;
      cout<<"ratio: "<<this->success_metropolis/this->n_metropolis<<endl;	  
      cout<<"------------------------------------------"<<endl;
      m->out();
      m->save();
      if (inputDmc)
	{
	  saveAddWalker();
	}
      else
	{
	  save();
	}
    }
  m->clear();
}


template<class comp>
void vmc<comp>::optimizationOut()
{
  
  unsigned int i;
  int j;
  int mpi_task;
  double param;
  double paramStep;
  double eOpt;
  ofstream f;
  vector<double> obs;
  mO->reduce(0);
  if (this->mpi_task == 0)
    {
      
      //paramStep=mO->estimateStep();
      
      //wave->setParameter(param+paramStep);
      
      //cout << "Param: "<<param<<endl;
      //cout << "paramStep: "<<paramStep<<endl;
      // print to file the optimization results
      save();
      param=mO->getParameter();
      eOpt=mO->getEnergy();
      
      f.open("optimization.dat",std::ios_base::app);
      
      f<< param <<" "<<eOpt<<endl;
      
      f.close();
      
      /*
	Prints out information on the last run
       */
      
      cout<<"----------------------------"<<endl;
      
      cout<<"step:"<<this->current_step<<endl;
      cout<<"ratio: "<<this->success_metropolis/this->n_metropolis<<endl;
      cout<<"energy: "<<eOpt<<endl;
      cout<<"param: "<<param<<endl;
      
      cout<<"------------------------------------------"<<endl;
      /* Predict the next optimal parameter */
      
      obs=mO->getMean();
      paramStep=optimizer.getParameterVariation(obs);
      
      if (( (param+paramStep) > 10) or (  (param+paramStep)<=0 ))
	{
	  paramStep=0;
	}
      
      
      
     
      /* Choose the next 3 set of possible optimal parameters  and broadcast to other jobs*/
      params[0]=param+ paramStep;
      params[1]=param + paramStep/10.;
      params[2]=param + paramStep/100.;
      
      
      
    }

  MPI_Bcast(&params[0],3,MPI_DOUBLE,0,MPI_COMM_WORLD);
  /* Set up correlated energy measurements for stabilization */
  
  cout << "params"<<endl;
  print_vector(params);
  mO->setParameter(params[0],mEnergyCorrelated->getWave(0));
  mO->setParameter(params[1],mEnergyCorrelated->getWave(1));
  mO->setParameter(params[2],mEnergyCorrelated->getWave(2));
  mO->clear();
  //mO->setParameter(param,wave);
  //mO->init();
  w->set(this);
  w->make_measurements(m,this->wave,this->current_step);
  
}

template<class comp>
void vmc<comp>::optimizationAfterCorrOut()
{
  
  unsigned int i;
  int j;
  int mpi_task;
  double param;
  ofstream f;
  vector<double> obs;
  int iMin;
  mEnergyCorrelated->reduce(0);
  
  if (this->mpi_task == 0)
    {
      
      //paramStep=mO->estimateStep();
      
      //wave->setParameter(param+paramStep);
      
      //cout << "Param: "<<param<<endl;
      //cout << "paramStep: "<<paramStep<<endl;
      // print to file the optimization results
      //f.open("optimization.dat",std::ios_base::app);
      
      //f << param <<" "<<eOpt<<endl;
      //f.close();

      
      obs=mEnergyCorrelated->getMean();
      cout << "Correlated energies"<<endl;
      print_vector(obs);
      
      iMin=minIndex(obs);
      param=params[iMin];
      cout << "param min index: "<<iMin<<endl;

      
    }
  MPI_Bcast(&param,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  mEnergyCorrelated->clear();
  mO->setParameter(param,wave);
  mO->init();
  w->set(this);
  w->make_measurements(m,this->wave,this->current_step);
  
}


template<class comp>
void vmc<comp>::run()
{
  if (optimize)
    {
     
      runOptimize();
    }
  else
    {
      runStandard();
    }
};

template<class comp>
void vmc<comp>::runStandard()
{
  int i,j;
  
  load();
  cout << "Starting..."<<endl;
  
  this->startTimers();
  // warmup
  for(i=0;i<this->warmupBlocks;i++)
    {
      for (j=0;j<this->stepsPerBlock;j++)
	{	  
	  warmup_step(); // performs a MC step
	}
    }
  
  // assumes termalized
  for(i=0;i<this->nBlocks;i++)
    {
      for (j=0;j<this->stepsPerBlock;j++)
	{
	  step(); // performs a MC step
	}
      // make outputs and reductions
      
      out();
      
    }
  
  this->stopTimers();
  this->printTimers();
  
};


#include "observables/optimizationMatrixLinearMethod.h"

template<class comp>
void vmc<comp>::runOptimize()
{
  int i,j;
  
  load();
  cout << "Starting Optimization run..."<<endl;
  
  this->startTimers();

  // initial warmup
  for(i=0;i<this->warmupBlocks;i++)
    {
      for (j=0;j<this->stepsPerBlock;j++)
	{	  
	  warmup_step(); // performs a MC step
	}
    }
  
  for(i=0;i<this->nBlocks;i++)
    { 
      // warmup
      for(j=0;j<warmupOptimizeSteps;j++)
	{
	  warmup_step(); // performs a MC step 
	}
      
      // assumes termalized
      for (j=0;j<this->stepsPerBlock;j++)
      	{
	  optimizationStep(); // performs a MC step 
      	}
      
      
      optimizationOut();
      // correlated energy measurements for stabilization
      for(j=0;j<correlatedEnergySteps;j++)
	{
	  correlatedEnergyOptimizationStep();
	}

      // choose the parameter with the minimum energy and reinitialize the wavefunction
      
      optimizationAfterCorrOut();
      
      this->success_metropolis=0;
      this->n_metropolis=0;
    }
  
  this->stopTimers();
  this->printTimers();
  
};

// prints the walker composition
template<class comp>
void vmc_walker<comp>::print()
{
  this->printGeneral();  
}

// template class walker<D1_t>;
// template class walker<spinor1D>;
// template class vmc_walker<D1_t>;
// template class vmc_walker<spinor1D>;
// template class vmc<D1_t>;
// template class vmc<spinor1D>;
