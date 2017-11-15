template<class comp>
template<class qmc_type>

#include "observables/optimize.h"
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
  mO.accumulate(w);
  
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
  xmlNodePtr cur;
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
      setOptimize(this->main_input->get_bool());
      
    }
  else
    {
      setOptimize(false);
    }
  
  if ( isOptimize() )
    {
      
      optimizePlan vmcOptplan;
      vector<double> parameters;
      statusGeneralizedEigenValue=0;
      
      vmcOptplan=buildOptimizePlan(this->getInputFileName());
      this->wave->getParameters(vmcOptplan,parameters);
      mO.setPlan(vmcOptplan ,this->wave);
      this->wave->getParameters(vmcOptplan,parameters);
      
      mO.setGradient(buildAllParticlesGradient1D(this->getInputFileName()));
      
      optParameters.resize(mO.getNParams());
      parametersProposal.resize(3);
      for(int k=0;k<parametersProposal.size();k++)
	{
	  parametersProposal[k].resize(optParameters.size());
	}
      
      // build the correlated energy estimator
      mEnergyCorrelated=new mEnergyCorrelated_t();
      mEnergyCorrelated->setGradient(buildAllParticlesGradient1D(this->getInputFileName()));
      
      for(int k1=0;k1<3;k1++)
	{
	  mEnergyCorrelated->addWavefunction(new wave_t(this->wave));
	}
      
      warmupOptimizeSteps=this->main_input->reset()->get_child("warmupOptimize")->get_value()->get_int();

      correlatedEnergySteps=this->main_input->reset()->get_child("correlatedEnergySteps")->get_value()->get_int();

      unCorrelationSteps=this->main_input->reset()->get_child("unCorrelationSteps")->get_value()->get_int();

     
      this->main_input->reset()->get_child("absErrOptimization");
      if (this->main_input->check())
	{
	  optimizationMode=absErrMode;
	  absErrorLimit=this->main_input->get_value()->get_real();
	 
	}
      else
	{
	  absErrorLimit=0;
	  optimizationMode=countStepsMode;
	  this->main_input->reset();
	} 
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
  
  int status;
  
  mO.transfer(0);
  
  if (this->mpi_task == 0)
    {
      vector<double> step;
      vector<double> means;
      vector<double> errors;
      
      save();
      
      this->wave->getParameters(mO.getPlan(),optParameters);
      mO.getMeanError(means,errors);
      
      /*
	Prints out information on the last run
       */
      
      cout<<"----------------------------"<<endl;
      
      cout<<"step:"<<this->current_step<<endl;
      cout<<"ratio: "<<this->success_metropolis/this->n_metropolis<<endl;
      cout<<"energy: "<<means[means.size()-1]<<"+-"<< errors[errors.size()-1] << endl;
      printf("Parameters\n");
      tools::print(optParameters);
      
      statusGeneralizedEigenValue=mO.getStep(step);
      #ifdef VERBOSE
      mO.print();
      printf("Step\n");
      if (statusGeneralizedEigenValue==0)
	{
	  printf("Sucessfull stepping!\n");
	}
      
      if (statusGeneralizedEigenValue<0)
	{
	  printf("Failed stepping!\n");
	}
      
      if (statusGeneralizedEigenValue>0)
	{
	  printf("Complex eigenvalues!\n");
	}
      #endif 
      //tools::print(step);
      assert(step.size()==optParameters.size());
      
      cout<<"------------------------------------------"<<endl;
      tools::add(optParameters,step,parametersProposal[0]);
  
      status=mO.getStep(step,abs(means[means.size()-1])*0.1);
      if (status==0)
	{
	  tools::add(optParameters,step,parametersProposal[1]);
	}
      else
	{
	  parametersProposal[1]=parametersProposal[0];
	}
      
      status=mO.getStep(step,abs(means[means.size()-1])*100);
      if (status==0)
	{
	  tools::add(optParameters,step,parametersProposal[2]);
	}
      else
	{
	  parametersProposal[2]=parametersProposal[0];
	}
      
      #ifdef VERBOSE
      printf("Proposed parameters [%i] \n",this->mpi_task);
      tools::print(parametersProposal[0]);
      tools::print(parametersProposal[1]);
      tools::print(parametersProposal[2]);
      #endif
      
      if (optimizationMode==absErrMode)
	{
	  
	  if (errors[errors.size()-1]>absErrorLimit)
	    {
	      statusGeneralizedEigenValue=-1;
	    }
	}
    }

  pTools::broadcast(statusGeneralizedEigenValue,0);
  pTools::broadcast(parametersProposal[0],0);
  pTools::broadcast(parametersProposal[1],0);
  pTools::broadcast(parametersProposal[2],0);
  
  if (statusGeneralizedEigenValue==0)
    {
      mEnergyCorrelated->getWaves()[0]->setParameters(mO.getPlan(),parametersProposal[0]);
  
      mEnergyCorrelated->getWaves()[1]->setParameters(mO.getPlan(),parametersProposal[1]);
  
      mEnergyCorrelated->getWaves()[2]->setParameters(mO.getPlan(),parametersProposal[2]);
    }
  
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
  int i,j,k;
  
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
      do
	{
	  
	  for (j=0;j<this->stepsPerBlock;j++)
	    {
	      for(k=0;k<this->unCorrelationSteps;k++)
		{
		  warmup_step(); // performs a MC step
		}
	      optimizationStep();
	    }
	  
	  optimizationOut();
	}
      
      while(statusGeneralizedEigenValue!=0);
	
      // correlated energy measurements for stabilization
      for(j=0;j<correlatedEnergySteps;j++)
	{
	  for(k=0;k<this->unCorrelationSteps;k++)
	    {
	      warmup_step(); // performs a MC step
	    }
	  stabilizationStep();
	}
      
      stabilizationOut();
      chooseNextOptimizationParameters();
      
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



template<class comp>
void vmc<comp>::stabilizationStep()
{
  w->update(this);
  
  this->wave->laplacianGradient((*w->state),w->e,w->e_f,w->getParticlesGradient());
  
  w->ev=this->potential_obj->evaluate(w->state);
  w->e+=w->ev;
  
  mEnergyCorrelated->accumulate(w);
}

template<class comp>
void vmc<comp>::stabilizationOut()
{
  if(correlatedEnergySteps==0)
    {
      indexMinEnergyProposal=0;
    }
  else
    {
      mEnergyCorrelated->transfer(0);
      if (this->mpi_task == 0)
	{
	  indexMinEnergyProposal=mEnergyCorrelated->getMinCorrelatedEnergy();
          #ifdef VERBOSE
	  mEnergyCorrelated->print();
      
	  printf("Chosen Parameter: %i\n",indexMinEnergyProposal);
      
          #endif
      
    }
      
    }
  mEnergyCorrelated->reset();
  
}

template<class comp>
void vmc<comp>::chooseNextOptimizationParameters()
{
  int changeParameter;
  changeParameter=-1;
  if (this->mpi_task==0)
    {
      if (statusGeneralizedEigenValue==0)
	{
	  if(indexMinEnergyProposal<3)
	    {
	      ofstream f;
	      vector<double> mean;
	      f.open("optimization.dat",fstream::app);
	      mO.getMean(mean);

	      f<<mean[mean.size()-1]<< " ";
	      for(int k=0;k<optParameters.size();k++)
		{
		  f<<optParameters[k]<<" ";
		}
	      f<<endl;
	      f.close();
	      changeParameter=0;
	      
	    }

	  
	}
    }
  
  pTools::broadcast(changeParameter,0);
  pTools::broadcast(indexMinEnergyProposal,0);
  
  if(changeParameter==0)
    {
      optParameters=parametersProposal[indexMinEnergyProposal];
      this->wave->setParameters(mO.getPlan(),optParameters);
      
      mO.reset();
      mO.setParameters(optParameters);
      w->set(this);
      
    }
	     
}
