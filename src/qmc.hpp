#include "tools.h"
#include "random.h"
#include "measures.h"
#include "timer.h"
#include "input.h"

template<class qmc_t>
potential<qmc_t>* build_potential(string filename,qmc_t* qmcO)
{
  xml_input* main_input;
  main_input=new xml_input;
  string kind;
  cout << filename<<endl;
  main_input->open(filename)->reset()->get_child("system")->get_child("oneBodyPotential");
  
  if ( main_input->check())
    {
      
      kind=main_input->get_attribute("kind")->get_string();
      if (kind=="harmonic")
	{
	  return buildSpeciesHarmonicPotential<qmc_t>(main_input,qmcO->geo);
	}
      else if (kind=="rabiCoupling")
	{
	  return buildRabiPotential<qmc_t>(main_input,qmcO->geo,qmcO->wave);
	}
      else
	{
	  cout<<"Unkown potential;"<<endl;
	  exit(0);
	}
    }
  else
    {
      return new empty_potential<qmc_t>(qmcO->geo);
    }
  
  delete main_input;
  
}



template<class comp>
qmc<comp>::qmc()
{  
  int i;
  xml_save=new xml_input;
  inputFileName="input.xml";
  
  xml_save->new_doc("save.xml");
  // adds a child to the total possible number of walkers
  xml_save->add_child("n_walkers",int_to_string(0));
  xml_save->add_child("walkers_file","");
  xml_save->add_child("current_step","");
  main_input=new xml_input;
  
  main_input->open("input.xml");
  
  if (main_input == NULL)
    {
      cout<<"Could not parse the input file.";
      
      exit(1);
      
    }

  n_metropolis=0;
  success_metropolis=0;
  if ( main_input->reset()->get_child("measures")->get_attribute("jumps") != NULL )
    {
      jumps=main_input->get_int();
    }
  else
    {
      jumps=0;
    }
  
  if (main_input->reset()->get_child("measures")->get_attribute("skip") != NULL)
  
    {
      skip=main_input->get_int();
    }
  else
    {
      skip=0;
    }
  
  running=false;
  current_step=0;
  
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_task);
  
  geo=new geometry_t("input.xml");
  
  rand=new random1(123456 + mpi_task);
  
  c=new counter(jumps,skip);
  
  for(i=0;i<14;i++)
    {
      timers.push_back(new timer());
    }
  
  timers[0]->setLabel("Walker Update");
  timers[1]->setLabel("Walker Measurements");
  timers[2]->setLabel("Metropolis");
  timers[3]->setLabel("Move");
  timers[4]->setLabel("Energy and Gradient");
  timers[5]->setLabel("Wavefunction");
  timers[6]->setLabel("Branching");
  timers[7]->setLabel("Sending walkers");
  timers[8]->setLabel("Receiving walkers");
  timers[9]->setLabel("Reducing");
  timers[10]->setLabel("Step Time");
  timers[11]->setLabel("Out Time");
  timers[12]->setLabel("Step init");
  timers[13]->setLabel("Walker bookkeeping");
  
  
  
  
  
  
  
  
  stepsPerBlock=main_input->reset()->get_child("stepsPerBlock")->get_value()->get_real();
   
  main_input->reset()->get_child("warmupBlocks");
  
  if (this->main_input->check() )
     {
       warmupBlocks=this->main_input->get_value()->get_int();
     }
   else
     {
       warmupBlocks=0;
     }
   
   nBlocks=this->main_input->reset()->get_child("nBlocks")->get_value()->get_int();

   core_populations.resize(mpi_tasks);
   n_metropolis=0;
   success_metropolis=0;
   delta_tau=main_input->reset()->get_child("method")->get_child("delta_tau")->get_value()->get_real();

}


template<class comp>
void qmc<comp>::saveGeneralQmc()
{
  xml_save=new xml_input; // new xml input file
  xml_save->new_doc("qmc_save");// creates a new document to be saved
  xml_save->reset()->add_child("success_step",real_to_string(success_metropolis));
  xml_save->reset()->add_child("current_step",real_to_string(this->current_step));
  xml_save->reset();
  xml_save->save("save.xml");
  //p->save(xml_save);
  delete xml_save;
  
};


template<class comp>
void qmc<comp>::loadGeneralQmc()
{
  if (check_file_exists("save.xml") )
    {
      xml_input* load_xml;
      load_xml=new xml_input;
      load_xml->open("save.xml");
  
      current_step=load_xml->reset()->get_child("current_step")->get_value()->get_real();
      max_steps+=current_step;
      delete load_xml;
    }
  
  
}

template<class comp>
void qmc<comp>::startTimers()
{
  t_init=MPI_Wtime();
}

template<class comp>
void qmc<comp>::stopTimers()
{
  t_end=MPI_Wtime();
}

template<class comp>
void qmc<comp>::printTimers()
{
  int i;
  double Tmax;

  Tmax=t_end - t_init;
  cout << "Wall time: "<< Tmax<<endl;
  cout << "task: "<<mpi_task<<endl;
  for(i=0;i<timers.size();i++)
    {
      cout << i << " " << timers[i]->get_total_time()/Tmax * 100 <<"%  "<< timers[i]->getLabel()<<endl;
      
    }
}

// template class qmc<D1_t>;
// template class qmc<spinor1D>;
