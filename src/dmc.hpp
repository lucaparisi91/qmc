#include "timer.h"

template<class comp>
dmc<comp>::dmc() : qmc<comp>()
{  
  load_wavefunctions<dmc<comp> >(this->main_input,waves,this);
  wave=new total_wavefunction<dmc<comp> >(this);
  wave->link_wavefunctions(this->main_input,waves,"main");
  wave->print_jastrows();
  
  xml_walkers=new xml_input;
  detailed_balance=0;
  // import variables from an input file
  inBuffer=new char[1000000];
  
  mean_walkers=this->main_input->reset()->get_child("method")->get_child("mean_walkers")->get_value()->get_int();
  
  delta_walkers=this->main_input->reset()->get_child("method")->get_child("delta_walkers")->get_value()->get_int();
  
  // number of steps per block
  
  detailed_balance=true;
  
  if (
      this->main_input->reset()->get_child("method")->get_attribute("kind")->get_string() == "svmc"
      )
    {
      smart_vmc=true;
      //mean_walkers=1;
    }
  else
    {
      smart_vmc=false;
    }
  
  delta_tau_or=this->delta_tau;
  e_t=0;
  ws=new walkers<comp>(this);
  
  
  // create a dock, for mpi comunications
  dmc_dock=new dock();
  
  //qmcMoverO=new qmcMover1Order<dmc<comp> >(*this->rand,this->delta_tau);

  qmcMoverO=buildQMCMover(this);
  
}

template<class comp>
template<class qmc_t>
void walkers<comp>::branch(int &i,qmc_t* dmc_obj)
{
  int j;
  walker_t* w;
  w=ws[i];
  int n_descendants_backup=0;
  // avoid a collapse if the number of walkers goes to zero
  if (w->descendants==0 )
    {
      // kill the walker
      // copy the last walker here...
      if (n==1)
	      {
	      i=i+1;
	      w->descendants=1;
	      return;
	      }
           *(w)=*(ws[n-1]);
	   
	    // ... and destroy the last element by resetting the length variable n
	    n=n-1;
	    //cout <<"b:"<<n<<":";
    }
  else
    {
    i=i+1;
  if (w->descendants > 1)
    {
      n_descendants_backup=w->descendants;
      w->descendants=1;
      // copy the walkers
      for (j=0;j<(n_descendants_backup-1);j++)
	{
	  n=n+1;
	  //cout << "b:"<<n<<":";
	  if (n >= 10000)
	    {
	      cout << "Population has grown too much."<<endl;
	      exit(1);
	    }
	  if ( n <= static_cast<signed int>(ws.size()))
	    {
	      
	      *ws[n-1]=*w;
	    }
	    else
	      {
		
		ws.push_back( new walker_t(dmc_obj));
		*ws[n-1]=*w;
	      }
	  
	}
      
    }
    }
}

template<class comp>
template<class qmc_type>
walkers<comp>::walkers(qmc_type* dmc_obj)
{
  
  work.resize(dmc_obj->geo->l_box);
  n=0;
  m=new measures_t("input.xml",dmc_obj);
}
template<class comp>
template<class qmc_t>
void walkers<comp>::generate_random(int n_to_add,qmc_t* dmc_obj)
{
  // add n walkers to the system
  int i;
  double e_f;
  ws.resize(0);
  n=0;
  dmc_obj->e_t=0;
  for (i=0;i<n_to_add;i++)
    {
      n=n+1;
      ws.push_back(new walker_t(dmc_obj));
      // generate a gaussian initial dsitribution
      ws[n-1]->state->gaussian(dmc_obj->rand);
      // apply periodic boundary condition
      dmc_obj->geo->all_particles_pbc(ws[n-1]->state);
      ws[n-1]->descendants=1;
    }
  assert((n) > 0);
  
  
}
template<class comp>
template<class qmc_t>
void walkers<comp>::generate_all_to(int n_to_add,double pos,qmc_t* dmc_obj)
{
  // add n walkers to the system
  int i;
  double e_f;
  ws.resize(0);
  n=0;
  dmc_obj->e_t=0;
  for (i=0;i<n_to_add;i++)
    {
      n=n+1;
      ws.push_back(new walker_t(dmc_obj));
      // generate a gaussian initial dsitribution
      ws[n-1]->state->set_all_positions(pos);
      // apply periodic boundary condition
      dmc_obj->geo->all_particles_pbc(ws[n-1]->state);
      ws[n-1]->descendants=1;
      
    }
  assert((n) > 0);
  
}

template<class comp>
template<class qmc_t>
void walkers<comp>::generate_uniform(int n_to_add,double pos,random1* randg,qmc_t* dmc_obj)
{
  // add n walkers to the system
  int i;
  double e_f;
  ws.resize(0);
  n=0;
  dmc_obj->e_t=0;
  for (i=0;i<n_to_add;i++)
    {      
      n=n+1;
      ws.push_back(new walker_t(dmc_obj));
      // generate a gaussian initial dsitribution
      ws[n-1]->state->set_uniform(pos,randg);
      // apply periodic boundary condition
      
      dmc_obj->geo->all_particles_pbc(ws[n-1]->state);
      ws[n-1]->descendants=1;
    }
  assert((n) > 0);
}

template<class comp>
template<class qmc_t>
void walkers<comp>::generate_gaussian(int n_to_add,double sigma,double position,random1* randg,qmc_t* dmc_obj)
{
  // add n walkers to the system
  int i;
  double e_f;
  ws.resize(0);
  n=0;
  dmc_obj->e_t=0;
  for (i=0;i<n_to_add;i++)
    {      
      n=n+1;
      ws.push_back(new walker_t(dmc_obj));
      // generate a gaussian initial dsitribution
      ws[n-1]->state->set_gaussian(randg,sigma,position);
      // apply periodic boundary condition
      dmc_obj->geo->all_particles_pbc(ws[n-1]->state);
      ws[n-1]->descendants=1;
    }
  assert((n) > 0);
}


// saves the state of the DMC to a file
template<class comp>
void dmc<comp>::save()
{
  int i,j;
  ofstream fp;
  ofstream fw;
  string filename;
  
  qmc<comp>::saveGeneralQmc();
  filename="walkers.xml";
  ws->save(filename);  
}

template<class comp>
void dmc<comp>::saveParallel(const string & filename)
{
  // build the output string
  char* outBuffer;
  ofstream f;
  string outString;
  int * stringLengths;
  int j,position;
  int size;
  stringLengths=new int[this->mpi_tasks];
  MPI_Status recvStatus[this->mpi_tasks];
  outString=ws->saveToString(filename);
  outBuffer=(char *)outString.c_str();
  size=outString.size();
  
  // find out all lengths
  MPI_Gather(&size,1,MPI_INT,stringLengths,1,MPI_INT,0,MPI_COMM_WORLD);
  //
  j=0;
  if (this->mpi_task==0)
    {
      inBuffer=strncpy(inBuffer,outBuffer,stringLengths[0]);
      position=stringLengths[0];
      
      for(j=1;j<(this->mpi_tasks);j++)
	{
	  MPI_Recv(&inBuffer[position],stringLengths[j],MPI_CHAR,j,167,MPI_COMM_WORLD,&recvStatus[j]);
	  position+=stringLengths[j];
	}
      inBuffer[position]='\0';      
      outString="<?xml version='1.0' ?><walkers>"+string(inBuffer) + "</walkers>";
  
      f.open(filename.c_str());
 
      f<<outString<<endl;
  
      f.close();
    }
  else
    {
      MPI_Send(outBuffer,size,MPI_CHAR,0,167,MPI_COMM_WORLD);
    }
  
}

template<class comp>
void dmc<comp>::load()
{
  int i,j,n;
  int i_walker,i_walker_tmp,i_particle;
  double x;
  int n_walkers;
  stringstream s;
  string line;
  xml_input* load_xml;
  ifstream load_file;
  double e_f;
  string kind;
  
  qmc<comp>::loadGeneralQmc();
  
  if (this->mpi_task==0)
    {
  
      if ( !check_file_exists("walkers.xml") )
    {
      cout<<"No file to load. Generating a random initial positions."<<endl;
      this->main_input->reset()->get_child("system")->get_child("initialCondition");
      if ( this->main_input->cur == NULL)
	{
	  
	  // generate an uniform distribution around the system
	  ws->generate_uniform(mean_walkers,this->geo->l_box,this->rand,this);
	}
      else
	{
	  
	  kind=this->main_input->get_attribute("kind")->get_string();
	  if(kind=="uniform")
	    {
	      
	      double length;
	      if (this->main_input->get_attribute("length")==NULL)
		{
		  length=this->geo->l_box;
		}
	      else
		{
		  length=this->main_input->get_real();
		}
	      
	      ws->generate_uniform(mean_walkers,length,this->rand,this);
	    }
	  else
	    {
	      if (kind=="gauss")
		{
		  double position=this->main_input->get_attribute("position")->get_real();
		  double sigma=this->main_input->get_attribute("sigma")->get_real();
		  ws->generate_gaussian(mean_walkers,sigma,position,this->rand,this);
		  
		}
	    }
	  
	 
	}
      
    }
  else
    {
      // creates a new xml imput
      load_xml=new xml_input;
      load_xml->open("walkers.xml");
      load_xml->get_child("walker");
      ws->ws.resize(0);
      ws->n=0;
  // load the walkers in file up to mean_walkers
  do
    {
      ws->n=ws->n+1;
      ws->ws.push_back(new walker_t(this));
      
      ws->ws[ws->n-1]->load(load_xml,this);
      
      //cout<< load_xml->get_name()<<endl;
      load_xml->get_next("walker");
    }
  while (load_xml->check() and (ws->n <mean_walkers) );
  
  delete load_xml;
  n_walkers=ws->n;
  for(i=n_walkers+1;i<=mean_walkers;i++)
    {
      ws->n=ws->n + 1;
      ws->ws.push_back(new walker_t(this));
      *(ws->ws[ws->n-1])=*(ws->ws[(i-n_walkers)%n_walkers]);
      
    }
  
    }
  
  e_t=0;
  
  qmcMoverO->reserve(ws->ws[0]->state->getNTot());
  for(i=0;i<ws->n;i++)
    {
      // measure the kinetic energy
      
      //this->wave->potential(this->potential_obj,ws->ws[i]->state,ws->ws[i]->ev);
      
      //ws->ws[i]->e+=ws->ws[i]->ev;
      //wave->set_drift_force(ws->ws[i]->state);
      
      ws->ws[i]->set(this);
      //ws->ws[i]->make_measurements(ws->m,this->wave,this->current_step);
      e_t=e_t + ws->ws[i]->e;
      
    }
  
  //exit(1);
  e_t=e_t/ws->n;
  
  n_walkers_c=ws->n;
  
  // set the current population of the walkers
  dmc_dock->pack_size=ws->ws[0]->get_pack_size();
  
  }

 
  // send the size of each pack to all walkers
  MPI_Bcast(&(dmc_dock->pack_size),1,MPI_INT,0,MPI_COMM_WORLD);
  cout << "pack_size "<< dmc_dock->pack_size << endl;
  
  dmc_dock->send_walkers_async(ws);
  dmc_dock->receive_walkers_async(ws,this);
  dmc_dock->send_walkers_async(ws);
  dmc_dock->print_populations();

  qmcMoverO->reserve(ws->ws[0]->state->getNTot());
  //dmc_dock->exchange_populations(ws->n);
  //dmc_dock->send_energy(e_t);
  //dmc_dock->recv_energy();
  
}
template<class comp>
// prints the state of the walkers object
void dmc_walker<comp>::print()
{
  this->printGeneral();
  cout<<"descendants: "<<descendants<<endl;
  
}

template<class comp>
void dmc<comp>::step()
{
  int i;
  double e_t_tmp;
  int n_walkers2,n_walkers1;
  MPI_Status stat;
  
  MPI_Request populations_send_requests[this->mpi_tasks];
  MPI_Request populations_rec_requests[this->mpi_tasks];
  MPI_Status populations_stat_requests[this->mpi_tasks];

  double start,stop;
  
  this->current_step++;
  //cout << current_step<< " "<<ws->n<<endl;
  ws->m->increment();
  n_walkers1=ws->n;
  e_t_tmp=0;
  assert((ws->n) >0);
  n_walkers2=0;
  
  if (this->n_metropolis > 0 and (this->success_metropolis>0) )
    {
      this->delta_tau=this->delta_tau_or*this->success_metropolis/(1.*this->n_metropolis);
    }
  
  for (i_walker=0;i_walker< ws->n;i_walker++)
    {
      this->timers[0]->start();
      ws->ws[i_walker]->update(this);
      this->timers[0]->stop();

      this->timers[1]->start();
      ws->ws[i_walker]->make_measurements(ws->m,this->wave,this->current_step);
      this->timers[1]->stop();
    }
  
  this->timers[8]->start();
  //dmc_dock->receive_walkers_async(ws,this);
  this->timers[8]->stop();
  
  for (;i_walker< ws->n;i_walker++)
    {
      this->timers[0]->start();
      ws->ws[i_walker]->update(this);
      this->timers[0]->stop();

      this->timers[1]->start();
      ws->ws[i_walker]->make_measurements(ws->m,this->wave,this->current_step);
      this->timers[1]->stop();
      
    }
  
  e_t_tmp=0;
  total_walkers=0;
  
  for(i_walker=0;i_walker<ws->n;i_walker++)
    {
      ws->ws[i_walker]->set_descendants(e_t,this);
      // update the temporary e_t
      
      e_t_tmp=e_t_tmp + (ws->ws[i_walker]->e)*(ws->ws[i_walker]->descendants);

      // increase the total number of walkers
      total_walkers=total_walkers + (ws->ws[i_walker]->descendants);
    }
  i_walker=0;
  
  this->timers[6]->start();
  // do not branch in a VMC calculation
   if (!smart_vmc)
      {
	
	while(i_walker<ws->n)
	  {
	    ws->branch(i_walker,this);
	  }
	
	//assert(n_walkers2 == ws->n);
      }
   
   e_t_tmp=e_t_tmp/ws->n;
   this->timers[6]->stop();
  //timers[3]->start();
  //MPI_Barrier(MPI_COMM_WORLD);
  //timers[3]->stop();
  //cout<<"pop "<< (stop -start)<<endl;
   // commits to a record the value measured over this time step
   ws->m->record(this->current_step);
   // decide how many walkers are present
   //start=MPI_Wtime();
   //timers[4]->start();
   //cout << n_walkers1<<" "<<mpi_task<<" "<<current_step<<endl;
   
   //timers[4]->stop();
   //cout<<"send_receive_determine: "<< stop -start<<endl;
   //cout << "pop----------"<< ws->n<<" "<<mpi_task<<endl;
   //timers[5]->start();
   
   //timers[5]->stop();
   //cout<<"pack: "<< stop -start<<endl;
   
   //start=MPI_Wtime();
   //timers[6]->start();
   
   this->timers[9]->start();
   //MPI_Allreduce( &e_t_tmp,&e_t,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   this->timers[9]->stop();
   
   e_t=e_t/this->mpi_tasks;
   //e_t=e_t_tmp;
   this->timers[7]->start();
   //dmc_dock->send_walkers_async(ws);
   
   this->timers[7]->stop();
   //timers[6]->stop();
   //dmc_dock->energies[mpi_task]=e_t_tmp;
   //dmc_dock->send_energy(e_t_tmp);
   //dmc_dock->recv_energy();
   //stop=MPI_Wtime();
   //cout<<"reduction: "<< stop -start<<endl;
   total_walkers=dmc_dock->total_walkers();
   
   population_control();
   
   //if (mpi_task==0)
   //  {
   //    cout << current_step<<endl;
   //  }
   
   //assert(dmc_dock->core_populations[mpi_task]==n_walkers1);
   
  
   //cout<<"n_walkers:"<<mpi_task<<" "<<n_walkers1<<" "<< dmc_dock->core_populations[mpi_task]<<endl;
  
   //start=MPI_Wtime();
   //dmc_dock->exchange_populations(n_walkers2 - dmc_dock->received_array.size() + dmc_dock->sent_array.size());
  
   //stop=MPI_Wtime();
   //timers[7]->start();
   //dmc_dock->send_recv_packs();
   //timers[7]->stop();
   //cout << "send/recv: "<< stop - start<<endl;
}
// installs a population control check
template<class comp>
void dmc<comp>::warmup_step()
{
  int i;
  double e_t_tmp;
  int n_walkers2,n_walkers1;
  MPI_Status stat;
  
  MPI_Request populations_send_requests[this->mpi_tasks];
  MPI_Request populations_rec_requests[this->mpi_tasks];
  MPI_Status populations_stat_requests[this->mpi_tasks];

  double start,stop;
  
  this->current_step++;
  //cout << current_step<< " "<<ws->n<<endl;
  ws->m->increment();
  n_walkers1=ws->n;
  e_t_tmp=0;
  assert((ws->n) >0);
  n_walkers2=0;
  
  if (this->n_metropolis > 0)
    {
      this->delta_tau=this->delta_tau_or*this->success_metropolis/this->n_metropolis;
    }
  
  for (i_walker=0;i_walker< ws->n;i_walker++)
    {
      ws->ws[i_walker]->update(this);
      //ws->ws[i_walker]-make_measurements();
    }
  

  e_t_tmp=0;
  total_walkers=0;
  
  for(i_walker=0;i_walker<ws->n;i_walker++)
    {
      ws->ws[i_walker]->set_descendants(e_t,this);
      // update the temporary e_t
      
      e_t_tmp=e_t_tmp + (ws->ws[i_walker]->e)*(ws->ws[i_walker]->descendants);

      // increase the total number of walkers
      total_walkers=total_walkers + (ws->ws[i_walker]->descendants);
      
      
    }
  
  i_walker=0;
  
  // do not branch in a VMC like calculation
   if (!smart_vmc)
      {
	while(i_walker<ws->n)
	  {
	    ws->branch(i_walker,this);
	  }
	//assert(n_walkers2 == ws->n);
	
      }
  
   
   
   //e_t_tmp=e_t_tmp/n_walkers2;
   e_t=e_t_tmp/ws->n;
   
   //timers[3]->start();
   //MPI_Barrier(MPI_COMM_WORLD);
   //timers[3]->stop();
   //cout<<"pop "<< (stop -start)<<endl;
   // commits to a record the value measured over this time step
   ws->m->record(this->current_step);
   // decide how many walkers are present
   //start=MPI_Wtime();
   //timers[4]->start();
   //cout << n_walkers1<<" "<<mpi_task<<" "<<current_step<<endl;
   //dmc_dock->gather_populations(ws->n);
   //dmc_dock->send_receive_determine();
   //timers[4]->stop();
   //cout<<"send_receive_determine: "<< stop -start<<endl;
   //cout << "pop----------"<< ws->n<<" "<<mpi_task<<endl;
   //timers[5]->start();
   //dmc_dock->allocate_packs();
   //dmc_dock->pack(ws);
   //dmc_dock->isend_packs();
   //dmc_dock->irecv_packs();
   //timers[5]->stop();
   //cout<<"pack: "<< stop -start<<endl;
   
   //start=MPI_Wtime();
   //timers[6]->start();
   //MPI_Allreduce( &e_t_tmp,&e_t,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   //timers[6]->stop();
   //dmc_dock->energies[mpi_task]=e_t_tmp;
   //dmc_dock->send_energy(e_t_tmp);
   //dmc_dock->recv_energy();
   
   //stop=MPI_Wtime();
   //cout<<"reduction: "<< stop -start<<endl;
   //total_walkers=dmc_dock->total_walkers();
   //e_t=e_t/mpi_tasks;
   //dmc_dock->send_walkers_async(ws);
   //total_walkers=dmc_dock->get_walkers();
   population_control();
   
   //if (mpi_task==0)
   //  {
   //    cout << current_step<<endl;
   //  }
   
   //assert(dmc_dock->core_populations[mpi_task]==n_walkers1);
   
  
   //cout<<"n_walkers:"<<mpi_task<<" "<<n_walkers1<<" "<< dmc_dock->core_populations[mpi_task]<<endl;
  
   //start=MPI_Wtime();
   //dmc_dock->exchange_populations(n_walkers2 - dmc_dock->received_array.size() + dmc_dock->sent_array.size());
  
   //stop=MPI_Wtime();
   //timers[7]->start();
   //dmc_dock->send_recv_packs();
   //timers[7]->stop();
   //cout << "send/recv: "<< stop - start<<endl;
}

// installs a population control check
template<class comp>
void dmc<comp>::population_control()
{
  
  if (total_walkers > (mean_walkers + delta_walkers))
    {
  
      e_t=e_t + 1./this->delta_tau*( log(mean_walkers*1./(mean_walkers + delta_walkers)));
    }
  
  else if(total_walkers < (mean_walkers - delta_walkers))
    {
      e_t=e_t + 1./this->delta_tau*( log(mean_walkers*1./(mean_walkers - delta_walkers))); 
    }
  
}
// launches the main cycle
template<class comp>
void dmc<comp>::run()
{
  int i,j;
  double start,stop;
  
  // load input files(only the master)
  load();
  
  cout<<"Starting... "<<endl;
  this->startTimers();
  // warmup
  for(i=0;i<this->warmupBlocks;i++)
    {
      for (j=0;j<this->stepsPerBlock;j++)
	{	  
	  warmup_step(); // performs a MC step
	}
    }
  
  // assues termalized
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
}
template<class comp>
void dmc<comp>::out()
{
  unsigned int i;
  int mpi_task;
  ws->m->reduce();
  
  if (this->mpi_task == 0)
    {
      cout<<"----------------------------"<<endl;
      ws->m->ms[0]->print();
      cout<<"walkers: "<< total_walkers <<endl;
      cout<<"step:"<<this->current_step<<endl;
      
      cout<<"ratio: "<<this->success_metropolis/this->n_metropolis<<endl;
      
      cout<<"------------------------------------------"<<endl;
      
      ws->m->out();
      //ws->m->save();
      ws->log();
      save();
    }
  
  saveParallel("walkers.xml");
  ws->m->clear();
}

template<class comp>
void walkers<comp>::log()
{
  ofstream log;
  double w;
  int i;
  //log.open("walkers.log",ios::app);
  // men weight
  //for(i=0;i<n;i++)
  //  {
  //    w+=ws[i]->weight;
  //    log << dmc_obj->current_step<< " "<< n << " " <<ws[i]->weight<<endl;
  //  }
  //w=w/n;
  
  //log.close();
  
}

template<class qmc_t>
void move1order(typename qmc_t::all_particles_t* p,qmc_t* qmc_obj)
{
  typedef typename qmc_t::all_particles_t all_particles_t;
  typedef typename all_particles_t::particles_t particles_t;
  int i,j;
  particles_t* p1;
  
  // drift move
  for(i=0;i<p->particle_sets.size();i++)
    {
      p1=p->particle_sets[i];
      for(j=0;j<p1->n;j++)
	{
	  p1->position_no_pbc[j]=p1->position_no_pbc[j] + qmc_obj->delta_tau*(p1->drift_force[j]);
	  
	}
    }
  
  // makes the gaussian move
  moveGaussian(p,qmc_obj->rand,qmc_obj->delta_tau);
  
}

// move the spin as under the influcence of an effective mass

template<class qmc_t>
void MoveSpinAsFreeParticle(typename qmc_t::all_particles_t* p,qmc_t* qmc_obj)
{
  typedef typename qmc_t::all_particles_t all_particles_t;
  typedef typename all_particles_t::particles_t particles_t;
  int i,j;
  double s;
  particles_t * p1;
  
  for(i=0;i<p->particle_sets.size();i++)
    {
      p1=p->particle_sets[i];
      qmc_obj->rand->gaussian(p1->work);
      
      for(j=0;j<p1->n;j++)
	{
	  // update the phase S
	  //---> gaussian
	  s=arg(p1->spinComp[j][0]); 
	  s+=sqrt(qmc_obj->delta_tau)*(p1->work[j]);
	  //----->drift
	  s+=real_part(p1->spinDrift[j]) + qmc_obj->delta_tau*real_part(p1->spinDrift[j]);
	  // update the components of the system
	  p1->spinComp[j][0]=std::polar(1.,s);
	  p1->spinComp[j][1]=polar(1.,-s);
	  
	}
    }  
}


// displacement computed at second order

template<class all_particles_t,class qmc_t>
void move2order(all_particles_t* p, all_particles_t* ptmp,qmc_t* qmc_obj)
{
  typedef typename all_particles_t::particles_t particles_t;
  int i,j;
  particles_t* p1,*p2;
  //first drift step
  for(i=0;i<p->particle_sets.size();i++)
    {
      p1=p->particle_sets[i];
      for(j=0;j<p1->n;j++)
  	{
  	  p1->position_no_pbc[j]=p1->position_no_pbc[j] + qmc_obj->delta_tau/2*p1->drift_force[j];
	  
  	}
    }
  // second drift step
  qmc_obj->wave->set_drift_force(p);

  for(i=0;i<p->particle_sets.size();i++)
    {
      p1=p->particle_sets[i];
      p2=ptmp->particle_sets[i];
      for(j=0;j<p1->n;j++)
  	{
  	  p1->position_no_pbc[j]=p2->position_no_pbc[j] + qmc_obj->delta_tau/2*(p1->drift_force[j] + p2->drift_force[j])/2;
	  
  	}
    }
  //third drift step
  qmc_obj->wave->set_drift_force(p);
  
  for(i=0;i<p->particle_sets.size();i++)
    {
      p1=p->particle_sets[i];
      for(j=0;j<p1->n;j++)
	{
	  p1->position_no_pbc[j]=p2->position_no_pbc[j] + qmc_obj->delta_tau*p1->drift_force[j];
	  
	  
	}
    }
  
  // performs a gaussian time tep
  for(i=0;i<p->particle_sets.size();i++)
    {
      p1=p->particle_sets[i];
      qmc_obj->rand->gaussian(p1->drift_force);
      for(j=0;j<p1->n;j++)
	{
	  p1->position_no_pbc[j]=p1->position_no_pbc[j] + sqrt(qmc_obj->delta_tau)*p1->drift_force[j];
	  
	}
    }
}
// rotate the spinors [ according to the green function quantum monte carlo]
// apply e^{-alpha*\sigma_x}(a b)^t to the spinors

template<class all_particles_t>
void rabiDMCPropagate(all_particles_t* p,double &alpha)
{
  int i,j;
  typedef typename all_particles_t::particles_t particles_t;
  particles_t* p1;
  
  for(i=0;i<p->particle_sets.size();i++)
    {
      p1=p->particle_sets[i];
      
      for(j=0;j<p1->n;j++)
	{
	  // first degree component
	  p1->spinComp[j][0]=cosh(alpha)*p1->spinComp[j][0] - sinh(alpha)*p1->spinComp[j][1];

	  p1->spinComp[j][1]=cosh(alpha)*p1->spinComp[j][1] - sinh(alpha)*p1->spinComp[j][0];
	  
	} 
    }
}




template<class comp>
void walkers<comp>::save(string & filename )
{
  int i;
  xml_input* xml_save;
  xml_save= new xml_input;
  xml_save->new_doc("walkers");
  // saves some settings
  for(i=0;i<n;i++)
    {
      xml_save->reset()->add_child("walker","");
      //cout <<xml_save->reset()->get_name();
      ws[i]->save(xml_save);
    }
  
  //cout << xml_save->reset()->toString()<<endl;
  //cout << "---------"<<endl;
  xml_save->save(filename.c_str());
  delete xml_save;
  
}
// save particles to string
template<class comp>
string walkers<comp>::saveToString(const string &filename)
{
  int i;
  xml_input* xml_save;
  string out;
  xml_save= new xml_input;
  xml_save->new_doc("walkers");
  // saves some settings
  out="";
  for(i=0;i<n;i++)
    {
      xml_save->reset()->add_child("walker","");
      //cout <<xml_save->reset()->get_name();
      ws[i]->save(xml_save);
      out=out +  xml_save->toString();
    }
  
  delete xml_save;
  return out;
}

template<class comp>
void walkers<comp>::saveAppend(xml_input * xml_save)
{
  int i=0;
  for(i=0;i<n;i++)
    {
      xml_save->reset()->add_child("waler","");
      ws[i]->save(xml_save);
    }
}
