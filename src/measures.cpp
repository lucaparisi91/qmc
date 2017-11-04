#include "tools.h"
#include "input.h"
#include "qmc.h"
#include "geometry.h"
#include <fstream>
#include <cstring>
#include <cassert>
#include <cmath>
#include "mesh.h"
#include "xml-input.h"
#include "wavefunction.h"
#include "mpi.h"
#include "dmc.h"
#include "vmc.h"
#include "measures.h"

using namespace std;
// using a certain namespace
measure_scalar* build_measure_scalar(xml_input* xml_m,string label)
{
  int setA,setB,jumps,skip;
  bool history;
  measure_scalar* m;
  
  if (xml_m->get_attribute("setA") != NULL)
    {
      setA=xml_m->get_int();
    }
  else
    {      
      setA=0;
    }
  if (xml_m->get_attribute("setB") != NULL)
    {
      setB=xml_m->get_int();
    }
  else
    {
      setB=0;
    }
  if (xml_m->get_attribute("history") != NULL)
    {
      history=xml_m->get_bool();
    }
  else
    {
      history="false";
    }
  
  if (xml_m->get_attribute("jumps") != NULL)
    {
      jumps=xml_m->get_int();
    }
  else
    {
      jumps=0;
    }
  if (xml_m->get_attribute("skip") != NULL)
    {
      skip=xml_m->get_int();
    }
  
  m=new measure_scalar(label);
  m->set_a=setA;
  m->set_b=setB;
  m->set_history(history);
  m->set(jumps,skip);  
  return m;
}

// build a local(on the walker) measurement of the center of mass
measure_vector_mult_index* build_measure_vector_mult_index(xml_input* xml_m,string label)
{
  int setA,setB,jumps,skip,bins;
  bool history;
  int nblocks;
  measure_vector_mult_index * m;
  if (xml_m->get_attribute("setA") != NULL)
    {
      setA=xml_m->get_int();
    }
  else
    {
      
      setA=0;
    }
  
    
  if (xml_m->get_attribute("jumps") != NULL)
    {
      jumps=xml_m->get_int();
    }
  else
    {
      jumps=0;
    }
  if (xml_m->get_attribute("skip") != NULL)
    {
      skip=xml_m->get_int();
    }
  else
    {
      skip=0;
      
    }
  // gets the total number of bins
  if (xml_m->get_attribute("bins") != NULL)
    {
      bins=xml_m->get_int();
    }
  else
    {
      bins=100;
    }
  
  if (xml_m->get_attribute("nblocks") != NULL)
    {
      nblocks=xml_m->get_int();
    }
  else
    {
      nblocks=1;
    }
  
  //md_=new center_of_mass_w<double>(bins));
  //md_->set_a=setA;
  // if the walker is not the only one requiring un update
  
  //m=new time_difference< center_of_mass_w<double> >(bins,label,md_.size()-1);
  m=new measure_vector_mult_index(bins,label,history);
  m->set(jumps,skip);
  m->n_block_jump=nblocks;
  m->set_a=setA;
  return m; 
}
// // append a measurement performed on a walker
// void append_walker_measures(vector<measure_dynamic*> &md,vector<measure*> &ms,xml_input* xml_md,bool only_walker)
// {

  
//   xml_md->reset()->get_child("measures")->get_first_child();
//   // checks out the measure to be performed
//   while(xml_md->check())
//     {
//        if (xml_md->get_name()=="winding_number")
// 	{
	  
// 	  append_winding_number(xml_md,ms,md,"winding_number",only_walker);
// 	}
     
//       xml_md->get_next();
//     }
// }
// build the space vector measurement

// used to jump and skip iterations
// check()->
//        0 still waiting to initialize blocks
//   
counter::counter(int jumps_, int skip_)
{
  i=0;
  termal=false;
  set(jumps_,skip_);
  status=0;
}
void counter::set(int jumps_,int skip_)
{
  jumps=jumps_;
  skip=skip_;
}

int counter::check()
{
  return status;
}

void counter::increment()
{
  i=i+1;
  // checks if jumps steps have been done
  if (!termal)
    {
      
      if (i > jumps)
	{
	  
	  i=1;
	  termal=true;
	  status=1;
	}
    }
  // if jumps steps have been done
  if (termal)
    {
      if (i>skip)
	{
	  // mark status, after skipping 'skip' steps and reset status
	  status=2;
	  i=0;
	}
      else
	{
	  status=1;
	}
    }
}

void measure_scalar::record(double step_)
{
  
  if (history and n>0 and (c->check()==2))
    {
      
      //cout << sum_m << " " <<n<<" "<<c->check()<<endl;
      recorded_measurements.push_back(mean());
      recorded_times.push_back(step_);
      n_recorded=0;
      clear();
    }
}

void measure_vector::record(double step_)
{
  
}
void measure_scalar::clear_history()
{
  recorded_times.resize(0);
  recorded_measurements.resize(0);
  n_recorded=0;
}

void measure_scalar::out_history()
{
  unsigned int i;
  ofstream f;
  f.open((label + string(".dat")).c_str(),std::ios_base::app);
  for(i=0;i<recorded_measurements.size();i++)
    {
      f <<recorded_measurements[i]<<" "<<recorded_times[i]<<endl;
    }
  clear_history();
  f.close();
}
// performs a scalar measurement
measure_scalar::measure_scalar(string label_)
{
  n=0;
  sum_m=0;
  n_t=0;
  sum_m_t=0;
  label=label_;
  n_recorded=0;
  history=false;
  c=new counter(0,0);
  
}
// adds a new scalar value to the measurements
void measure_scalar::add(double m,double step)
{
  
  sum_m=sum_m + m;
  n=n+1;
  
  last_step=step;
  
}

// reset all measurements
void measure_scalar::clear()
{
  sum_m=0;
  n=0;
}

// prints out the measurement of the system
void measure_scalar::out()
{
  ofstream f;
  unsigned int i;
  
  f.open((label + string(".dat")).c_str(),std::ios_base::app);
  if (history==true)
	{
	  if (recorded_measurements.size() > 0)
	    {
	      for(i=0;i<recorded_measurements.size();i++)
		{
		  f << recorded_times[i]<< " "<< recorded_measurements[i] <<endl;
		}
	      
	    }
	}
  else
    {
      if (n>0)
	{
	  f << last_step << " "<< sum_m/n<<endl;
	}
    }
      f.close();
      dec.out(label);
      out_t();
}

void measure_scalar::out_t()
{
  ofstream f;
  
  if (n_t > 0)
    {
      f.open((label + string("_t.dat")).c_str());
      f << sum_m_t/n_t<<" "<<dec.get_error()<<" "<<is_converged()<<endl;
      f.close();
      
    }
  
}

void measure_scalar::reduce(int root)
{
  unsigned int i;
  int tasks=0;
  int task=0;
  
  double sum_m2;
  int n2;
  double mean_sum;
  vector<double> recorded_measurements2;
  vector<double> recorded_times2;
  
  MPI_Comm_size(MPI_COMM_WORLD,&tasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&task);
  
  recorded_measurements2.resize(recorded_measurements.size());
  recorded_times2.resize(recorded_times.size());
  n2=0;
  //cout << recorded_measurements.size();
  
  
  //cout << ";"<< recorded_measurements.size()<<";";
  
  if (history)
    {
      assert(recorded_measurements.size() > 0);
      MPI_Reduce(&recorded_measurements.front(),&recorded_measurements2.front(), recorded_measurements.size(), MPI_DOUBLE, MPI_SUM, root,MPI_COMM_WORLD);
    }
  else
    {
      
      MPI_Reduce(&sum_m, &sum_m2, 1,
	     MPI_DOUBLE, MPI_SUM, root,
	     MPI_COMM_WORLD);
      
      MPI_Reduce(&n,&n2,1,MPI_INT,MPI_SUM,root,MPI_COMM_WORLD);
      
    }
  
  if (task==root)
    {
      
      if (history)
	{
	  for(i=0;i<recorded_measurements.size();i++)
	    {
	      recorded_measurements[i]=recorded_measurements2[i]/tasks;
	    }
	}
      else
	{
	  n=n2;
	  // if measurements were actually performed
	  if (n>0)
	    {
	      sum_m=sum_m2;
	      dec.add(sum_m2/n2);// add the measurements to be decorrelated
	      sum_m_t=sum_m_t + sum_m;
	      n_t=n_t + n;
	    }
	  
	}
	
      
    }
  
}

void measure_vector::reduce(int root)
{
  int tasks;
  int task;
  unsigned int i;
  int n2;
  
  vector<double> sum2,mean;
    
  MPI_Comm_size(MPI_COMM_WORLD,&tasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&task);
  
  sum2.resize(sum.size());
  mean.resize(sum.size());
  // compute the mean value of the measured vector
  
  MPI_Reduce(&sum.front(),&sum2.front(), sum.size(), MPI_DOUBLE, MPI_SUM, root,MPI_COMM_WORLD);
  
  MPI_Reduce(&n,&n2,1,MPI_INT,MPI_SUM,root,MPI_COMM_WORLD);
  
  if (task == root ) 
    {
      iBlock++;
      n_tot=n_tot + n2;
      n=n2;
      for(i=0;i<sum.size();i++)
	{
	  sum[i]=sum2[i];
	  
	  sum_tot[i]=sum_tot[i] + sum[i];
	  
	  sum_block[i]=sum_block[i]+sum[i]/n2;
	  if ( iBlock%n_block_jump==0)
	    {
	      sum2_tot[i]=sum2_tot[i]+pow(sum_block[i]/n_block_jump,2);
	      sum_block[i]=0;
	    }
	    
	}
      dec.add(sum,n);    
      
    }
  
}

void measure_vector_mult_index::reduce(int root)
{
  int tasks;
  int task;
  unsigned int i;
  vector<int> n2;
  vector<double> sum2;
  vector<double> ns2;
  MPI_Comm_size(MPI_COMM_WORLD,&tasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&task);
  
  sum2.resize(sum.size());
  ns2.resize(ns.size());
  // compute the mean value of the measured vector
  
  MPI_Reduce(&sum.front(),&sum2.front(), sum.size(), MPI_DOUBLE, MPI_SUM, root,MPI_COMM_WORLD);
  MPI_Reduce(&ns.front(),&ns2.front(),ns.size(),MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
  
  if (task == root ) 
    {
      for(i=0;i<sum.size();i++)
	{
	  if (ns2[i]>0)
	    {
	      ns[i]=ns2[i];
	      sum[i]=sum2[i];     
	      sum_tot[i]=sum_tot[i] + sum[i];
	      ns_tot[i]+=ns[i];
	      dec.add(sum[i]/ns[i],i);
	    }
	}
      
    }
  
  
}

// // reduce a vector with multiple indices
// void measure_vector_mult_index::reduce(int root)
// {
//   int tasks;
//   int task;
//   unsigned int i;
//   vector<int> n2;
  
//   vector<double> sum2,mean;
//   MPI_Comm_size(MPI_COMM_WORLD,&tasks);
//   MPI_Comm_rank(MPI_COMM_WORLD,&task);
  
//   sum2.resize(sum.size());
//   mean.resize(sum.size());
  
//   // compute the mean value of the measured vector
  
//   //for(i=0;i<sum.size();i++)
//   //  {
//   //    if (ns[i] > 0)
//   //	{
//   //	  mean[i]=sum[i]/ns[i];
//   //	  //cout << i << " " << mean[i]<<endl;
//   //	}
//   //    else
//   //	{
//   //	  mean[i]=0;
//   //	}
//   //  }
  
//   //MPI_Reduce(&mean.front(),&sum2.front(), mean.size(), MPI_DOUBLE, MPI_SUM, root,MPI_COMM_WORLD);
  
//   for(i=0;i<sum.size();i++)
//   {
//     if (ns[i]>0)
//       {
// 	ns_tot[i]=ns_tot[i] + ns[i];
	
// 	sum_tot[i]=sum_tot[i] + sum[i];
// 	dec.add(sum[i]/ns[i],i);
//       }
//     //sum[i]=0;
//     //ns[i]=0;
//   }
  
  
  
// }




measure_vector_mult_index::measure_vector_mult_index(int len_,string label_,bool history_) : measure_vector(len_,label_,history_) 
{
  int i;
  ns.resize(len_);
  ns_tot.resize(len_);
  ns_block.resize(len_);
  nb.resize(len_);
  sum_of_blocks.resize(len_);
  n_blocks=0;
  for(i=0;i<len_;i++)
    {
      ns[i]=0;
      ns_tot[i]=0;
      ns_block.resize(len_);
      sum_of_blocks[i]=0;
      nb[i]=0;
    }
}

// template<>
// measures<rqmc>::measures(rqmc* qmc_obj_)
// {
//   xml_input* main_input;
//   double n_particles;
//   bool history;
//   unsigned int i,id;
//   int bins;
//   int slices;
//   string label;
//   int set_a,set_b;
//   double max=0;
//   qmc_obj=qmc_obj_;
//   main_input=new xml_input;
//   main_input->open("input.xml");
//   bins=0;
//   id=0;
  
//   //skip=main_input->reset()->get_child("stepsPerBlock")->get_value()->get_int() - 1;
  
//   /*
//   main_input->reset()->get_child("measures")->get_attribute("jumps");  
//   if (main_input->check())
//     {
//       jumps=main_input->get_int();
//     }
//   else
//     {
//       jumps=0;
//     }
//   */
//   //jumps=0;
//   n_particles=qmc_obj->geo->l_box;
  
//   ms.push_back(
// 	       new rqmcMeasurement<wave_t,get_energy_dmc<walker_t,total_wavefunction> >( new get_energy_dmc<walker_t,total_wavefunction>("energy"),1)
// 	       );
//     ms.push_back(
// 	       new rqmcMeasurement<wave_t,get_f_energy_dmc<walker_t,total_wavefunction> >( new get_f_energy_dmc<walker_t,total_wavefunction>("force_energy"),1)
// 		 );
  
//   i=1;
  
//   main_input->reset()->get_child("measures")->get_first_child();
  
//   while(main_input->check())
//     {
      
//       i=i+1;
//       if (main_input->get_attribute("slices")!=NULL)
// 	{
// 	  slices=main_input->get_int();
// 	}
      
//       if (main_input->get_name() == "density")
// 	{
// 	  if (main_input->get_attribute("label") != NULL)
// 	    {
// 	      label=main_input->get_string();
// 	    }
// 	  else
// 	    {
// 	      label="density";
// 	    }
	  
// 	  // add a new density measurement
// 	  ms.push_back(
// 		       new rqmcMeasurement<wave_t, density<walker_t,total_wavefunction> >
// 		       (
// 			new density<walker_t,total_wavefunction>(build_space_vector(main_input,label,-n_particles/2,n_particles/2))
// 			,1 )
// 		       );
// 	}
//       if (main_input->get_name() == "pair_correlation")
// 	{
// 	  if (main_input->get_attribute("max") != NULL)
// 	    {
// 	      max=main_input->get_real();
// 	    }
// 	  else
// 	    {
// 	      max=n_particles/2;
// 	    }
// 	  if (main_input->get_attribute("label")!=NULL)
// 	    {
// 	      label=main_input->get_string();
// 	    }
// 	  else
// 	    {
// 	      label="pair_correlation";
// 	    }
// 	  // add a new pair correlation to the system
// 	  ms.push_back(
// 		       new rqmcMeasurement<wave_t, pair_correlation_m<walker_t,total_wavefunction>  >
// 		       (
// 		       new pair_correlation_m<walker_t,total_wavefunction>(build_space_vector(main_input,label,0,max))
// 		       ,slices)
// 			);
// 	  // gets the winding number of the system
// 	  if (main_input->get_name() == "winding_number")
// 	    {
// 	      cout << "Winding number not yet supported!"<<endl;
// 	      exit(0);
// 	    }
// 	}
      
//       // if (main_input->get_attribute("time_slices")!= NULL)
//       // 	{
//       // 	  slices=main_input->get_int();
//       // 	  ms[ms.size()-1]->set_slices(slices);
//       // 	}
      
      
//       main_input->get_next();
	
//     }
  
//   for(i=0;i<ms.size();i++)
//     {
//       ms[i]->load();
//     }
  
// }

// mean value of a scalar measurement 
double measure_scalar::mean()
{
assert(n>0);
  
  return sum_m/n;
 
}

double measure_scalar::mean_t()
{
  assert(n_t>0);
  return (sum_m_t/n_t);
}
// perform measurements on some vector
measure_vector::measure_vector(int len_,string label_,bool history_=true) : dec(len_)
{
  int i;
  
  len=len_;
  history=history_;
  sum.resize(len,0);
  sum_tot.resize(len,0);
  sum_block.resize(len,0);
  sum2_tot.resize(len,0);
  iBlock=0;
  n=0;
  n_tot=0;
  label=label_;
  
  c=new counter(0,0);
  n_block_jump=1;
}

void measure_vector::add(vector<double> &value,double step)
{
   int i;
   
   for(i=0;i<len;i++)
    {
      add(value[i],i,step);
    }
  increment_index();
}
void measure_vector::add(double value,int i, double step)
{
  sum[i]=sum[i] + value;
  
  last_step=step;
  
}

void measure_vector::increment_value(double value,int i, double step)
{
  sum[i]=sum[i] + value;
  last_step=step;  
}

void measure_vector::increment_index()
{
  
  n=n+1;
}
// increments the index of a certain system
void measure_vector_mult_index::increment_index()
{
  int i;
  for (i=0;i<len;i++)
    {
      ns[i]=ns[i]+1;
    }
  
}
void measure_vector_mult_index::increment_weight(double weight)
{
  int i;
  for (i=0;i<len;i++)
    {
      increment_index(i,weight);
    }
  
}
// record a measure taken in a certain position
void measure_vector_mult_index::record_block(double time_step)
{
  ofstream f;
  int i,nb_tmp;
  f.open((label + string(".record.dat")).c_str(),std::ios_base::app);
  nb_tmp=0;
  for(i=0;i<len;i++)
    {
      if (ns_block[i]!=0)
	{
      f << i <<" "<< sum_block[i]/ns_block[i]<<" "<< time_step << endl;
      nb_tmp++;
	}
    }
  //cout <<"nb: "<< nb_tmp<<endl;
  f.close();
  
}

void measure_vector_mult_index::increment_index(int i)
{  
      ns[i]=ns[i]+1;
      //cout << "inc"<<i;
}

void measure_vector::clear()
{
  int i;
  n=0;
  for(i=0;i<len;i++)
    {
      sum[i]=0;
    }
  
}

// compute a multiple index
void measure_vector_mult_index::clear()
{
  int i;
  n=0;
  //cout<<"clear"<<endl;
  for(i=0;i<len;i++)
    {
      sum[i]=0;
      ns[i]=0;
    }
  
}

vector<double> measure_vector::mean()
{
  int i;
  vector<double> value;
  value.resize(len);
  for(i=0;i<len;i++)
    {
      assert(n>0);
      value[i]=sum[i]/n;
    }
  return value;
}

vector<double> measure_vector_mult_index::mean()
{
  int i;
  vector<double> value;
  value.resize(len);
  for(i=0;i<len;i++)
    {
      assert(ns[i]>0);
      value[i]=sum[i]/ns[i];
    }
  return value;
}



void measure_vector::out()
{
  ofstream f;
  int i;
  
  if (history)
    {
      f.open((label + string(".dat")).c_str(),std::ios_base::app);
      for(i=0;i<len;i++)
	{
	  if (n > 0)
	    {
	      f <<i<< " "<< sum[i]/n<<" "<<last_step<<endl;
	    }
      
	}
      f.close();
    }
  out_t();
}
void measure_scalar::save()
{
  ofstream f;
  f.open((label + string(".save.dat")).c_str());
  f<<sum_m_t<<" "<<n_t<<endl;
  f.close();
  dec.save(label);
  //out_t();
}

void measure_scalar::load()
{
  ifstream f;
  
  f.open((label + string(".save.dat")).c_str());
  if (f.good())
	 {
  
	   f>>sum_m_t;
	   f>>n_t;
	 }
    
  f.close();
  
  dec.load(label);
  
}

void measure_vector::save()
{
  ofstream f;
  int i;
  f.open((label + string(".save.dat")).c_str());
  f<<len<<endl;
  for(i=0;i<len;i++)
    {
      f<<i<<" "<<sum_tot[i]<<" "<<sum2_tot[i]<<" "<<n_tot<<endl;
      
    }
  dec.save(label);
  f.close();
  
}

void measure_vector_mult_index::save()
{
  ofstream f;
  int i;
  f.open((label + string(".save.dat")).c_str());
  f<<len<<" "<<n_blocks<<endl;
  for(i=0;i<len;i++)
    {
      f<<i<<" "<<sum_tot[i]<<" "<<sum2_tot[i]<<" "<<ns_tot[i]<<" "<<nb[i]<<" "<<sum_of_blocks[i]<<endl;
    }
    
  f.close();
  dec.save(label);
}

void measure_vector_mult_index::load()
{
  ifstream f;
  int i;
  int len_;
  
  f.open((label + string(".save.dat")).c_str());
  if (f.good())
   
    {
       f>>len_;
       f>>n_blocks;
       assert(len==len_);
      for(i=0;i<len;i++)
	{
	  f>>i;
	  f>>sum_tot[i];
	  f>>sum2_tot[i];
	  f>>ns_tot[i];
	  f>>nb[i];
	  f>>sum_of_blocks[i];
	}
    }
  f.close();
  dec.load(label);
}
  
void measure_vector::load()
{
  int len_;
  ifstream f;
  int i;
  
  
  
  f.open((label + string(".save.dat")).c_str());
  if (f.good())
   
    {
       f>>len_;
       assert(len==len_);
      for(i=0;i<len;i++)
	{
	  f>>i;
	  f>>sum_tot[i];
	  f>>sum2_tot[i];
	  f>>n_tot;
	  
	}
    }
  dec.load(label);
  f.close();
  
}

void measure_vector::out_t()
{
  ofstream f;
  int i;
  double error;

  dec.out(label);
  if (n_tot > 0 & n_tot%n_block_jump==0)
    {
      f.open((label + string("_t.dat")).c_str());
      for(i=0;i<len;i++)
	{
	  error=dec.get_error(i);
	  f <<i<< " "<< sum_tot[i]/n_tot<<" "<<error<<" "<<dec.is_converged(i)<<endl;
      
	}
      f.close();
      
    }
  
}

void measure_vector_mult_index::out_t()
{
  ofstream f;
  int i;
  double error;
  
  dec.out(label);
  f.open((label + string("_t.dat")).c_str());
  for(i=0;i<len;i++)
    {
      if (ns_tot[i]>0)
   	{
   	  f <<i<< " "<< sum_tot[i]/ns_tot[i]<<" "<<dec.get_error(i)<<" "<<dec.is_converged(i)<<endl;
   	}
      
    }
  
  f.close();
  
}

// averages a vector with an other vector
/*
void measure_scalar::average_with(measure* ms_)
{
  
  unsigned int i;
  measure_scalar* ms;
  ms=dynamic_cast<measure_scalar*>(ms_);
  if (ms != 0)
    {
      sum_m=sum_m + ms->sum_m;
      sum_m_t=sum_m_t + ms->sum_m_t;
      error_m=error_m + ms->error_m;
      n=n+ms->n;
      n_t=n_t + ms->n_t;
      
      if (last_step < ms->last_step)
	{
	  last_step=ms->last_step;
	}
      if (history)
	{
	  if (recorded_measurements.size()==0)
	    {
	      recorded_measurements.resize(ms->recorded_measurements.size());
	      recorded_times.resize(ms->recorded_measurements.size());
	   
	      for(i=0;i<recorded_measurements.size();i++)
		{
		  recorded_times[i]=ms->recorded_times[i];
		  recorded_measurements[i]=0;
		}
	   
	    }
	  
	  assert(recorded_measurements.size() == ms->recorded_measurements.size());
	  assert(recorded_measurements.size() == recorded_times.size());
	
	  for(i=0;i<recorded_measurements.size();i++)
	    {
	      recorded_measurements[i]=recorded_measurements[i] + ms->recorded_measurements[i];
	      
	      //cout<<recorded_times[i]<<" "<<ms->recorded_times[i]<<endl;
	      
	      assert(recorded_times[i]==ms->recorded_times[i]);
	    }
	  n_recorded=n_recorded+ms->n_recorded;
	}
    }
  else
    {
      cout<<"Downcasting unsuccessfull"<<endl;
      exit(1);
    }
}
*/
/*
void measure_vector::average_with(measure* ms_)
{
  
  int i;
  measure_vector* ms;
  ms=dynamic_cast<measure_vector*>(ms_);
  if (ms != 0)
    {
  assert(len==ms->len);
  assert(len==sum.size());
  last_step=ms->last_step;
  
  for(i=0;i<len;i++)
    {
      sum[i]=sum[i] + ms->sum[i];
      sum_tot[i]=sum_tot[i] + ms->sum_tot[i];
      n=n+ms->n;
      n_tot=n_tot + ms->n_tot;
    }
    }
  else
    {
      cout<<"Downcasting unsuccessfull"<<endl;
      exit(1);
    }
  
}
*/
// average over different measurements
/*
void measures::average_with(measures* m2)
{
  int i;
  for(i=0;i<ms.size();i++)
    {
      ms[i]->average_with(m2->ms[i]);
    }
  
  
}
*/

// overlap_measure::overlap_measure(string label,string wave_o_name,qmc* qmc_obj) : measure_scalar(label,0)
// {
//   wave_o=new total_wavefunction(qmc_obj);
//   wave_o->link_wavefunctions(qmc_obj->main_input,qmc_obj->waves,wave_o_name);
  
// }

// performs the required measurements on the whole system

void measure_scalar::print()
{
  
  cout<<label<<"(current): "<<mean()<<endl;
    cout<<label<<": "<<mean_t()<<"+\\-"<<get_error()<<" ";
  if (is_converged())
    {
      cout<<"Converged.";
    }
  else
    {
      cout<<"Not converged";
    }
  cout<<endl;
}

void measure_vector::print()
{
  int i;
  vector <double> tmp;
  
  tmp=mean();
  cout<<"label: ";
  for(i=0;i<len;i++)
    {
      cout<<tmp[i]<<",";
    }
  cout<<endl;
}


void measure_scalar::increment()
{
  c->increment();
}
void measure_vector::increment()
{
  c->increment();
}
void measure_scalar::set(int jumps,int skip)
{
  c->set(jumps,skip);
}

void measure_vector::set(int jumps,int skip)
{
  c->set(jumps,skip);
}

int measure_vector::check()
{
  return c->check();
}
int measure_scalar::check()
{
  return c->check();
}

void measure_scalar::set_history(bool hist)
{
  history=hist;
}

void measure_vector::set_history(bool hist)
{
  history=hist;
}

void measure_vector::clear_history()
{
  
}

// average the winding number
// template<class qt>
// void center_of_mass_w<qt>::time_difference_average(typename qt::all_particles_t * state,typename qt::wave_t* wave,vector<double> &sum,vector<double> &ns)
// {  
//   int j;
//   this->add(wave->center_of_mass_no_pbc(state,this->set_a));
  
//   for (j=0;j<=this->last;j++)
//     {
//       sum[j]=sum[j] + pow(this->ms[this->last] - this->ms[this->last - j],2);
//       ns[j]=ns[j] + 1;
//     }
  
// }


// performs measurements on the system
//template<class T>
// template<class walker_t,class wave_t>
// void time_difference<T>::make_measurement(walker_t* state,wave_t* wave)
// {
//   int j,i;
//   double w=0;
  
//   //T* wm=static_cast<T*>(state->md[i_wm]);
  
//   //cout <<wm->ms[wm->last]<< "p: " << dmc_obj->current_step<<endl;
//   //cout<<dmc_obj->i_walker<<":"<<wm->last<<" :";
//   //for(i=0;i<wm->bins;i++)
//   //  {
//   //    cout<<wm->ms[i]<<",";
//   //  }
//   //cout<<">>"<<endl;
  
  
//   // adds the recorded value
//   //cout <<wm->ms[wm->last]<< "a: " << dmc_obj->current_step<<endl;
//   //cout << dmc_obj->ws-n<<endl;
//   //cout <<dmc_obj->i_walker<<" "<< dmc_obj->ws->ws[dmc_obj->i_walker]->e<<endl;
//     // if (wm->filled) 
//     //    {
//     //     i=wm->last;
//     //     for(j=len-1;j>=0;j--)
//     // 	{
//     // 	  i=(i+1)%(wm->bins);
//     // 	  sum[j]=sum[j] + pow(wm->ms[wm->last] - wm->ms[i],2);
//     // 	  ns[j]=ns[j]+1;
//     // 	}
//     // 	// increment index
//     //    }
//     // else
//     //   {
//   //i=wm->last;
  
  
  
//   //}
    
// }

// template<class T>
// time_difference<T>::time_difference(int len_,string label_, int i_wm_):
//   measure_vector_mult_index(len_,label_,false)
// {
//   int i=0;
//   i_wm=i_wm_;
//   ns_tot.resize(len);
  
//   for(i=0;i<len;i++)
//     {
//       ns_tot[i]=0;
//     }
  
// }

// makes a dynamic scalar measurement

measure_scalar& measure_scalar::operator()(const measure_scalar & m2)
{
  measure_scalar(m2.label);
  n=m2.n;
  n_t=m2.n_t;
  last_step=m2.last_step;
  sum_m_t=m2.sum_m_t;
  sum_m=m2.sum_m;
  error_m=m2.error_m;
  c=new counter((*c));
  label=m2.label;
  dec=m2.dec;
  history=m2.history;
  n_recorded=m2.n_recorded;
  recorded_measurements=m2.recorded_measurements;
  recorded_times=m2.recorded_times;
  return *this;
}
counter& counter::operator()(counter &c)
{
  status=c.status;
  skip=c.skip;
  termal=c.termal;
  i=c.i;
  jumps=c.jumps;
}
// copy assignement operator
measure_vector& measure_vector::operator()(const measure_vector &m2)
{
  measure_vector(m2.len,m2.label,m2.history);
  len=m2.len;
  c=new counter(*c);
  last_step=m2.last_step;
  history=m2.history;
  n=m2.n;
  vector_n=m2.vector_n;
  n_tot=m2.n_tot;
  sum=m2.sum;
  sum_tot=sum_tot;
  sum2_tot=sum2_tot;
  sum_block=sum_block;
  dec=m2.dec;
}

// template class measures<dmc<D1_t> >;
// template class measures<dmc<spinor1D> >;
// template class measures<vmc<D1_t> >;
// template class measures<vmc<spinor1D> >;
// template class measure_dynamic_scalar<dmc<D1_t> >;
// template class measure_dynamic_scalar<dmc<spinor1D> >;

measure_vector* build_measure_vector(xml_input* xml_m,string label)
{
  int setA,setB,bins,jumps,skip;
  bool history;
  measure_vector* m; // pointer to the space vector measurement
  
  if (xml_m->get_attribute("setA") != NULL)
    {
      setA=xml_m->get_int();
    }
  else
    {
      setA=0;
    }
  if (xml_m->get_attribute("setB") != NULL)
    {
      setB=xml_m->get_int();
    }
  else
    {
      setB=0;
    }
  if (xml_m->get_attribute("history") != NULL)
    {
      history=xml_m->get_bool();
    }
  else
    {
      history="false";
    }
  if (xml_m->get_attribute("bins") != NULL)
    {
      bins=xml_m->get_int();
    }
  else
    {
      bins=100;
    }
  
  if (xml_m->get_attribute("jumps") != NULL)
    {
      jumps=xml_m->get_int();
    }
  else
    {
      jumps=0;
    }
  if (xml_m->get_attribute("skip") != NULL)
    {
      skip=xml_m->get_int();
    }
  else
    {
      skip=0;
    }
  
  m=new measure_vector(bins,label,history);
  m->set_a=setA;
  m->set_b=setB;
  m->set_history(history);
  m->set(jumps,skip);
  return m;
  
}

double getDeltaQ(double l_box,double qMax,int bins)
{
  
  double deltaQ;
  deltaQ=int(qMax/(2*M_PI*bins/l_box))*2*M_PI/l_box;
  
  if (deltaQ==0)
    {
      deltaQ=2*M_PI/l_box;
    }
  
  return deltaQ;
  
}

// builds the Q vector for the static structure factor
vector<double> build_q_vector(int bins,double l_box,double qMax)
{
  vector<double> qs;
  double deltaQ;
  
   int i;
   
   deltaQ=getDeltaQ(l_box,qMax,bins);
   // change the number of bins
   qs.resize(bins);
   
   qs[0]=2*M_PI/l_box;
   
   for(i=1;i<bins;i++)
     {
       qs[i]=qs[0]+i*deltaQ;
     }
   return qs;
   
}

