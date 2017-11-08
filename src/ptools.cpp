#include "mpi.h"
#include "traits.h"
#include "dmc.h"
#include "ptools.h"
#include <cmath>
#include "tools.h"
#include "cassert"

#include <cstdlib>
#include <iostream>

using namespace std;
// packs data in the buffer

//packs some information in the system count times
void packed_data::pack(double & in,int count)
{
  MPI_Pack(&in,count,MPI_DOUBLE,buf,size,&position,MPI_COMM_WORLD);
}

void packed_data::pack(vector<double> & in)
{
  pack(in[0],in.size() );
}
//packs some information in the system count times
void packed_data::pack(int& in,int count)
{
  MPI_Pack(&in,count,MPI_INT,buf,size,&position,MPI_COMM_WORLD);
}

void packed_data::pack(bool& in,int count)
{
  MPI_Pack(&in,count,MPI::BOOL,buf,size,&position,MPI_COMM_WORLD);
}



void packed_data::pack(complex<double> &in,int count)
{
  MPI_Pack(&in,count,MPI_DOUBLE_COMPLEX,buf,size,&position,MPI_COMM_WORLD);
}

void packed_data::unpack(complex<double> &in,int count)
{
  MPI_Unpack(buf,size,&position,&in,count,MPI_DOUBLE,MPI_COMM_WORLD);
}

void packed_data::unpack(bool& in,int count)
{
  MPI_Unpack(buf,size,&position,&in,count,MPI::BOOL,MPI_COMM_WORLD);
}

//packs some information in the system count times
void packed_data::unpack(double& in,int count)
{
 MPI_Unpack(buf,size,&position,&in,count,MPI_DOUBLE,MPI_COMM_WORLD);
}

void packed_data::unpack(vector<double> & in)
{
  unpack(in[0],in.size());
}
//packs some information in the system count times
void packed_data::unpack(int& in,int count)
{
  MPI_Unpack(buf,size,&position,&in,count,MPI_INT,MPI_COMM_WORLD);
}

// reserves some space for the buffer
void packed_data::reserve(int n)
{
  int i;
  char* buf_tmp;
  
  buf_tmp= new char[size];
  for(i=0;i<size;i++)
    {
      buf_tmp[i]=buf[i];
    }
  delete [] buf;
  
  buf=new char[n];
  size= n > size ? size : n; 
  for(i=0;i<size;i++)
    {
      buf[i]=buf_tmp[i];
    }
  size=n;
  delete [] buf_tmp;
  
}

void packed_data::reserveReset(int n )
{
  if ( n > size)
    {
      delete [] buf;
      buf=new char[n];
    }
  size=n;
}
// asyncronously send
void packed_data::isend(int dest,int tag,MPI_Request &req)
{
  MPI_Isend(buf,pack_size,MPI_PACKED,dest,tag,MPI_COMM_WORLD,&req);
}
//asyncronously receive
void packed_data::irecv(int src,int tag,MPI_Request &req)
{
  MPI_Irecv(buf,pack_size,MPI_PACKED,src,tag,MPI_COMM_WORLD,&req);
}
void packed_data::send(int dest,int tag)
{
  MPI_Send(buf,pack_size,MPI_PACKED,dest,tag,MPI_COMM_WORLD);
}
//asyncronously receive the walkers
void packed_data::recv(int src,int tag,MPI_Status &stat)
{
  MPI_Recv(buf,pack_size,MPI_PACKED,src,tag,MPI_COMM_WORLD,&stat);
}

packed_data::packed_data()
{
  size=1;
  buf=new char[1];
}

//-----------------------------------------------

void dock::allocate_packs()
{
  int i;
  MPI_Request req;
  
  n_recv_packs=0;
  // allocates the required memory
  
  for(i=0;i<received_array.size();i++)
    {
      n_recv_packs++;
      if (n_recv_packs > received_packed_walkers.size())
	{
	  
	  received_packed_walkers.push_back(new packed_data());
	  received_packed_walkers[n_recv_packs - 1]->reserve(pack_size);
	  received_packed_walkers[n_recv_packs - 1]->pack_size=pack_size;
	  recv_requests.push_back(req);
	  
	}
      
    }
  // set to zero the total number of walkers
  n_send_packs=0;
  for(i=0;i<sent_array.size();i++)
    {
      n_send_packs++;
      if (n_send_packs > sent_packed_walkers.size())
	{
	  sent_packed_walkers.push_back(new packed_data());
	  sent_packed_walkers[n_send_packs - 1]->reserve(1000000);
	  send_requests.push_back(req);
	  sent_packed_walkers[n_send_packs - 1]->pack_size=pack_size;
	  
	}
      
    }
  
}

// send packs over to a certain system
void dock::isend_packs()
{
  int i;
  
  for(i=0;i<received_array.size();i++)
    {
      received_packed_walkers[i]->isend(received_array[i],100,recv_requests[i]);
      
    }
  
}

// receive packs
void dock::irecv_packs()
{
  int i;
  
  for(i=0;i<sent_array.size();i++)
    {
      sent_packed_walkers[i]->irecv(sent_array[i],100,send_requests[i]);
      
    }
  
}
// send and receive packs at the same time
void dock::send_recv_packs()
{
  int i,j1,j2;
  MPI_Status stat;
  j1=0;
  j2=0;
  
  //for (i=0;i<all_sent.size();i++)
    //{
      // cout << all_sent[i] << " "<<all_recv[i]<<";";
      //}
  //cout <<":"<<mpi_task<<":"<<sent_packed_walkers.size()<<" "<< received_packed_walkers.size()<< " "<< sent_array.size()<< " " <<received_array.size()<< endl;
  for(i=0;i<all_sent.size();i++)
    {
      
      if (all_sent[i]==mpi_task)
	{
	  
	  
	  received_packed_walkers[j1]->send(all_recv[i],100);
	  j1++;
	}
	if (all_recv[i]==mpi_task)
	{
	  
	  sent_packed_walkers[j2]->recv(all_sent[i],100,stat);
	  j2++;
	}
      
    }
  
    
  
  //assert(sent_packed_walkers.size()==j1==j2);
  
}

void dock::wait_packs()
{
  int i;
  MPI_Status stat;
  
  for(i=0;i<send_requests.size();i++)
    {
      MPI_Wait(&send_requests[i],&stat);
    }
  for(i=0;i<recv_requests.size();i++)
    {
      MPI_Wait(&recv_requests[i],&stat);
    }
//requires first determining sendings and receivings
}

int dock::get_walkers()
{
  int i;
  int nw;
  nw=0;
  for(i=0;i<core_populations.size();i++)
    {
      nw+=core_populations[i];
    }
  return nw;
}
void dock::send_receive_determine()
{
  int n_walkers=0;
  double diff=0;
  int i_send,i_rec,i;
  int tag;
  vector<int> diffs(mpi_tasks);
  
  sent_array.resize(0);
  all_sent.resize(0);
  all_recv.resize(0);
  received_array.resize(0);
  tags_sent.resize(0);
  tags_rec.resize(0);
  n_walkers=0;
  tag=0;
  for(i=0;i<mpi_tasks;i++)
    {
      n_walkers=n_walkers + core_populations[i];    
    }
  diff=0;
  
  //cout << "n_walkers: "<<n_walkers<< " ";
  
  for(i=0;i<mpi_tasks;i++)
    {
      diffs[i]= core_populations[i] - optimal_n(n_walkers,mpi_tasks,i);
      //cout <<core_populations[i] << " "<< optimal_n(n_walkers,mpi_tasks,i)<<"#"<<mpi_task<<"#";
      diff=diff + abs(diffs[i]);
    }
  //cout<<endl;
  //cout << "diff: "<< diff<<" ;";
  // as long as the sum of differencies is different from zero
  
  while(diff != 0)
    {
      // index of the vector to send
      i_send=vec_max_int(diffs);
      // index of the vector to receive
      i_rec=vec_min_int(diffs);

      assert(diffs[i_send] > 0);
      assert(diffs[i_rec] < 0);
      
      diffs[i_send]=diffs[i_send] - 1;
      diffs[i_rec]=diffs[i_rec] + 1;
         
      assert(i_send != i_rec);
      tag=tag + 1;
      if (i_send == mpi_task)
	{
	  received_array.push_back(i_rec);
	  tags_rec.push_back(tag);
	}
      if (i_rec == mpi_task)
	{
	  sent_array.push_back(i_send);
	  tags_sent.push_back(tag);
	}
      all_sent.push_back(i_send);
      all_recv.push_back(i_rec);
      diff=diff - 2;
    }
 }
// creates a docking application
dock::dock()
{
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_task);
  core_populations.resize(mpi_tasks);
  sent_array.resize(0);
  received_array.resize(0);
  tags_sent.resize(0);
  tags_rec.resize(0);
  sent_packed_walkers.resize(0);
  send_requests.resize(0);
  recv_requests.resize(0);
  n_send_packs=0; // number of packs to send
  n_recv_packs=0; // number of packs to receive
  
  populations_send_requests=new MPI_Request[mpi_tasks];
  populations_rec_requests=new MPI_Request[mpi_tasks];
  populations_status=new MPI_Status[mpi_tasks];
  energies_send_requests=new MPI_Request[mpi_tasks];
  energies_recv_requests=new MPI_Request[mpi_tasks];
  energies_recv_status=new MPI_Status[mpi_tasks];
  energies_send_status=new MPI_Status[mpi_tasks];
  
  energies.resize(mpi_tasks);
}

void dock::exchange_populations(int n_walkers)
{
  int i;
  //cout<<"exch:"<< mpi_task<< " "<<n_walkers<<endl;
  core_populations[mpi_task]=n_walkers;
  for(i=0;i<mpi_tasks;i++)
    {
      MPI_Isend(&(core_populations[mpi_task]), 1, MPI_INT,i,500, MPI_COMM_WORLD,&populations_send_requests[i]);
      MPI_Irecv(&(core_populations[i]), 1, MPI_INT,i,500, MPI_COMM_WORLD,&populations_rec_requests[i]);
    }
  //MPI_Allgather(&n_walkers,1,MPI_INT,&core_populations.front(),1,MPI_INT,MPI_COMM_WORLD);
  //cout<<mpi_task<< "cores: ";
  
  // for(i=0;i<mpi_tasks;i++)
  //  {
  //    cout <<core_populations[i]<<",";
  //  }
  //cout<<endl;
}

void dock::gather_populations(int n_walkers)
{
  core_populations[mpi_task]=n_walkers;
  MPI_Allgather(&n_walkers,1,MPI_INT,&core_populations.front(),1,MPI_INT,MPI_COMM_WORLD);
  
}
// wait for the walker to be populated
void dock::send_energy(double &e)
{
  
  int i=0;
  for(i=0;i<mpi_tasks;i++)
    {
      if (i!=mpi_task)
	{
      MPI_Isend(&e, 1, MPI_DOUBLE,i,600, MPI_COMM_WORLD,&energies_send_requests[i]);
	}
      
    }
}
// receives the energy e from other systems
void dock::recv_energy()
{
  int i=0;
  for(i=0;i<mpi_tasks;i++)
    {
      if (i!=mpi_task)
	{
	  MPI_Irecv(&energies[i], 1, MPI_DOUBLE,i,600, MPI_COMM_WORLD,&energies_recv_requests[i]);
	}
    }
  
  
}


void dock::wait_energy()
{
  int i;
  for(i=0;i<mpi_tasks;i++)
    {
      if (i!= mpi_task)
	{
	  MPI_Wait(&energies_send_requests[i],&energies_send_status[i]);
	  MPI_Wait(&energies_recv_requests[i],&energies_recv_status[i]);
	}
    }
  
}
double dock::energy()
{
  int i=0;
  double e=0;
  
  for(i=0;i<mpi_tasks;i++)
    {
      e+=energies[i];
    }
  return e/mpi_tasks;
}

void dock::wait_populations()
{
  int i;
  for(i=0;i<mpi_tasks;i++)
    {
      if (i!= mpi_task)
	{
	  MPI_Wait(&energies_send_requests[i],&energies_send_status[i]);
	  MPI_Wait(&energies_recv_requests[i],&energies_recv_status[i]);
	}
    }

  MPI_Waitall(mpi_tasks,populations_send_requests,populations_status);
  MPI_Waitall(mpi_tasks,populations_rec_requests,populations_status);
}

int optimal_n(int n_walkers,int tasks, int task)
{
   int n_optimal;
   int n_remainder;
   
   n_optimal=0;
   n_remainder=n_walkers%tasks;
   if (task < n_remainder)
     {
       n_optimal=1;
     }
   n_optimal=n_optimal + n_walkers/tasks;
   return n_optimal;
}
// returns the number of total walkers
int dock::total_walkers()
{
  int i;
  int total_walkers;
  total_walkers=0;
  for(i=0;i<mpi_tasks;i++)
    {
      total_walkers=total_walkers + core_populations[i];
    }
  return total_walkers;
}

namespace pTools
{
template<>
int get_pack_size<double>(const double & x)
{
  int size;
  MPI_Pack_size(1,MPI_DOUBLE,MPI_COMM_WORLD,&size);
  return size;
};

template<>
int get_pack_size<int>(const int & x)
{
  int size;
  MPI_Pack_size(1,MPI_INT,MPI_COMM_WORLD,&size);
  return size;
};
  //returns the size of the system
template<>
int get_pack_size<bool>(const bool & x)
{
  int size;
  MPI_Pack_size(1,MPI::BOOL,MPI_COMM_WORLD,&size);
  return size;
}

template<>
// returns the size of the algorithm
int get_pack_size< vector<double> >(const vector<double> & x)
{
  int size;
  MPI_Pack_size(x.size(),MPI_DOUBLE,MPI_COMM_WORLD,&size);
  return size;
};

template<>
int get_pack_size<vector<int> >(const vector<int> & x)
{
  int size;
  MPI_Pack_size(x.size(),MPI_INT,MPI_COMM_WORLD,&size);
  return size;
};
  // get pack size for the system
  
template<>
int get_pack_size< vector<complex<double> > >(const vector<complex<double> > & x)
{
  int size;
  MPI_Pack_size(x.size(),MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,&size);
  return size;
};

template<>
int get_pack_size< vector<vector<complex<double> > > >(const vector< vector<complex<double> > > & x)
{
  int size;
  MPI_Pack_size(get_counts(x),MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,&size);
  return size;
}
  
}

void dock::print_populations()
{
  int i;
  cout << "Populations"<<endl;
  for(i=0;i<core_populations.size();i++)
    {
      cout << i << " "<<core_populations[i] << endl;
    }
}

void pTools::transferSum(vector<double> &sum,int root)
{
  // saves the sum of tasks in the same input vector if root
  // otherwise reset to zero
  vector<double> sumMPIStore;
  sumMPIStore.resize(sum.size());
  
  MPI_Reduce(&sum.front(),&sumMPIStore.front(), sum.size(), MPI_DOUBLE, MPI_SUM, root,MPI_COMM_WORLD);
  //copy back into input vector
  
  for(int i=0;i<sum.size();i++)
    {
      sum[i]=sumMPIStore[i];
    } 
  
}

void pTools::transferSum(int &n,int root)
{
  // saves the sum of tasks in the same input variable if root otherwise reset to zero
  int MPIStore;
  
  MPIStore=0;
  
  MPI_Reduce(&n,&MPIStore, 1, MPI_INT, MPI_SUM, root,MPI_COMM_WORLD);
  
  n=MPIStore;
}
