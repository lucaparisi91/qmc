#include "random.h"
#include "cassert"
#include "input.h"
#include "system.h"
#include "tools.h"
#include "measures.h"
#include "geometry.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include "wavefunction.h"
#include "omp.h"
#include "xml-input.h"
#include "mpi.h"
#include "ptools.h"
#include "timer.h"
#include "dmc.h"
#include <string.h>
#include <string>

using namespace std;

// returns the optimal number of walkers for the system
// performs some branching at a certain interger position


// void walker::send(gatherer *ws_global)
// {
//   if (ws_global->n >= static_cast<signed int>(ws_global->ws.size()))
//     {
//       ws_global->ws.push_back(new walker(dmc_obj));
//     }
//   ws_global->n=ws_global->n+1;
//   // copy the walker contents
//   *(ws_global->ws[ws_global->n-1])=*this;
// }

// void walker::send_mpi(int dest,int tag,vector<MPI_Request>& send_requests )
// {
//   int mpi_task;
//   vector<MPI_Request> state_requests;
//   send_requests.resize(3);
//   MPI_Isend(&wavefunction_value, 1, MPI_DOUBLE, dest, tag , MPI_COMM_WORLD,&send_requests[0]);
//   MPI_Isend(&descendants, 1, MPI_INT, dest, tag + 1 , MPI_COMM_WORLD,&send_requests[1]);
//   MPI_Isend(&e, 1, MPI_DOUBLE, dest, tag + 2 , MPI_COMM_WORLD,&send_requests[2]);
//   state->send_mpi(dest,tag + 3,state_requests);
//   // add to the the particles request
//   send_requests.insert(send_requests.end(),state_requests.begin(),state_requests.end());
  
  
  
//   //MPI_Send(&drift_force.front(), drift_force.size(), MPI_DOUBLE, dest, tag + 4, MPI_COMM_WORLD);
//   //MPI_Send(&work2.front(), work2.size(), MPI_DOUBLE, dest, tag + 5, MPI_COMM_WORLD);
//   MPI_Comm_rank(MPI_COMM_WORLD,&mpi_task);
//   //cout<< "send:"<< mpi_task << " " << dest<< " " << state->position_no_pbc[0]<<endl;
// }
// // receives a mpi message with tag 'tag' from source 'dest'

// void walker::receive_mpi(int dest,int tag,vector<MPI_Request> & recv_requests )
// {
//   int mpi_task;
//   MPI_Status stat;
//   vector<MPI_Request> state_requests;
//   recv_requests.resize(3);
//   MPI_Irecv(&wavefunction_value, 1, MPI_DOUBLE, dest, tag , MPI_COMM_WORLD,&recv_requests[0]);
//   MPI_Irecv(&descendants, 1, MPI_INT, dest, tag + 1 , MPI_COMM_WORLD,&recv_requests[1]);
//   MPI_Irecv(&e, 1, MPI_DOUBLE, dest, tag + 2 , MPI_COMM_WORLD,&recv_requests[2]);
//   state->receive_mpi(dest,tag + 3,state_requests);
//   recv_requests.insert(recv_requests.end(),state_requests.begin(),state_requests.end());
  
//   //MPI_Recv(&drift_force.front(), drift_force.size(), MPI_DOUBLE, dest, tag + 4, MPI_COMM_WORLD,&stat);
//   //MPI_Recv(&work2.front(), work2.size(), MPI_DOUBLE, dest, tag + 5, MPI_COMM_WORLD,&stat);
  
//   //MPI_Comm_rank(MPI_COMM_WORLD,&mpi_task);
//   //cout<< "rec:"<< mpi_task << " " << dest<< " " << state->position_no_pbc[0]<<endl;
  
// }
// // receives something from a gatherer
// void walker::receive(gatherer* g)
// {
//   // makes sure the number of walkers is less than a certain amount
//   assert(g->n<= static_cast<int>(g->ws.size()));
//   assert(g->n>0);
//   //copy the last walkers in stack on the local walkers
//   *this=*(g->ws[g->n -1]);
//   g->n=g->n-1;
// }
// // send the excess walkers

// void walkers::send(gatherer *g)
// {
//   int n_excess;
//   int i;
  
//   //cout << g->n_walkers<<endl;
  
//   n_excess= n - optimal_n(g->n_walkers,omp_get_num_threads(),omp_get_thread_num());
  
//   if (n_excess>0)
//     {
//       for(i=0;i<n_excess;i++)
// 	{
// 	  ws[n - 1 - i]->send(g);
// 	}
//       n=n-n_excess;
//     }
// }

// void walkers::receive(gatherer *g)
// {
//   int n_excess;
  
//   int i;
  
//   n_excess=n - optimal_n(g->n_walkers,omp_get_num_threads(),omp_get_thread_num());
//   if (n_excess<0)
//     {
//       n_excess=-n_excess;
//       for(i=0;i<n_excess;i++)
// 	{
// 	  n=n+1;
// 	  if (n > static_cast<signed int>(ws.size()))
// 	    {
// 	      ws.push_back(new walker(dmc_obj));
// 	    }
// 	  ws[n - 1]->receive(g);
// 	}
      
//     }
// }

// gatherer::gatherer()
// {
//   e_t=0;
//   n=0;
//   n_walkers=0;
//   saving=false;
//   ws.resize(0);
// }


bool metropolis(double log_ratio,random1* rand)
{ 
  double random_number;
  if (log_ratio >=0)
    {
      return true;
    }
  else
    {
      random_number=rand->uniform();
      assert(random_number>0);
      assert(random_number<=1);
      if ( log_ratio > log(random_number))
	{
	  return true;
	}
      else
	{
	  return false;
	}
    }
}

// determines which copies to send/recieve over the network

// sends the excess walkers to mpi

// void dmc::send_mpi()
// {
//   unsigned int i;
//   vector <int> tasks_to_send;
//   vector <int> n_copies;
//   vector<MPI_Request>walker_requests;
//   //cout<< "step: "<< current_step <<  " "<< mpi_task<<endl;
   
//   /*
//   for(i=0;i<received_array.size();i++)
//     {
      
//       // see if i already have sent tasks
//       found=false
//       for(j=0;j<tasks_to_send.size();j++)
// 	{
// 	  if (tasks_to_send[j] ==received_array[i])
// 	    {
// 	      n_copies[j]=n_copies[j] + 1;
// 	      found=true;
// 	    }
// 	  // if the walker was not found on the system
// 	  if (! found)
// 	    {
// 	      n_copies.push_back(1);
// 	      tasks_to_send.push_back(received_array[i]);
// 	    }
// 	}
//      }
//   j=0;
//   i=0;
//   k=0
//     while( k < received_array.size())
//     {
//       mult=ws->ws[ws->n - j]->n_descendants;
//       if (mult <= n_copies[i])
// 	{
// 	  k=k + mult;
// 	  n_copies[i]=n_copies[i] - mult;
// 	  // send the walker over mpi
// 	}
//       else
// 	{
// 	  n_copies[i]=0;
// 	  k=k + n_copies[i];
// 	  ws->ws[ws->n - j]->n_descendants=n_copies[i];
// 	  // send the walker over mpi
	  
// 	}
//     }
  
// }

//   */
  
  
//   for (i=0;i<received_array.size();i++)
//     {
//       ws->ws[ws->n - 1 - i]->send_mpi(received_array[i],tags_rec[i],walker_requests);
//       walkers_send_requests.insert(walkers_send_requests.end(),walker_requests.begin(),walker_requests.end());
      
//     }
//   /*
//   if (received_array.size() > 0)
//     {
//       cout<<"s:"<<mpi_task<<":";
//       for(i=0;i< received_array.size();i++)
// 	{
// 	  cout<<received_array[i];
// 	}
//       cout<<":step:"<<current_step<<endl;
//     }
//   */
//   // update the population of walkers
//   ws->n=ws->n - received_array.size();
  
//     }

// void dmc::receive_mpi()
// {
//   unsigned int i;
//   vector <int> tasks_to_send;
//   vector <int> n_copies;
//   walker* tmp;
//   vector<MPI_Request> walker_requests; // requests issued by a walker
//   //cout<< "step: "<< current_step <<  " "<< mpi_task<<endl;

//   // add empty walkers equal to the number of walkers sent
  
  
//   for(i=0;i<n_recv;i++)
//     {
//       // adds a single walker
//       ws->n=ws->n+1;
      
//       if (ws->n > static_cast<signed int>(ws->ws.size()))
// 	{
// 	  ws->ws.push_back( new walker(this));
// 	}
//     }
//   /*
//   if (sent_array.size() > 0)
//     {
//       cout<<"r:"<<mpi_task<<":";
//       for(i=0;i< sent_array.size();i++)
// 	{
// 	  cout<<sent_array[i];
// 	}
//       cout<<":step:"<<current_step<<endl;
//     }
//   */
//   // loop over a number of cpus less than the total size
//   // must be less than the size of the array
  
//   for (i=0;i<n_recv;i++)
//     {
//       tmp=ws->ws[ws->n - 1 - i];
//       ws->ws[ws->n - 1 - i]=ws_recv[i];
//       ws_recv[i]=tmp;
      
//     }  
  
// }


// void dmc::receive_mpi_tmp()
// {
  
//   unsigned int i;
//   vector <int> tasks_to_send;
//   vector <int> n_copies;
//   vector<MPI_Request> walker_requests; // requests issued by a walker
//   //cout<< "step: "<< current_step <<  " "<< mpi_task<<endl;

//   // add empty walkers equal to the number of walkers sent
  
//   n_recv=0;
  
//   for(i=0;i<sent_array.size();i++)
//     {
//       // adds a single walker
//       n_recv=n_recv + 1;
      
//       if (n_recv  > static_cast<signed int>(ws_recv.size()))
// 	{
// 	  ws_recv.push_back( new walker(this));
// 	}
//     }

  
//   for (i=0;i<sent_array.size();i++)
//     {
//       ws_recv[i]->receive_mpi(sent_array[i],tags_sent[i],walker_requests);
//       walkers_recv_requests.insert(walkers_recv_requests.end(),walker_requests.begin(),walker_requests.end());
//     }
  
//   //ws->n=ws->n - sent_array.size();
  
// }

// load the walkers from where i saved them

// compute the logarithm of the number of walkers

