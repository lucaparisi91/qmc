#ifndef PTOOLS_H
#define PTOOLS_H

#include "mpi.h"
#include "vector"
#include "cstring"
#include <complex>
using namespace std;
template<class T> class walker;
template<class T> class walkers;
// containes and build the message to be sent to the processor
class packed_data
{
 public:
  int pack_size;
  int size; // the size of the buffer
  char* buf;// the buffer containing all the variables
  int position;
  packed_data();
  
  void pack(int& in,int count);//packs data in the buffer
  void pack(vector<double> &in);
  void unpack(vector<double> &in);
  void pack(int& in){pack(in,1);}
  void unpack(int& in,int count);//packs data in the buffer
  void unpack(int& in){unpack(in,1);}
  void pack(double& in,int count);
  void pack(double& in){pack(in,1);}
  void pack(bool& in,int count);
  void pack(bool& in){pack(in,1);}
  void unpack(bool& in,int count);
  void unpack(bool& in){unpack(in,1);}
  void unpack(double& in,int count);
  void unpack(double& in){unpack(in,1);}
  void pack(complex<double> &in,int count);
  void pack(complex<double>& in){pack(in,1);}
  void unpack(complex<double> &in,int count);
  void unpack(complex<double>& in){unpack(in,1);}
  
  template<class T>
  void pack(vector<T> & vec)
  {
    
    pack( vec[0],vec.size());
  }
  
  template<class T>
  void unpack(vector<T> & vec)
  {
    
    unpack( vec[0],vec.size());
  }

  template<class T>
  void pack(vector<vector<T> > & data)
  {
    for(int i=0;i<data.size();i++)
      {
	pack(data[i]);
      }
  }
  
  template<class T>
  void unpack(vector<vector<T> > & data)
  {
    for(int i=0;i<data.size();i++)
      {
	unpack(data[i]);
      }
  }
  
  void reserve(int n); // allocate n slots in the buffer
  void reserveReset(int n); // allocate n slots in the buffer
  void isend(int dest,int tag,MPI_Request &req);
  void irecv(int src,int tag,MPI_Request &req);
  void send(int dest,int tag);
  void recv(int dest,int tag,MPI_Status &stat);
};
// a dock class

class dock
{
 public:
  int mpi_task;
  int mpi_tasks;
  int pack_size;
  dock();
  MPI_Request* populations_send_requests;
  MPI_Request* populations_rec_requests;
  MPI_Request* energies_send_requests;
  MPI_Request* energies_recv_requests;
  MPI_Status* populations_status;
  MPI_Status* energies_send_status;
  MPI_Status* energies_recv_status;
  vector<double> energies;
  
  vector<int> core_populations;
  vector<int> sent_array;
  vector<int> received_array;
  vector<int> tags_sent;
  vector<int> tags_rec;
  vector<int> all_sent;
  vector<int> all_recv;
  vector<packed_data*> sent_packed_walkers;//packs to send
  vector<MPI_Request> send_requests;
  vector<MPI_Request> recv_requests;
  unsigned int n_send_packs;
  unsigned int n_recv_packs;
  
  vector<packed_data*> received_packed_walkers;//packs to received_array
  template<class walker_t> int get_pack_size(walker_t * ws);
  void allocate_packs();
  void send_packs();
  void recv_packs();
  void isend_packs();
  void irecv_packs();
  // print the populations of the walkers
  void print_populations();
  void send_recv_packs();
  double energy();
  void wait_energy();
  void send_energy(double &e);
  void recv_energy();
  template<class walkers_t> void pack(walkers_t* ws);
  template<class walkers_t,class qmc_t> void unpack(walkers_t* ws,qmc_t* dmc_obj);
  void isend();// send packages asyncrhonously
  void irecv();//send packages in an asynchronous manner
  void send_receive_determine();
  void exchange_populations(int n_walkers);
  void wait_populations();
  void wait_packs();
  int get_walkers();
  int total_walkers();
  void exchange_sizes();
  void gather_populations(int n_walkers);
  
  
  template<class walkers_t> void send_walkers_async(walkers_t* ws)
  {
    gather_populations(ws->n);
    
    send_receive_determine(); // determine which packets to send
    allocate_packs();
    pack(ws);
    isend_packs(); 
    irecv_packs();
    
  };
  template<class walkers_t,class qmc_t> void receive_walkers_async(walkers_t* ws,qmc_t* dmc_obj)
  {
    wait_packs();
    unpack(ws,dmc_obj);
  };
 
};
// returns the optimal number of walkers
int optimal_n(int n_walkers,int tasks, int task);


namespace pTools
{
  template<class T> int get_pack_size(const T & x);
  void transferSum(vector<double> & vec,int root);
  void transferSum(int &n,int root);
  
  void broadcast(vector<double> &vec,int root);
  void broadcast(vector<vector<double> > &vec,int root);
  void broadcast(int &b,int root);
  
}
#include "ptools.hpp"

#endif
