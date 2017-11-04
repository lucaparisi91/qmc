#include <cstdlib>
#include <cassert>

template<class walkers_t>
void dock::pack(walkers_t* ws)
{
  int i;
  MPI_Request req;
  int mpi_task;
    // pack the data
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_task);
  
  for(i=0;i<received_array.size();i++)
    {
      assert((ws->n - 1 - i) >= 0);
      received_packed_walkers[i]->position=0;
      //cout << "p:"<<mpi_task<<"->"<<received_array[i]<<" "<<ws->ws[ws->n - 1 - i]->e<<endl;
      ws->ws[ws->n - 1 - i]->pack(received_packed_walkers[i]);
    }
  ws->n=ws->n - received_array.size();  
}
template<class walkers_t,class qmc_t>
void dock::unpack(walkers_t *ws,qmc_t* dmc_obj)
{
  int i;
  MPI_Request req;
  int mpi_task;
  typedef typename walkers_t::walker_t walker_t;
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_task);
  
  // creates the directory of new walkers
  for(i=0;i<sent_array.size();i++)
    {
      ws->n=ws->n+1;
      if (ws->n > ws->ws.size())
	{
	  ws->ws.push_back(new walker_t(dmc_obj));
	}
  
    }
  for(i=0;i<sent_array.size();i++)
    {
      sent_packed_walkers[i]->position=0;
      
      ws->ws[ws->n - 1 - i]->unpack(sent_packed_walkers[i]);
      //cout << "u:"<<mpi_task<<"<-"<<sent_array[i]<<" "<<ws->ws[ws->n - 1 - i]->e<<endl;
    }
  
}
