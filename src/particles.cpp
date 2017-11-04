#include "system.h"
#include "input.h"
#include <iostream>
#include <sstream>
#include <string>
#include "qmc.h"
#include "geometry.h"
#include "mpi.h"
#include "random.h"
#include "dmc.h"
#include "ptools.h"
#include "xml-input.h"
#include "tools.h"
#include <cstdlib>
#include <cassert>

using namespace std;
// uses a different kind of namespace
template<class T>
particles<T>::particles(int n_)
{
  resize(n_);
}

template<class T>
particles<T>::particles()
{
  
}


template<class T>
void particles<T>::resize(int n_)
{
  n=n_;
  position.resize(n);
  position_no_pbc.resize(n);
  work.resize(n);
  drift_force_derivative.resize(n);
}

template<class T>
particles<T>& particles<T>::operator=(particles<T> &p)
{
  // if (n != p.n)
  //  {
  //    n=p.n;
  //    p.position.resize(n);
  //  }
  int i;
  
  for(i=0;i<n;i++)
    {
      position[i]=p.position[i];
      position_no_pbc[i]=p.position_no_pbc[i];
      
      drift_force_derivative[i]=p.drift_force_derivative[i];
    }
  
  return *this;
}
template<class T>
int particles<T>::get_pack_size()
{
  int size;
  size=pTools::get_pack_size(n) + pTools::get_pack_size(position) + pTools::get_pack_size(position_no_pbc);
  
  return size;
}
template<class T>
// packs all the particles in a common buffer
void particles<T>::pack(packed_data* packed_particles)
{
  packed_particles->pack(n,1);
  packed_particles->pack(position.front(),position.size());
  packed_particles->pack(position_no_pbc.front(),position_no_pbc.size());
  
}
// unpack particles sets
template<class T>
void particles<T>::unpack(packed_data* packed_particles)
{
  packed_particles->unpack(n,1);
  packed_particles->unpack(position.front(),position.size());
  packed_particles->unpack(position_no_pbc.front(),position_no_pbc.size());
  
}

// packs a set of particles in a certain container

template<class all_particles_t>
int all_particles<all_particles_t>::get_pack_size()
{
  int i;
  int size;
  // pack with particles and get the pack size
  size=0;
  for(i=0;i<particle_sets.size();i++)
    {
      size+=particle_sets[i]->get_pack_size();
    }
  return size;
}
template<class all_particles_t>
void all_particles<all_particles_t>::pack(packed_data* packed_particles)
{
  int i;
  // packed of particles
  for(i=0;i<particle_sets.size();i++)
    {
      particle_sets[i]->pack(packed_particles);
    }  
}
// unpack some particles somewhere
template<class all_particles_t>
void all_particles<all_particles_t>::unpack(packed_data* packed_particles)
{
  int i;
  // packed of particles
  for(i=0;i<particle_sets.size();i++)
    {
      particle_sets[i]->unpack(packed_particles);
    }  
}

template<class all_particles_t>
int all_particles<all_particles_t>::getN()
{
  return particle_sets.size();
}

template<class all_particles_t>
void all_particles<all_particles_t>::getNs(vector<int> & ns)
{
  ns.resize(particle_sets.size());
  for(int i=0;i<particle_sets.size();i++)
    {
      ns[i]=particle_sets[i]->getN();
    }
  
}

// // send the particle object over mpi
// void particles::send_mpi(int dest,int tag, vector<MPI_Request> &send_requests )
// {
  
//   send_requests.resize(4);
//   MPI_Isend(&n,1,MPI_INT,dest,tag + 2,MPI_COMM_WORLD,&send_requests[2]);
  
//   MPI_Isend(&position.front(), position.size(), MPI_DOUBLE, dest, tag, MPI_COMM_WORLD,&send_requests[0]);
//   MPI_Isend(&position_no_pbc.front(), position_no_pbc.size(), MPI_DOUBLE, dest, tag + 1, MPI_COMM_WORLD,&send_requests[1]);
//   MPI_Isend(&drift_force.front(), drift_force.size(), MPI_DOUBLE, dest, tag + 3, MPI_COMM_WORLD,&send_requests[3]);
  
// }

// // receives the particle objects over MPI trough a certain tag
// void particles::receive_mpi(int src,int tag,vector<MPI_Request> &recv_requests )
// {
//   MPI_Status stat;
//   recv_requests.resize(4);
//   MPI_Irecv(&n,1,MPI_INT,src,tag + 2,MPI_COMM_WORLD,&recv_requests[2]);
//   MPI_Irecv(&position.front(), position.size(), MPI_DOUBLE, src, tag, MPI_COMM_WORLD,&recv_requests[0]);
//   MPI_Irecv(&position_no_pbc.front(), position_no_pbc.size(), MPI_DOUBLE, src, tag + 1, MPI_COMM_WORLD,&recv_requests[1]);
//   MPI_Irecv(&drift_force.front(), drift_force.size(), MPI_DOUBLE, src, tag + 3, MPI_COMM_WORLD,&recv_requests[3]);
// }

// updates the position of the particles
// prints a certain number
/*
void all_particles::send_mpi(int dest,int tag)
{
  int i;
  
  for(i=0;i<particle_sets.size();i++)
    {
      particle_sets[i]->send_mpi(dest,tag);
    }
}

void all_particles::receive_mpi(int src,int tag)
{
  int i;

  for(i=0;i<particle_sets.size();i++)
    {
      particle_sets[i]->receive_mpi(src,tag);
    }
}
*/

template<class particles_t>
all_particles<particles_t>& all_particles<particles_t>:: operator=(all_particles<particles_t> & p)
{
  int i;
  particle_sets.resize(p.particle_sets.size());
  for(i=0;i<particle_sets.size();i++)
    {
      (*(particle_sets[i]))=(*(p.particle_sets[i]));
      //if (i==1)
      //	{
      //	  cout << p.particle_sets[i]->position_no_pbc[0] - particle_sets[i]->position_no_pbc[0]<<endl;
      //	}
    }
}
template<class particles_t>
all_particles<particles_t>::all_particles()
{
  particle_sets.resize(0);
  
}
// gaussian system of particles
template<class particles_t>
void all_particles<particles_t>::gaussian(random1* rand_obj)
{
  int i;
  for(i=0;i<particle_sets.size();i++)
    {
      particle_sets[i]->gaussian(rand_obj);
    }
}
// set a unfirom distribution over [-pos/2,pos/2]
template <class T>
void particles<T>::set_uniform(double pos,random1* randg)
{
  int i;
  double step;
  step=pos/position_no_pbc.size();
  for(i=0;i<position_no_pbc.size();i++)
    {
      position_no_pbc[i]=-pos/2 + i*step + 0.1*randg->uniform();
      
    }
  
}

template<class T>
void particles<T>::set_gaussian(random1* randg,double sigma,double position)
{
  int i;
  randg->gaussian(position_no_pbc);
  for(i=0;i<position_no_pbc.size();i++)
    {
      position_no_pbc[i]=position_no_pbc[i]*sqrt(sigma)+position;
    }
  
}


template<class particles_t>
void all_particles<particles_t>::set_gaussian(random1* randg,double alpha,double position)
{
  int i;
  for(i=0;i<particle_sets.size();i++)
    {
      particle_sets[i]->set_gaussian(randg,alpha,position);
    }
}

template<class particles_t>
void all_particles<particles_t>::set_uniform(double pos,random1* randg)
{
  int i;
  for(i=0;i<particle_sets.size();i++)
    {
      particle_sets[i]->set_uniform(pos,randg);
    }
}
template<class particles_t>
void all_particles<particles_t>::set_all_positions(double pois)
{
  int i;
  for(i=0;i<particle_sets.size();i++)
    {
      particle_sets[i]->set_all_positions(pois);
    }
}
// returns some gaussian values
template<class T>
void particles<T>::gaussian(random1* rand_obj)
{
  int i;
  rand_obj->gaussian(position_no_pbc);
}

// set the position of all the walkers in the syste
template<class T>
void particles<T>::set_all_positions(double pois)
{
  int i;
  for(i=0;i<position_no_pbc.size();i++)
    {
      position_no_pbc[i]=pois;
    }
}
template<class particles_t>
void all_particles<particles_t>::save(xml_input* xml_save_particles)
{
  int i;
  xmlNodePtr cur;
  string str;
  cur=xml_save_particles->cur;
  // adds a child to the particle root node
  
  for(i=0;i<particle_sets.size();i++)
    {
      xml_save_particles->cur=cur;
      xml_save_particles->add_child("particles","");
      particle_sets[i]->save(xml_save_particles);
    }
  
  xml_save_particles->cur=cur;
    

}
// load particles from an input file
template<class particles_t>
void all_particles<particles_t>::load(xml_input* xml_save_particles)
{
  int i;
  xmlNodePtr cur;
  
  cur=xml_save_particles->cur;
  // adds a child to the particle root node
  xml_save_particles->get_child("particles");
  //cout<<"name: "<< xml_save_particles->get_name()<<endl;
  
  
  
  for(i=0;i<particle_sets.size() ;i++)
    {    
      particle_sets[i]->load(xml_save_particles);
      xml_save_particles->get_next("particles");
    }
  
  xml_save_particles->cur=cur;
  
}

// loads the particles object from file
template<class T>
void particles<T>::load(xml_input* xml_save_particles)
{
  int i;
  xmlNodePtr cur;
   // string stream to save the particles
  string content;
  
  cur=xml_save_particles->cur;
  //cout << xml_save_particles->get_name()<<endl;
  
  assert(n==xml_save_particles->get_attribute("n")->get_int());
  content=xml_save_particles->get_value()->get_string();
  
  stringstream strs(content);
  for(i=0;i<n;i++)
    {
      strs >> position_no_pbc[i];
    }
  
  xml_save_particles->cur=cur;
}

void spinors::load(xml_input* xml_save_particles)
{
  int i;
  xmlNodePtr cur;
   // string stream to save the particles
  string content;
  
  cur=xml_save_particles->cur;
  //cout << xml_save_particles->get_name()<<endl;
  
  assert(n==xml_save_particles->get_attribute("n")->get_int());
  content=xml_save_particles->get_value()->get_string();
  
  stringstream strs(content);
  for(i=0;i<n;i++)
    {
      strs >> position_no_pbc[i];
      strs >> spinComp[i][0];
      strs >> spinComp[i][1];
    }
  
  xml_save_particles->cur=cur;
  
}
// saves the particles to a certain system
template<class T>
void particles<T>::save(xml_input* xml_save_particles)
{
  int i;
  xmlNodePtr cur;
  stringstream strs;
  string str;
  cur=xml_save_particles->cur;
  strs<<endl;
  for(i=0;i<n;i++)
    {
      strs << position_no_pbc[i] << endl;
      
    }
  
  xml_save_particles->set(strs.str());
  xml_save_particles->set("n",int_to_string(n));
}
// saves the particles to a certain system
void spinors::save(xml_input* xml_save_particles)
{
  int i;
  xmlNodePtr cur;
  stringstream strs;
  string str;
  cur=xml_save_particles->cur;
  strs<<endl;
  for(i=0;i<n;i++)
    {
      strs << position_no_pbc[i] << endl;
      strs << spinComp[i][0] << endl;
      strs << spinComp[i][1] << endl;
      
    }
  
  xml_save_particles->set(strs.str());
  xml_save_particles->set("n",int_to_string(n));
}

// inits the creation of a new set of particles
template<class particles_t>
void all_particles<particles_t>::init(xml_input* xml_main_input)
{
  
  int tmp;
  particle_sets.resize(0);
  xml_main_input->reset()->get_child("system")->get_child("particles");
  
  do
    {
      xml_main_input->get_attribute("n");
      particle_sets.push_back(new particles_t(xml_main_input->get_int()));  
      xml_main_input->get_next("particles");
    }
  while(xml_main_input->check());
  
  xml_main_input->reset();
  
}
// init the creation of particles in the system
template<class T>
void particles<T>::init(int n)
{
  work.resize(n);
  position.resize(n);
  position_no_pbc.resize(n);
  
}
template<class T>
void particles<T>::gaussian_move(random1* rand_obj,const double &sigma_diff)
{
  int i;
  rand_obj->gaussian(work);
  for(i=0;i<n;i++)
    {
      position_no_pbc[i]= position_no_pbc[i] + sigma_diff * work[i];
    }
}
template<class particles_t>
void all_particles<particles_t>::gaussian_move(random1* rand_obj,const double & sigma_diff)
{
  int i;
  
  for(i=0;i<particle_sets.size();i++)
    {
      particle_sets[i]->gaussian_move(rand_obj,sigma_diff);
    }
    
}
template<class particles_t>
// print the particle positions
void all_particles<particles_t>::print()
{
  int i;
  for(i=0;i<particle_sets.size();i++)
    {
      cout<< "--set "<<i<<"--"<<endl;
      particle_sets[i]->print();
      cout << "----"<<endl;
    }
  
}
// print particles in a set
template<class T>
void particles<T>::print()
{
  int i=0;
  cout<<n<<" particles in set"<<endl;
  for(i=0;i<n;i++)
    {
      cout<<i<<" "<<position[i]<<endl;
    }
}


// resets the drift force for the first derivative
template<class particles_t>
void all_particles<particles_t>::reset_drift_force_derivative()
{
  
  int i;
  for(i=0;i<particle_sets.size();i++)
    {
      
      particle_sets[i]->reset_drift_force_derivative();
    }
}

template class all_particles<particles<double> >;
template class all_particles<particles< complex<double> > >;
template class all_particles<spinors>;
