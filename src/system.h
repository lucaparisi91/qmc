#ifndef SYSTEM_H
#define SYSTEM_H
#include <vector>
#include "mpi.h"
#include "tools.h"
#include <complex>

using namespace std;
class packed_data;
class random1;
class xml_input;


template<class value_t>
class particles
{
 public:
  double m;// mass of the particle
  int n; //number of particles
  vector<double> work;
  vector<double> position;
  vector<double> position_no_pbc;
  
  vector<value_t> drift_force_derivative;
  particles(int n_);
  particles();
  void resize(int n_);
  
  int get_pack_size();
  
  inline void reset_drift_force_derivative()
  {
    fill(drift_force_derivative.begin(),drift_force_derivative.end(),0);
  }
  
  double& operator[](size_t i){return position[i];}

  const double& operator[](size_t i) const {return position[i];}

  double& getNoBC(size_t i){return position_no_pbc[i];}
  
  const double& getNoBC(size_t i) const {return position_no_pbc[i];}
  
  
  int getN() const {return position.size();}; // returns the  number of particles
  
  void set_all_positions(double pois);
  particles& operator=(particles & p);
  void send_mpi(int dest,int tag,vector<MPI_Request> &requests);
  void receive_mpi(int src,int tag,vector<MPI_Request> &requests);
  void move(random1 * rand_obj,const double & delta_tau); // makes the QMC move
  void pack(packed_data*);
  void unpack(packed_data*);
  void gaussian(random1* rand_obj); // fills with a gaussian
  void gaussian_move(random1* rand_obj,const double &sigma_diff);
  void save(xml_input* save_xml); // adds a particles node
  void load(xml_input* load_xml);
  void init(int n);
  void print();
  void set_gaussian(random1* randg,double alpha,double position);
  void move(vector<double> &displacement);
  void set_uniform(double pos,random1* randg);
};

class spinorsRealDmc : public particles<double>
{
public:
  
  typedef tinyVector<double,2> spinComp_t;

  spinorsRealDmc(int n_) : particles<double>::particles(n_)
  {
    spinComp.resize(n_); 
  }

  spinorsRealDmc& operator=(spinorsRealDmc & s)
  {
    particles<double>::operator=(s);
    spinComp=s.spinComp;
    
  };
  
public:
  
  vector<spinComp_t> spinComp;
};

class spinors : public particles<complex<double> >
{
public:
  spinors(int n_) : particles<complex<double> >::particles(n_)
  {
    int i;
    vector<complex<double> > a(2);
    
    for(i=0;i<n_;i++)
      {
	spinComp.push_back(a);
	spinDrift.push_back(0);
      };
    
  };
  
  spinors& operator=(spinors & s)
  {
    particles<complex<double> >::operator=(s);
    spinComp=s.spinComp;
  };
  
  vector<complex<double> > spinDrift; // vector of the drift force
  vector<vector < complex<double> > > spinComp;
  void save(xml_input* xml_save_particles);
  void load(xml_input* xml_save_particles);
  
};

// the set of all particles in a wavefunction
template<class T>
class all_particles
{
  
 public:
  typedef T particles_t; 
  vector<particles_t*> particle_sets;
  void send_mpi(int dest,int tag);
  void receive_mpi(int src,int tag);
  int get_pack_size();
  void move(random1* rand_obj,const double & delta_tau);// makes the dmc move
  void gaussian_move(random1* rand_obj,const double & sigma_diff);
  all_particles<particles_t>& operator=(all_particles<particles_t> & p);
  all_particles();
  void gaussian(random1* rand_obj);
  void save(xml_input* xml_save); // adds a particle node
  void load(xml_input* xml_save);// loads particles
  void pack(packed_data* packed_data);
  void unpack(packed_data* packed_data);
  void init(xml_input* xml_input); // allocates the initial required size
  void print();
  void reset_drift_force();
  void reset_drift_force_derivative();
  void set_gaussian(random1* randg,double alpha,double position);
  void set_all_positions(double pois);
  void set_uniform(double pos,random1* rang);
  int getN();
  void getNs(vector<int> &ns);
  // returns the total number of particles
  int getNTot()
  {
    int n;
    n=0;
    for(int i=0;i<getN();i++)
      {
	n+=(*this)[i].getN();
      }
    
    return n;
  }
  
  particles_t & operator[](size_t i){return *particle_sets[i];};
  const particles_t & operator[](size_t i) const {return *particle_sets[i];};
  
  
  
private:
};
 
#endif
