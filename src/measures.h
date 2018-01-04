#ifndef MEASURES_H
#define MEASURES_H

#include <vector>
#include <string>
#include "traits.h"
#include "mesh.h"
#include "circular_vector.h"
#include <iostream>
#include "xml-input.h"
#include <complex>
#include "exception_classes.h"
#include "ptools.h"
#include <fstream>

using namespace std;
class xml_input;
class mesh;
class measure_scalar;
class measure_vector;

// build a certain kind of object from ground
class factory
{
  template<class m>
  m* makeMeasurement(xml_input* mNode); // performs a measurement on the system
};


class counter
{
 public:
  int i;
  bool termal;
  int jumps;
  int skip;
  int status;
  counter(int jumps,int skip); //a counter to count the number of steps
  counter& operator()(counter & c);
  void increment(); //increments the counter
  int check();
  void set(int jumps,int skip);
};

// interface to measurements
template<class T>
class estimator_decorrelator
{
public:
  estimator_decorrelator();
  void add(const T& m);
  void add(const T& m,const int& i);
  void variances();
  void out(const string &filename);
  void save(const string &filename);
  void load(const string &filename);
  void check_variances();
  T get_error(){return error; }
  bool is_converged(){return converged;}
private:
  double eps;//acceptaince ratio for the convergence
  int depth; // the depth of the binary tree
  T zero_el;
  vector<T> sums; // mean values
  vector<T> sum_squares;// sum of the squares
  vector<int> n;
  vector<T> last;
  vector <int> waiting;
  T error;
  bool converged;
  vector<T> vars;
  
};

template<class T>
class estimator_decorrelator_vector
{
public:
  estimator_decorrelator_vector(int n);
  void add(const T& m,const int &j);
  void add(const vector<T>& m);
  void add(const T& m,const int& j,const int&i);
  void add(const vector<T> & m,int n)
{
  int i;
  for(i=0;i<m.size();i++)
    {
      add(m[i]/n,i);
    }
}

  void variances();
  void out(const string &filename);
  void save(const string &filename);
  void load(const string &filename);
  void check_variances();
  T get_error(int i){return errors[i];};
  bool is_converged(int i){return converged[i];}
private:
  vector<int> depth; // the depth of the binary tree
  T zero_el;
  vector<vector<T> > sums; // mean values
  vector<vector<T> > sum_squares;// sum of the squares
  vector<vector<int> > n;
  vector<vector<T> > last;
  vector<vector<int> > waiting;
  vector<vector<T> > vars;
  vector<T> errors;
  vector<bool> converged;
  T eps;
  
};

class measure
{
 public:
  string label;
  int set_a;
  int n_block_jump;
  int n_blocks;
  int set_b; // used only if required in two body measurements
  virtual void set_history(bool hist)=0;
  virtual void out()=0;
  virtual void clear_history()=0;
  virtual void record(double step)=0;
  virtual void print()=0;
  virtual void load()=0;
  
  virtual void clear()=0;
  
  virtual void save()=0;
  virtual void increment()=0;
  virtual void set(int jumps_,int skip_)=0;
  virtual int check()=0;
  virtual void add(double e , double step){cout<<"No add."<<endl;};
  virtual void reduce(int root)=0;

};

class measure_scalar : public measure
{
 public:
  int n;
  double n_t;
  double last_step;
  double sum_m_t;// keeps a mean not reset by a block
  double sum_m; //keeps a mean value of the measurements
  double error_m; // current estimate on the error of the measurements
  counter* c;
  vector<double> recorded_measurements;
  vector<double> recorded_times;
  measure_scalar& operator()( const     measure_scalar & m2);
  int n_recorded;
  measure_scalar(string label_);
 public:
  void add(double m, double step);
  void clear();//erase all information
  void clear_history();
  bool history;
  int check();
  void increment();
  void set_history(bool hist);
  void out_history();
  void out();//write to file the last measurement
  void out_t();
  void set(int jumps_,int skip_);
  void reduce(int root);
  double mean();
  void save();
  void print();
  void load();
  double get_error(){return dec.get_error();}
  bool is_converged(){return dec.is_converged();}
  virtual void record(double step);// record the mean value and the last time step in an array
  double mean_t();
  //void average_with(measure* ms);
  
  protected:
  estimator_decorrelator<double> dec;
 };
// continuosly record the time measurements
// a vector of measurement

class measure_vector : public measure
{
public:
  int len;
  counter* c;
  void increment();
  int check();
  double last_step;
  bool history;
  int n;
  bool vector_n;
  int iBlock;
  int n_tot;
  vector<double> sum;// used for the mean sum
  vector<double> sum_tot;
  vector<double> sum2_tot;//used for the mean squared sum over blocks
  vector<double> sum_block; // performs the sum over all blocks
  
  //int n_block_jump;// number of blocks over which to compute the variances
  vector<double> getMean()
  {
    vector<double> res;
    res=sum;
    for(unsigned int i=0;i<res.size();i++)
      {
	res[i]=res[i]/n;
      }
    
    return res;
  }
  
  void out();
  void reduce(int root);
  void out_t();
  void record(double step);
  void clear();
  measure_vector& operator()(const measure_vector &m2);
  void clear_history();
  void print();
  void set(int jumps,int skip);
  void add(vector<double> &vec,double step);// add a whole vector
  void add(double value,int i,double step);
  void increment_value(double value,int i,double step);
  void increment_index();
  void increment_index(int n2){n+=n2;}
  void save();
  void set_history(bool hist);
  void load();
  
  vector<double> mean();
  vector<double> mean_t();
  measure_vector(int len_,string label_,bool history_);
protected:
  estimator_decorrelator_vector<double> dec;
};

vector<double> build_q_vector(int bins,double l_box,double qMax);

// add a measurement on the pair correlation vector
template<class mt >
class space_measure : public mt
{
 public:
  mesh *grid;
  void add(double x,double value,double step);
  void add(double value,double x,double w,double step){mt::add(value,grid->index(x),w,step);}
  space_measure(double a,double b,int len_,string label_,bool history_);
  void increment_value(double value,double x,double step);
  
  double get_step(){return grid->get_step();}
  
};


template<class mode>
class outMode
{
public:
  void out(vector<double> &ms,vector<double> &ns,string filename)
  {
    ofstream f;
    unsigned int i=0;
    f.open(filename.c_str(),std::ios_base::app);
    
    for(i=0;i<ms.size();i++)
      {
	if (ns[i]>0)
	  {
	    f << i <<" "<< ms[i]/ns[i]<< endl;
	    //cout << ns[i]<<endl;
	  }
      }
  }
  
  void out(double &s,double &ns,string filename)
  {
    ofstream f;
    int i=0;
    f.open(filename.c_str(),ios::app);
    f << s/ns<<endl;
    f.close();
  }
  
};

class measure_vector_mult_index : public measure_vector
{
 public:
  measure_vector_mult_index(int len_,string label_,bool history_);
  vector <double> ns;
  vector<double> ns_block; 
  vector<double> ns_tot;
  vector<double> nb;
  vector<double> sum_of_blocks;
  void increment_index(int i);
  void increment_index();
  void increment_weight(double weight);
  void increment_index(int i,double weight){ns[i]+=weight;}
  vector<double> mean();
  void reduce(int root); // makes a reduction operation
  void clear();
  void add(double x,int i,double step){ns[i]++;sum[i]+=x;}
  void add(double x,int i,double w,double step){ns[i]+=w;sum[i]+=x*w;}
  void save();
  void load();
  void out_t();
  void out(){if (history) {outObj.out(sum,ns,label + ".dat") ;} else    {out_t();} };
  void record_block(double time_step);
  
  outMode<keep_history_t> outObj;
  
};

// dynamic measurement on the system

template<class qt>
class measure_dynamic
{
  typedef typename qt::all_particles_t all_particles_t;
  typedef typename qt::wave_t wave_t;
  typedef typename qt::walker_t walker_t;
 public:
  //virtual void make_measurement(all_particles* state,total_wavefunction* wave)=0;
  measure_dynamic(){set_a=0;};
  measure_dynamic(const int &set){set_a=set;}
  virtual void print()=0;
  virtual void reset()=0;
  virtual int isFilled() {throw notYetSsupported("isFilled");};
  virtual void pack(packed_data* packed_data)=0;
  virtual void unpack(packed_data* packed_data)=0;
  virtual int get_pack_size()=0;
  virtual measure_dynamic<qt>& operator=(measure_dynamic<qt>&)=0;
  int set_a;
  virtual void time_difference_average(all_particles_t * state,wave_t* wave,vector<double> &sum,vector<double> &ns,double mes)
  {
    throw notYetSsupported("measure_dynamic::time_difference_average");
  }
    
  virtual double average(){throw notYetSsupported("distance_center_of_mass_d::time_difference_average");};
  
  virtual void average(vector<double> &){throw notYetSsupported("distance_center_of_mass_d::time_difference_average");};
  
  virtual void add(const double &){throw notYetSsupported("Adding a scalar is not supported");};

  virtual void add(vector<double>&){throw notYetSsupported("Adding a vector is not supported");};

  virtual vector<double>& currentVector(){throw notYetSsupported("Getting a vector is not supported");};

  virtual void incrementIndex(){throw notYetSsupported("increent index");};

  
  
  virtual int get_n(){throw notYetSsupported("get_n is not supported");};

  virtual void make_measurement(walker_t* w,wave_t* wave){throw notYetSsupported("No measurement supported");};
  
};

// dynamic measurements
template<class comp>
class measure_dynamic_scalar : public measure_dynamic<comp>
{
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename comp::wave_t wave_t;
  
 public:
  
  virtual void pack(packed_data* packed_data);
  virtual void unpack(packed_data* packed_data);
  virtual int get_pack_size();
  measure_dynamic_scalar(const int &bins_,const int &set);
  
  virtual double average()
  {
    int i;
    double sum;
    sum=0;
   
    for(i=0;i<this->ms.size();i++)
      {
	
	sum+=this->ms[i];
      }
    
    return sum/this->ms.size();  
  }
  
  virtual void time_difference_average(all_particles_t * state,wave_t* wave,vector<double> &sum,vector<double> &ns,double mes)
  {
    int j;

    if (this->filled)
      {
	this->reset();
      }
    
    this->add(mes);
  
    for (j=0;j<=this->last;j++)
      {
	sum[j]=sum[j] + pow(this->ms[this->last] - this->ms[this->last - j],2);
	ns[j]=ns[j] + 1;
	
      }
    
    
  }
  
  
  
  int last;
  int bins;
  bool filled; // whatever all the bins have been filled
  vector<double> ms;//measurements
  //void make_measurement(all_particles* state,total_wavefunction* wave);//
  //measure_dynamic_scalar<T>& operator=(measure_dynamic_scalar<T>&);
  void add(double m1);
  void clear();
  void reset();
  void print();
  bool isFilled(){return filled;};
  measure_dynamic<comp>& operator=(measure_dynamic<comp> &);
  //void clear();
  
};
// compute the difference between center of masses

#include "storage.h"

template<class comp,class storage_t>
class futureWalker : public measure_dynamic<comp>
{
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename comp::wave_t wave_t;
  typedef typename storage_t::scalar_t scalar_t;
  
public:
  
  void add(scalar_t m1)
  {
    storage+=m1;
  }
  
  futureWalker() : measure_dynamic<comp>(0)
  {
    
  }
  futureWalker(int ns) : measure_dynamic<comp>(0),storage(ns)
  {
    
  }
  
  template<class mesh_t>
  futureWalker(int ns,mesh_t mesh) : measure_dynamic<comp>(0),storage(ns,mesh)
  {
    
  }
  

  
  virtual void unpack(packed_data* md_pack)
  {
    storage.unpack(md_pack);
  }
  virtual void pack(packed_data* md_pack)
  {
    storage.pack(md_pack);
  }
  virtual int get_pack_size()
  {
    return storage.get_pack_size();
  }
  
  virtual void reset()
  {
    storage.reset();
  }
  
  int isFilled()
  {
    return storage.isFilled();
  }

  void increment_index()
  {
    storage.increment_index();
  }
  
  measure_dynamic<comp>& operator=(measure_dynamic<comp> & m2)
  {
    
    storage.clone(static_cast<futureWalker<comp,storage_t>  & >(m2).storage);
  }
  
  // get the number of particles
  
  
  // get the sum of the measurements
  scalar_t get_sum(){return storage.get_sum();}
  // get the mean over particles
  int get_n(){return storage.get_n();}
  void clear(){storage.clear();};
  
protected:
  // contains the data to be copied
  storage_t storage;
  
};

template<class comp,class storage_t>
class futureWalkerScalar : public futureWalker<comp,storage_t>
{
public:
  futureWalkerScalar() : futureWalker<comp,storage_t>()
  {
    
  }
  
  double average(){return this->storage.average();};

  void print()
  {
    if(this->storage.get_n()>0)
      {
	cout <<"value: "<< this->storage.average()<<endl;	
      }
  }
  
};




template<class comp,class T>
class futureWalkerVector : public futureWalker<comp,vectorStorage<T> >
{
  
public:
  typedef vectorStorage<T> storage_t;
  
  futureWalkerVector(int ns) : futureWalker<comp,storage_t>::futureWalker(ns)
  {
    
  }
  void print()
  {
    cout << "vectorFutureWalker"<<endl;
  }
  void increment_value(T x,int j,double step)
  {
    this->storage.increment_value(x,j,step);
  }
  
  void average(vector<T> & v)
  {
    this->storage.average(v);
  }
  
};

template<class comp,class T,class mesh_t>
class futureWalkerSpaceVector : public futureWalker<comp,vectorSpaceStorage<T,mesh_t> >
{
  
public:
  typedef vectorSpaceStorage<T,mesh_t> storage_t;
  
  futureWalkerSpaceVector(int ns,mesh_t m_): futureWalker<comp,storage_t>::futureWalker(ns,m_)
  {
    
  }
  
  void increment_value(T x,T y,double step)
  {
    this->storage.increment_value(x,y,step);
  }
  
  double get_step(){return this->storage.get_step();}

  void average(vector<T> & v)
  {
    this->storage.average(v);
  }

  void print()
  {
    cout << "vectorSpaceFutureWalker"<<endl;
  }
};


template<class comp,class T>
class structureFactorWalker : public futureWalkerVector<comp,T >
{
public:
  typedef T measure_scalar_t;
  typedef typename comp::walker_t walker_t;
  typedef typename comp::wave_t wave_t;
  
  structureFactorWalker(vector<double> qs_,int setA_) : futureWalkerVector<comp,double >::futureWalkerVector(qs_.size())
  {  
    qs=qs_;
    setA=setA_;
    work.resize(qs.size());
  }
  
  void make_measurement(walker_t* w,wave_t* wave)
  {
    
    wave->structure_factor(w->state,&(this->storage),qs,setA,work);
    
  };
  
private:
  vector<T> qs;
  int setA;
  vector<complex<double> > work;
};


template<class comp,class T>
class structureFactorDoubleWalker : public futureWalkerVector<comp,T >
{
public:
  typedef T measure_scalar_t;
  typedef typename comp::walker_t walker_t;
  typedef typename comp::wave_t wave_t;
  
  structureFactorDoubleWalker(vector<double> qs_,int setA_,int setB_) : futureWalkerVector<comp,double >::futureWalkerVector(qs_.size())
  {  
    qs=qs_;
    setA=setA_;
    setB=setB_;
    work.resize(qs.size());
  }
  
  void make_measurement(walker_t* w,wave_t* wave)
  {
    
    wave->structure_factor(w->state,&(this->storage),qs,setA,setB,work);
    
  };
  
private:
  vector<T> qs;
  int setA;
  int setB;
  vector<complex<double> > work;
};

template<class comp,class T,class mesh_t>
class pairCorrelationSymmFuture : public futureWalkerSpaceVector<comp,T,mesh_t >
{
public:
  typedef T measure_scalar_t;
  typedef typename comp::walker_t walker_t;
  typedef typename comp::wave_t wave_t;
  
  pairCorrelationSymmFuture(mesh_t grid,int setA_,int bins,double max_) : futureWalkerSpaceVector<comp,double,mesh_t>::futureWalkerSpaceVector(bins,grid)
  {
    setA=setA_;
    max=max_;
  }
  
  void make_measurement(walker_t* w,wave_t* wave)
  {
    wave->pair_correlation_symm(w->state,&(this->storage),max,setA); 
  };

  
  
  
private:
  double max;
  int setA;
  
};

template<class comp,class T,class mesh_t>
class pairCorrelationAsymmFuture : public futureWalkerSpaceVector<comp,T,mesh_t >
{
public:
  typedef T measure_scalar_t;
  typedef typename comp::walker_t walker_t;
  typedef typename comp::wave_t wave_t;
  
  pairCorrelationAsymmFuture(mesh_t grid,int setA_,int setB_,int bins,double max_) : futureWalkerSpaceVector<comp,double,mesh_t>::futureWalkerSpaceVector(bins,grid)
  {
    setA=setA_;
    max=max_;
    setB=setB_;
  }
  
  void make_measurement(walker_t* w,wave_t* wave)
  {
    wave->pair_correlation_asymm(w->state,&(this->storage),max,setA,setB);
  };

  
  
private:
  double max;
  int setA;
  int setB;
};

template<class comp,class T>
class structureFactorSpinWalker : public futureWalkerVector<comp,T >
{
public:
  typedef T measure_scalar_t;
  typedef typename comp::walker_t walker_t;
  typedef typename comp::wave_t wave_t;
  
  structureFactorSpinWalker(vector<double> qs_,int setA_,int setB_) : futureWalkerVector<comp,double >::futureWalkerVector(qs_.size())
  {  
    qs=qs_;
    setA=setA_;
    setB=setB_;
    work.resize(qs.size());
  }
  
  void make_measurement(walker_t* w,wave_t* wave)
  {
    
    wave->structure_factorSpin(w->state,&(this->storage),qs,setA,setB,work);
    
  };
  
private:
  vector<T> qs;
  int setA;
  int setB;
  vector<complex<double> > work;
};

// a vector measurement from scalar measurement

// template for the parameters
template<class comp>
class center_of_mass_w : public measure_dynamic_scalar<comp>
{
  typedef typename comp::all_particles_t all_particles_t;
  typedef typename comp::wave_t wave_t;
 public:


  center_of_mass_w(const int &bins,const int &set) : measure_dynamic_scalar<comp>(bins,set){};
  void clear();

  
  
};

template<class comp>
void append_walker_measures(vector<measure_dynamic<comp>*> &md,vector<measure*> &ms,xml_input* xml_md,bool only_walker);
// applies and stores everything in a certain slide

// the interface from which measurement are accessed
template<class walker_t,class wave_function_t>
class measurementInterface
{
public:
  virtual void save()=0;
  virtual void out()=0;
  virtual void clear()=0;
  virtual int check()=0;
  virtual void print()=0;
  virtual void reduce(int i)=0;
  virtual void clear_history()=0;
  virtual void load()=0;
  virtual void increment()=0;
  virtual void record(double step)=0;
  // perfrom the measurement
  virtual void make_measurement(walker_t* w,wave_function_t* wave)=0;
  
};



class rqmc_walker;
template<class wave_t>
class rqmcMeasurementInterface
{
public:
  virtual void save()=0;
  virtual void out()=0;
  virtual void make_measurement(vector<rqmc_walker*> p,wave_t* w)=0;
  virtual void load()=0; // load a certain wavefunction for the system
  virtual void set_slices(int i)=0;
};

template<class wave_t,class m>
class rqmcMeasurement : public rqmcMeasurementInterface<wave_t>
{
public:
  rqmcMeasurement(m* m1,int slices_);
  virtual void make_measurement(vector<rqmc_walker*> p,wave_t* w);
  virtual void save(){for(int i=0;i<slices;i++){ms[i]->save();};};
  virtual void out(){for(int i=0;i<slices;i++){ms[i]->out();};};
  virtual void set_slices(int i){slices=i;ms.resize(i);}; // set the  number of slices for the operator
  virtual void load(){for(int i=0;i<slices;i++){ms[i]->load();};}
  
protected:
  vector<m*> ms; // a vector of pointer to single measurements
  int slices;
  
};

// a partial implementation for all kind of measurements. To be inherthed by single measurement classes (energy, density...)
template<class walker_t,class wave_function_t,class estimator_t>
class measurement : public measurementInterface<walker_t,wave_function_t>
{
public:
  measurement(estimator_t* ms_){ms=ms_;};
  measurement(){};
  measurement<walker_t,wave_function_t,estimator_t>& operator()(const measurement<walker_t,wave_function_t,estimator_t> &m)
  {
    ms(*m);
  };
  
  virtual void print(){ms->print();};
  virtual void save(){ms->save();};
  virtual void out(){ms->out();};
  virtual int  check(){return ms->check();};
  virtual void clear(){ms->clear();};
  virtual void clear_history(){ms->clear_history();};
  virtual void record(double step){ms->record(step);};
  virtual void load(){ms->load();};
  virtual void reduce(int i){ms->reduce(i);};
  virtual void increment(){ms->increment();};
  virtual void make_measurement(walker_t* w,wave_function_t* wave)=0;
protected:
  estimator_t* ms;  
};




// actually performs the measurements
template<class walker_t,class wavefunction_t>
class measurement_scalar : public measurement<walker_t,wavefunction_t,measure_scalar>
{
public:
  measurement_scalar(string label){this->ms=new measure_scalar(label);}
  measurement_scalar(measure_scalar* es_) : measurement<walker_t,wavefunction_t,measure_scalar>::measurement(es_) {};
  
};

// get the DMC energy from a value stored in the walker
template<class walker_t,class wave_t>
class get_energy_dmc : public measurement_scalar<walker_t,wave_t>
{
public:
  
  get_energy_dmc(string label) :  measurement_scalar<walker_t,wave_t>::measurement_scalar(label){};
  void make_measurement(walker_t* w,wave_t* wave){this->ms->add(w->e,0); }
};

template<class walker_t,class wave_t>
class get_f_energy_dmc : public measurement_scalar<walker_t,wave_t>
{
public:
  get_f_energy_dmc(string label) :  measurement_scalar<walker_t,wave_t>::measurement_scalar(label){};
  void make_measurement(walker_t* w,wave_t* wave){this->ms->add(w->e_f,0);};
  
};

template<class walker_t,class wave_t>
class center_of_mass : public measurement_scalar<walker_t,wave_t>
{
public:
  center_of_mass(measure_scalar* mScal) :  measurement_scalar<walker_t,wave_t>::measurement_scalar(mScal){};
  void make_measurement(walker_t* w,wave_t* wave){this->ms->add(wave->center_of_mass_no_pbc(w->state,this->ms->set_a),0);};
private:
  
};

template<class walker_t,class wave_t>
class center_of_mass_difference : public measurement_scalar<walker_t,wave_t>
{
public:
  center_of_mass_difference(measure_scalar* mScal) :  measurement_scalar<walker_t,wave_t>::measurement_scalar(mScal){};
  void make_measurement(walker_t* w,wave_t* wave){this->ms->add(wave->center_of_mass_no_pbc(w->state,this->ms->set_a) - wave->center_of_mass_no_pbc(w->state,this->ms->set_b) ,0);};
  
};

template<class walker_t,class wave_t>
class center_of_mass_difference2 : public measurement_scalar<walker_t,wave_t>
{
public:
  center_of_mass_difference2(measure_scalar* mScal) :  measurement_scalar<walker_t,wave_t>::measurement_scalar(mScal){};
  void make_measurement(walker_t* w,wave_t* wave)
  {
    this->ms->add(
		  pow(   wave->center_of_mass_no_pbc(w->state,this->ms->set_a) - wave->center_of_mass_no_pbc(w->state,this->ms->set_b)
			 ,2),0);
  }
  
};



template<class walker_t,class wave_function_t>
class density : public measurement<walker_t,wave_function_t,space_measure<measure_vector> >
{
public:
  density(space_measure<measure_vector>* ms_) : measurement<walker_t,wave_function_t,space_measure<measure_vector> >(ms_) {}; 
  virtual void make_measurement(walker_t* w,wave_function_t* wave){wave->density(w->state,this->ms);};
};

/*
template<class walker_t,class wave_function_t>
class hessian : public measurement<walker_t,wave_function_t,measure_vector >
{
public:
  hessian(measure_vector* ms_) : measurement<walker_t,wave_function_t,measure_vector>(ms_)
  {
    hessianAndSMatrix.resize(9);
    wavefunctionsToOptimize.resize(1);
    wavefunctionsToOptimize[0]=1;
  };
  
  virtual void make_measurement(walker_t* w,wave_function_t* wave)
  {
    double eD;
    double e;
    double phi;
    double phi2;
    
    wave->kinetic_energy_derivative(w->state,eD,wavefunctionsToOptimize,phi,phi2);
    
    e=w->e;
    
    //eD=eD+w->e;
    //cout << phi*(eD-e)<<endl;
    //cout << phi << endl;
    hessianAndSMatrix[0]=e;
    hessianAndSMatrix[1]=e*phi;
    hessianAndSMatrix[2]=phi;
    hessianAndSMatrix[3]=phi*(eD);
    hessianAndSMatrix[4]=phi*phi*e;
    hessianAndSMatrix[5]=phi*phi*eD;
    hessianAndSMatrix[6]=phi*phi;
    hessianAndSMatrix[7]=phi2;
    hessianAndSMatrix[8]=phi2*e;
    
    //print_vector(hessianAndSMatrix);
    //cout << endl;
    //cout <<w->e<<" "<<eD << " "<<phi << endl;
    //exit(0);
    // accumulates the matrix elements
    this->ms->add(hessianAndSMatrix,0);
    //cout << eD<<" "<<phi<<endl;
  };
  
  virtual void out()
  {
    
  };
  int getWavefunctionIndex(int j)
  {
    return wavefunctionsToOptimize[j];
  }
  double estimateStep()
  {
    double step;
    vector<double> HS;
    HS=this->ms->getMean();
    //print_vector(HS);
    step=optO.estimateStep( HS );
    cout << step << endl;
    return step;
    
  }
  
  double getEnergy()
  {
    return this->ms->getMean()[0];
  }
  
  
private:  
  vector<double> hessianAndSMatrix;
  vector<int> wavefunctionsToOptimize;
  NewtonOptimizer<1> optO;
};

*/

template<class walker_t,class wave_t>
class oneBodyDensityMatrixOffdiagonal : public measurement<walker_t,wave_t,space_measure<measure_vector_mult_index> >
{
public:
  oneBodyDensityMatrixOffdiagonal(space_measure<measure_vector_mult_index>* ms_,int nMCM_,double Max_) : measurement<walker_t,wave_t,space_measure< measure_vector_mult_index> >(ms_){nMCM=nMCM_;Max=Max_;};
  
  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    wave->oneBodyDensityMatrixOffdiagonal(w->state,this->ms,nMCM,Max);
    
  }
  
private:
  int nMCM;
  double Max;
};

// is a small DMC driver for the system
template<class tm>
class miniDMCDriver
{
  typedef typename tm::walker_t walker_t;
  typedef typename tm::all_particles_t all_particles_t;
public:
  miniDMCDriver(tm* dmc_obj_,int nEvolve_)
  {
    dmc_obj=dmc_obj_;
    w=new walker_t(dmc_obj);
    reset_weight();
    nEvolve=nEvolve_; //  number of evolving objects  
  }
  
  void set(walker_t* w2)
  {
    // copy the value of the second walker
    (*w)=(*w2);
  }
  // evolve over a certain number of steps
  inline void evolve()
  {
    evolve(nEvolve);
  }
  
  void set(all_particles_t* p)
  {
    w->set_particles(p);
    w->set(dmc_obj);
    
  }
  
  
  void setPosition(int set,int i,double rnew)
  {
    w->state->particle_sets[set]->position_no_pbc[i]=rnew;
    w->state->particle_sets[set]->position[i]=rnew;
  }
  void setPosition(int set,int i,double rnewNoPbc,double rnewPbc)
  {
    w->state->particle_sets[set]->position_no_pbc[i]=rnewNoPbc;
    w->state->particle_sets[set]->position[i]=rnewPbc;
  }
  void evolve(int n)
  {
    int i;
    
    for(i=0;i<n;i++)
      {
	// update the phantom walker
	w->update(dmc_obj);
	// accumulate the weight of the walker
	weight*=w->get_weight(dmc_obj->e_t,dmc_obj);
	
      }
    
    
  }
  
  // returns the wight of a certain configuration
  inline double get_weight()
  {
    return weight;
  }
  
  // reset the weight
  inline void reset_weight()
  {
    weight=1;
  }
  
private:
  int nEvolve;
  tm *dmc_obj;
  walker_t* w;
  double weight;
};


template<class tm>
class oneBodyDensityMatrixOffdiagonalFutureWalkers : public measurement< typename tm::walker_t,typename tm::wave_t,space_measure<measure_vector_mult_index> >
{
  typedef typename tm::walker_t walker_t;
  typedef typename tm::wave_t wave_t;
  
public:
  oneBodyDensityMatrixOffdiagonalFutureWalkers(space_measure<measure_vector_mult_index>* ms_,int nMCM_,double Max_,int nFutureWalkers_,tm* qmc_obj ) : measurement<walker_t,wave_t,space_measure< measure_vector_mult_index> >(ms_),driver(qmc_obj,nFutureWalkers_){nMCM=nMCM_;Max=Max_;nFutureWalkers=nFutureWalkers_;};
  
  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    wave->oneBodyDensityMatrixOffdiagonalFutureWalkers(w->state,this->ms,nMCM,Max,&driver);
  }
private:
  miniDMCDriver<tm> driver;
  int nMCM;
  double Max;
  int nFutureWalkers; // number of walkers in the future
};

// measurement of pair correlations
template<class walker_t,class wave_t>
class pair_correlation_m : public measurement<walker_t,wave_t,space_measure<measure_vector> >
{
public:
  pair_correlation_m(space_measure<measure_vector>* ms_,int setA_,int setB_) : measurement<walker_t,wave_t,space_measure<measure_vector> >(ms_),setA(setA_),setB(setB_) {};

  virtual void make_measurement(walker_t* w,wave_t* wave);
private:
  int setA;
  int setB;
  
};

// compute the symmetri structure factor
template<class walker_t,class wave_t>
class structure_factor_symm : public measurement<walker_t,wave_t,measure_vector >
{
public:
  
  structure_factor_symm(measure_vector* ms_,int bins,double deltaQ,double l_box,int setA_) : measurement<walker_t,wave_t,measure_vector >(ms_)
  {
    int i;
    // change the number of bins
    qs.resize(bins);
    qs[0]=2*M_PI/l_box;
    for(i=1;i<bins;i++)
      {
	qs[i]=qs[0]+i*deltaQ;
      }
    
    setA=setA_;
    
  };

  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    wave->structure_factor_symm(w->state,this->ms,qs,setA);
  }
private:
  
  vector<double> qs;
  int setA;
};

template<class walker_t,class wave_t>
class structure_factor_single_f : public measurement<walker_t,wave_t,measure_vector>
{
public:
  
  structure_factor_single_f(measure_vector* m_,int iW_,int nFutureWalkers,int bins) : measurement<walker_t,wave_t,measure_vector>::measurement(m_)
  {
    iW=iW_;
    filled=false;
    nMeasures=nFutureWalkers;
    structureFactorTmp.resize(bins);
  }
  
private:
  int iW;
  int setA;
  int nMeasures;
  int nC;
  vector<double> structureFactorTmp;
  bool filled;
  
  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    // add to the walker the measurement just done
    
    w->md[iW]->make_measurement(w,wave);
    
    if (w->md[iW]->get_n()>=nMeasures)
      {
	
	if (w->md[iW]->isFilled())
	  {
	    
	    w->md[iW]->average(structureFactorTmp);
	    this->ms->add(structureFactorTmp,0);
	    
	  }
	w->md[iW]->reset();
      }
    
   
    
    
  }
};

template<class walker_t,class wave_t>
class structure_factor_single_complex : public measurement<walker_t,wave_t,measure_vector >
{
public:
  
  structure_factor_single_complex(measure_vector* ms_,int bins,double deltaQ,double l_box,int setA_) : measurement<walker_t,wave_t,measure_vector >(ms_)
  {
    int i;
    // change the number of bins
    qs.resize(bins);
    qs[0]=2*M_PI/l_box;
    for(i=1;i<bins;i++)
      {
	qs[i]=qs[0]+i*deltaQ;
      }
    
    setA=setA_;
    work.resize(bins);
    
  };

  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    wave->structure_factor(w->state,this->ms,qs,setA,work);
  }
  
private:
  
  vector<double> qs;
  vector<complex< double> >  work;
  int setA;
  
};

template<class walker_t,class wave_t>
class structure_factor_double_complex : public measurement<walker_t,wave_t,measure_vector >
{
public:
  
  structure_factor_double_complex(measure_vector* ms_,int bins,double deltaQ,double l_box,int setA_,int setB_) : measurement<walker_t,wave_t,measure_vector >(ms_)
  {
    int i;
    // change the number of bins
    qs.resize(bins);
    qs[0]=2*M_PI/l_box;
    for(i=1;i<bins;i++)
      {
	qs[i]=qs[0]+i*deltaQ;
      }
    
    setA=setA_;
    setB=setB_;
    work.resize(bins);
    
  };

  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    wave->structure_factor(w->state,this->ms,qs,setA,setB,work);
  }
  
private:
 
  vector<double> qs;
  vector<complex< double> >  work;
  int setA;
  int setB;
};



template<class walker_t,class wave_t>
class structure_factor_spin_complex : public measurement<walker_t,wave_t,measure_vector >
{
public:
  
  structure_factor_spin_complex(measure_vector* ms_,int bins,double deltaQ,double l_box,int setA_,int setB_) : measurement<walker_t,wave_t,measure_vector >(ms_)
  {
    int i;
    // change the number of bins
    qs.resize(bins);
    qs[0]=2*M_PI/l_box;
    for(i=1;i<bins;i++)
      {
	qs[i]=qs[0]+i*deltaQ;
      }
    
    setA=setA_;
    setB=setB_;
    work.resize(bins);
    
  };

  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    wave->structure_factorSpin(w->state,this->ms,qs,setA,setB,work);
  }
  
private:
 
  vector<double> qs;
  vector<complex< double> >  work;
  int setA;
  int setB;
};

template<class walker_t,class wave_t>
class structure_factor_spin_complexOrbitals : public measurement<walker_t,wave_t,measure_vector >
{
public:
  
  structure_factor_spin_complexOrbitals(measure_vector* ms_,int bins,double deltaQ,double l_box,int setA_,int setB_) : measurement<walker_t,wave_t,measure_vector >(ms_)
  {
    
    int i;
    // change the number of bins
    qs.resize(bins);
    qs[0]=2*M_PI/l_box;
    for(i=1;i<bins;i++)
      {
	qs[i]=qs[0]+i*deltaQ;
      }
    
    setA=setA_;
    setB=setB_;
    work.resize(bins);
    
  };
  
  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    structureFactorSpin((*w->state)[setA],this->ms,qs,work);
  }
  
  
private:
 
  vector<double> qs;
  vector<complex< double> >  work;
  int setA;
  int setB;
};

template<class walker_t,class wave_t>
class structure_factor_spin_complexOrbitalsForwardWalking : public measurement<walker_t,wave_t,measure_vector >
{
public:
  
  structure_factor_spin_complexOrbitalsForwardWalking(measure_vector* ms_,int bins,double deltaQ,double l_box,int setA_,int setB_,int indexStorage_) : measurement<walker_t,wave_t,measure_vector >(ms_)
  {
    
    int i;
    // change the number of bins
    qs.resize(bins);
    qs[0]=2*M_PI/l_box;
    for(i=1;i<bins;i++)
      {
	qs[i]=qs[0]+i*deltaQ;
      }
    
    setA=setA_;
    setB=setB_;
    work.resize(bins);
    indexStorage=indexStorage_;
    
  };
  
  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    if (w->md[indexStorage]->isFilled())
      {
	this->ms->add(w->md[indexStorage]->currentVector(),0);
      }
    
    structureFactorSpin((*w->state)[setA],w->md[indexStorage]->currentVector(),qs,work);
    w->md[indexStorage]->incrementIndex();
    
    
  }
  
private:
  
  int indexStorage;
  vector<double> qs;
  vector<complex< double> >  work;
  int setA;
  int setB;
};


template<class walker_t,class wave_t>
class structure_factor_density_complexOrbitalsForwardWalking : public measurement<walker_t,wave_t,measure_vector >
{
public:
  
  structure_factor_density_complexOrbitalsForwardWalking(measure_vector* ms_,int bins,double deltaQ,double l_box,int setA_,int setB_,int indexStorage_) : measurement<walker_t,wave_t,measure_vector >(ms_)
  {
    
    int i;
    // change the number of bins
    qs.resize(bins);
    qs[0]=2*M_PI/l_box;
    for(i=1;i<bins;i++)
      {
	qs[i]=qs[0]+i*deltaQ;
      }
    
    setA=setA_;
    setB=setB_;
    work.resize(bins);
    indexStorage=indexStorage_;
    
  };
  
  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    if (w->md[indexStorage]->isFilled())
      {
	this->ms->add(w->md[indexStorage]->currentVector(),0);
      }
    
    structureFactorDensity((*w->state)[setA],w->md[indexStorage]->currentVector(),qs,work);
    w->md[indexStorage]->incrementIndex();
    
  }
  
private:
  
  int indexStorage;
  vector<double> qs;
  vector<complex< double> >  work;
  int setA;
  int setB;
};


template<class walker_t,class wave_t>
class structure_factor_asymm : public measurement<walker_t,wave_t,measure_vector >
{
public:
  
  structure_factor_asymm(measure_vector* ms_,int bins,double deltaQ,double l_box,int setA_,int setB_) : measurement<walker_t,wave_t,measure_vector >(ms_)
  {
    int i;
    // change the number of bins
    qs.resize(bins);
    qs[0]=2*M_PI/l_box;
    
    for(i=1;i<bins;i++)
      {
	qs[i]=qs[0]+i*deltaQ;
      }
    
    setA=setA_;
    setB=setB_;
  };

  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    wave->structure_factor_asymm(w->state,this->ms,qs,setA,setB);
  }
  
private:
  vector<double> qs;
  int setA;
  int setB;
};

template<class walker_t,class wave_t>
class center_of_mass_difference_future_walker : public measurement_scalar<walker_t,wave_t>
{
public:
  center_of_mass_difference_future_walker(measure_scalar* m_,int setA_,int setB_,int iW_,int nFutureWalkers) : measurement_scalar<walker_t,wave_t>::measurement_scalar(m_){iW=iW_;setA=setA_;setB=setB_;filled=false;nMeasures=nFutureWalkers;nC=0;}
  
private:
  int iW;
  int setA;
  int setB;
  int nMeasures;
  int nC;
  bool filled;
  
  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    // add to the walker the measurement just done
    w->md[iW]->add(
		    wave->center_of_mass_no_pbc(w->state,setA) - wave->center_of_mass_no_pbc(w->state,setB));
    
    if (w->md[iW]->get_n()>=nMeasures)
      {
	if (w->md[iW]->isFilled())
	  {
	    this->ms->add(w->md[iW]->average(),0);
	  }
	w->md[iW]->reset();
	nC=0;
      }
    
    
    
  }
};

template<class walker_t,class wave_t>
class center_of_mass_differenceSquared_future_walker : public measurement_scalar<walker_t,wave_t>
{
public:
  center_of_mass_differenceSquared_future_walker(measure_scalar* m_,int setA_,int setB_,int iW_) : measurement_scalar<walker_t,wave_t>::measurement_scalar(m_){iW=iW_;setA=setA_;setB=setB_;accumulator=0;nMeasures=0;}

  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    
    // add to the walker the measurement just done
    if (w->md[iW]->isFilled())
      {
        accumulator+=w->md[iW]->average();
	nMeasures+=1;
      }
    
    w->md[iW]->add(pow( wave->center_of_mass_no_pbc(w->state,setA) - wave->center_of_mass_no_pbc(w->state,setB),2));
    
  }
  
private:
  int iW;
  int setA;
  int setB;
  double accumulator;
  int nMeasures;
  virtual void record(double step)
  {
    if (nMeasures>0)
      {
	this->ms->add(accumulator/nMeasures,step);
      }
    nMeasures=0;
    accumulator=0;
  }
  
};

// performs a time difference type of measurement
template <class walker_t,class wave_t>
class winding_number_time : public measurement<walker_t,wave_t,measure_vector_mult_index> 
{
  bool filled;
  int iWm;
public:
  winding_number_time(int iWm_,measure_vector_mult_index * ms_,int set_a_) : measurement<walker_t,wave_t,measure_vector_mult_index>::measurement(ms_){iWm=iWm_;filled=false;this->ms->set_history(true);set_a=set_a_;};
  virtual void out()
  {
    if (filled)
      {	
	measurement<walker_t,wave_t,measure_vector_mult_index>::out();
	this->clear();
	filled=false;
      }
    
  }
  
  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    w->md[iWm]->time_difference_average(w->state,wave,this->ms->sum,this->ms->ns,wave->center_of_mass_no_pbc(w->state,set_a));
    
    filled=w->md[iWm]->isFilled();
  }
  virtual void clear()
  {
    //cout << "clear"<<endl;
  };
private:
  int set_a;
};

template <class walker_t,class wave_t>
class winding_number_time_sum : public measurement<walker_t,wave_t,measure_vector_mult_index> 
{
  bool filled;
  int iWm;
public:
  winding_number_time_sum(int iWm_,measure_vector_mult_index * ms_,int set_a_,int set_b_) : measurement<walker_t,wave_t,measure_vector_mult_index>::measurement(ms_){iWm=iWm_;filled=false;this->ms->set_history(true);set_a=set_a_;set_b=set_b_;};
  
  virtual void out()
  {
    if (filled)
      {	
	measurement<walker_t,wave_t,measure_vector_mult_index>::out();
	this->clear();
	filled=false;
      }
    
  }
  
  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    w->md[iWm]->time_difference_average(w->state,wave,this->ms->sum,this->ms->ns,wave->center_of_mass_no_pbc(w->state,set_a) + wave->center_of_mass_no_pbc(w->state,set_b));
    
    filled=w->md[iWm]->isFilled();
  }
  virtual void clear()
  {
    //cout << "clear"<<endl;
  };
private:
  int set_a;
  int set_b;
};

template <class walker_t,class wave_t>
class winding_number_time_diff : public measurement<walker_t,wave_t,measure_vector_mult_index> 
{
  bool filled;
  int iWm;
public:
  
  winding_number_time_diff(int iWm_,measure_vector_mult_index * ms_,int set_a_,int set_b_) : measurement<walker_t,wave_t,measure_vector_mult_index>::measurement(ms_){iWm=iWm_;filled=false;this->ms->set_history(true);set_a=set_a_;set_b=set_b_;};
  
  virtual void out()
  {
    if (filled)
      {	
	measurement<walker_t,wave_t,measure_vector_mult_index>::out();
	this->clear();
	filled=false;
      }
    
  }
  
  virtual void make_measurement(walker_t* w,wave_t* wave)
  {
    w->md[iWm]->time_difference_average(w->state,wave,this->ms->sum,this->ms->ns,wave->center_of_mass_no_pbc(w->state,set_a) - wave->center_of_mass_no_pbc(w->state,set_b));
    
    filled=w->md[iWm]->isFilled();
  }
  
  virtual void clear()
  {
    //cout << "clear"<<endl;
  };
private:
  int set_a;
  int set_b;
};

template<class walker_t,class wave_t,class qmcKind>
class  winding_number_creator
{
  typedef measurementInterface<walker_t,wave_t> measurementInterface_t;
public:
  bool isSupported(){return false;};
  void append(vector< measurementInterface_t * > &ms,int id,xml_input* main_input,string label)
  {
    
  };
};

template<class tm>
class measures
{
public:
  typedef typename tm::qmcKind qmcKind;
  typedef typename tm::walker_t walker_t;
  typedef typename tm::wave_t wave_t;
  typedef measurementInterface<walker_t,wave_t> measurementInterface_t;
  typedef tm qmc_t;
  typedef walker_t measure_obj_t;
  vector< measurementInterface_t * > ms;
  
  //void average_with(measures* m);
  void out();
  void save();
  void increment();
  void clear();
  void reduce();
  void record(double step_);
  void make_measurements(measure_obj_t* w,wave_t* wave);
  measures(string filename,qmc_t* qmc_obj);
  
  winding_number_creator<walker_t,wave_t,qmcKind> winding_number_creator_obj;
  
};


// center_of_mass_w<double>* build_center_of_mass_w(xml_input* xml_m);

// template <class md_t>
// void build_md(vector<md_t*> &md)
// {
//   xml_input* main_input;
//   main_input=new xml_input;
//   main_input->open("input.xml");
//   main_input->reset()->get_child("measures")->get_first_child();
  
//   md.resize(0);
//   // while checking on other things
//   while(main_input->check())
//     {
//       if (main_input->get_name() == "winding_number")
// 	{
// 	  md.push_back(build_center_of_mass_w(main_input));
// 	}
//       main_input->get_next();    
//     }
// }

template<class walker_t,class wave_t>
class magnetizationMeasurement : public measurement_scalar<walker_t,wave_t>
{
public:
  magnetizationMeasurement(measure_scalar* mScal) :  measurement_scalar<walker_t,wave_t>::measurement_scalar(mScal){};
  
  void make_measurement(walker_t* w,wave_t* wave){this->ms->add(getMagnetization((*w->state)[0]),0); }
};

measure_scalar* build_measure_scalar(xml_input* xml_m,string label);

measure_vector_mult_index* build_measure_vector_mult_index(xml_input* xml_m,string label);


measure_vector* build_measure_vector(xml_input* xml_m,string label);


#include "measures.hpp"
#include "observables/optimizationObservablesLinearMethod.h"
#endif
