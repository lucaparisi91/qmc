#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <string>
#include "xml-input.h"
#include "jastrow.h"
#include "system.h"
#include "potential.h"
#include "traits.h"
#include "tools.h"
#include "qmc.h"
#include "exception_classes.h"
#include "measures.h"
#include "particles/gradientParticles.h"

// class of the wavefuntion
class space_vector_measure;
class deviation_measure_vector_time;
class measure_scalar;
class overlap_measure;
class measure_energy_difference;


// the wavefunction class previously studied
template<class tm>
class wavefunction
{
public:
  typedef tm qmc_t;
  typedef typename qmc_t::particles_t particles_t;
  typedef typename qmc_t::all_particles_t all_particles_t;
  typedef typename qmc_t::potential_t potential_t;
  typedef typename qmc_t::value_t value_t;
  typedef allParticlesGradient1D grad_t;
  
  qmc_t* qmc_obj;
  int src_particle_set;
  int target_particle_set;
  
  //void clean_state();
  // tmp variables for measurements
  value_t energy_first_term;
  double log_wavefunction_value;
  double wavefunction_derivative;
  int sign; // used in fixed diffusion monte carlo
  double phase; // used in fixed phase diffusion monte carlo
  string label;
  
  wavefunction(qmc_t * qmc_obj);
  void getLabel(){return label;};
  void setLabel(string label_){label=label_;};
  virtual double evaluate_derivative(all_particles_t *p){throw notYetSsupported("evaluate_derivative");};
  virtual double evaluate_derivative_second(all_particles_t *p,double phi){throw notYetSsupported("evaluate_derivative_second");};
  virtual wavefunction<tm>* clone(){throw notYetSsupported("Clone wavefunction " + label);};
  
  void setOptParameter(int i){optParameter=i;}
  int getOptParamater(){return optParameter;}
  virtual double log_evaluate(all_particles_t *p){};//evaluates the logarithm of the jastrow
  
  virtual void print(int i)=0; // prints information about the wavefunction to file (jastrows)
  virtual double potential(empty_t *,all_particles_t* state)=0;
  virtual double potential(rabiCoupling *,all_particles_t* state)=0;
  
  virtual double one_particle_log_evaluate(all_particles_t * p,double r,int iP, int set_a){throw notYetSsupported("one_particle_log_evaluate");};
  
  /* 
     Optimization
  */
  virtual void setParameter(double x,int i){throw notYetSsupported("setParameter");};
  
  virtual double getParameter(int i){throw notYetSsupported("getParameter");};

  double getParameter(){return getParameter(optParameter);}
  double setParameter(double x){setParameter(x,optParameter);}
  
  // mark the wavefunction for later optimization
  void setOptimized(bool toOptimize_){toOptimize=toOptimize_;};
  //virtual double derivative1Param(all_particles_t* p){return 0;};
  //virtual double derivative2Param(all_particles_t* p){return 0;};
  bool isOptimized(){return toOptimize;}

  /*------------------------------------------------------------
    Wavefunction observables
    -----------------------------------------------------------*/
  
  // compute the gradient of the wavefunction
  // Just adds to the gradient. Does not reset the initial value
  virtual void gradient(const all_particles_t & p,grad_t & grad){throw notYetSsupported("gradient");};
  // compute the laplacian of the wavefunction
  virtual value_t laplacian(all_particles_t & p,grad_t & grad){throw notYetSsupported("laplacian");};
  // evaluates the laplacian, the gradient and the logarithm of the value of the wavefunction

  /*
    - evaluates \grad^2 \psi - (\grad psi)^2
    - evaluates \grad \psi 
*/
  
  virtual void laplacianMinusGradientSquared(const all_particles_t & p,grad_t & grad,value_t & e){throw notYetSsupported("logEvaluate");};
  
  // evaluate the value of the wavefunction
  virtual value_t logEvaluate(all_particles_t &p){throw notYetSsupported("logEvaluate");};
protected:
  wavefunction(){};
private:
  
  bool toOptimize;
  int optParameter;
  
};


// defines a bill jastrow form wavefunction
template<class jastrow_t,class tm>
class bill_jastrow_wavefunction : public wavefunction< tm >
{
public:
  typedef typename wavefunction<tm>::particles_t particles_t;
  typedef typename wavefunction<tm>::qmc_t qmc_t;
  typedef typename wavefunction<tm>::potential_t potential_t;
  typedef typename wavefunction<tm>::all_particles_t all_particles_t;
  typedef typename wavefunction<tm>::value_t value_t;
  virtual double getParameter(int i){throw notYetSsupported("get_parameter bill_jastrow_wavefunction " + this->label );};
  virtual void setParameter(double x,int i){throw notYetSsupported("setParameter bill_jastrow_wavefunction");};
  bill_jastrow_wavefunction(qmc_t* qmc_obj_,xml_input* xml_wave,string filename);

  
  
  void overlap(const double &w1,const double &w2,measure_scalar* m);
  virtual void print(int i); 
  virtual double log_evaluate(all_particles_t *p)=0;
  
  
  virtual double potential(empty_t *p,all_particles_t* state) {return 0;};
  
  virtual double potential(rabiCoupling *ps,all_particles_t* state){throw notYetSsupported("potential bill_jastrow_wavefunction");};
  
protected:
  
  jastrow_t jastrowc;
  bill_jastrow_wavefunction(jastrow_t &jastrowo) : jastrowc(jastrowo) {} ; 
  void copyTo(bill_jastrow_wavefunction<jastrow_t,tm> * wave2)
  {
    
    wave2->qmc_obj=this->qmc_obj;
    wave2->src_particle_set=this->src_particle_set;
    wave2->target_particle_set=this->target_particle_set;
    wave2->label=this->label;
    
    
  }
};

template<typename  jastrow_t,class comp>
class bill_jastrow_wavefunction_two_body_symmetric : public bill_jastrow_wavefunction< jastrow_t,comp>
{
  
  typedef typename wavefunction<comp>::qmc_t qmc_t;
  typedef typename wavefunction<comp>::particles_t particles_t;
  typedef typename wavefunction<comp>::all_particles_t all_particles_t;
  typedef typename wavefunction<comp>::value_t value_t;
  typedef typename wavefunction<comp>::grad_t grad_t;
  
 public:
  
  bill_jastrow_wavefunction_two_body_symmetric(qmc_t* qmc_obj_,xml_input* xml_wave,string filename) : bill_jastrow_wavefunction< jastrow_t,comp>(qmc_obj_,xml_wave,filename) {}
  double log_evaluate(all_particles_t *p);
  virtual double one_particle_log_evaluate(all_particles_t * p,double r,int iP, int set_a);
  
  
  virtual bill_jastrow_wavefunction_two_body_symmetric<jastrow_t,comp> * clone()
  {
    bill_jastrow_wavefunction_two_body_symmetric<jastrow_t,comp> * wave2;

    wave2=new bill_jastrow_wavefunction_two_body_symmetric(this->jastrowc);
    this->copyTo(wave2);
    return wave2;
  }
  
  virtual double getParameter(int i){return this->jastrowc.getParameter(i);};
  
  virtual void setParameter(double x,int i){this->jastrowc.setParameter(x,i);}
  
  void laplacianMinusGradientSquared(const all_particles_t & p, grad_t & grad,value_t & e);
  
  virtual void gradient(const all_particles_t & p,grad_t & grad);
  
protected:
  
  bill_jastrow_wavefunction_two_body_symmetric(jastrow_t & jastrowo) : bill_jastrow_wavefunction< jastrow_t,comp>(jastrowo)
  {}    
  
};

template<typename  jastrow_t,typename comp>
class bill_jastrow_wavefunction_one_body : public bill_jastrow_wavefunction< jastrow_t,comp>
{
 public:
  typedef typename wavefunction<comp>::qmc_t qmc_t;
  typedef typename wavefunction<comp>::particles_t particles_t;
  typedef typename wavefunction<comp>::all_particles_t all_particles_t;
  typedef typename wavefunction<comp>::grad_t grad_t;
  typedef typename wavefunction<comp>::value_t value_t;
  
  bill_jastrow_wavefunction_one_body(qmc_t* qmc_obj_,xml_input* xml_wave,string filename) : bill_jastrow_wavefunction< jastrow_t,comp>(qmc_obj_,xml_wave,filename) {}
  double log_evaluate(all_particles_t *p);
  
  virtual double one_particle_log_evaluate(all_particles_t * p,double r,int iP, int set_a);
  // saves the unormalized drift force derivatives
  // void set_drift_force_derivative(all_particles_t* p,vector<double> &drift_force_derivative );
  // compute the first term of the kinetic energy(not normalized)
  
  void kinetic_energy_derivative_first_term(all_particles_t *p , double &,const double wavefunctionValue);
  // returns the value of the first derivative
  virtual double evaluate_derivative(all_particles_t *p);
  virtual double evaluate_derivative_second(all_particles_t *p,double phi);
  virtual void setParameter(double x,int i){this->jastrowc.setParameter(x,i);}
  
  virtual double getParameter(int i){return this->jastrowc.getParameter(i);};
  
  virtual void laplacianMinusGradientSquared(const all_particles_t & p, grad_t & grad,value_t & e);

  void gradient(const all_particles_t & p, grad_t & grad);
  virtual  bill_jastrow_wavefunction_one_body<jastrow_t,comp> * clone()
  {
    
     bill_jastrow_wavefunction_one_body<jastrow_t,comp> * wave2;
     
     wave2=new  bill_jastrow_wavefunction_one_body(this->jastrowc);
     
     this->copyTo(wave2);
     return wave2;
  }

protected:
   
  bill_jastrow_wavefunction_one_body(jastrow_t & jastrowo) : bill_jastrow_wavefunction< jastrow_t,comp>(jastrowo)
  {}    
};

template<typename  jastrow_t,typename comp>
class bill_jastrow_wavefunction_two_body_asymmetric : public bill_jastrow_wavefunction< jastrow_t,comp>
{
  
  typedef typename wavefunction<comp>::particles_t particles_t;
  typedef typename wavefunction<comp>::qmc_t qmc_t;
  typedef typename wavefunction<comp>::all_particles_t all_particles_t;
  typedef typename wavefunction<comp>::value_t value_t;
  typedef typename wavefunction<comp>::grad_t grad_t;
 public:
  bill_jastrow_wavefunction_two_body_asymmetric(qmc_t* qmc_obj_,xml_input* xml_wave,string filename) : bill_jastrow_wavefunction< jastrow_t,comp>(qmc_obj_,xml_wave,filename) {}
  double log_evaluate(all_particles_t *p);
  
  //virtual double derivative1Param(all_particles_t* p);
  //virtual double derivative2Param(all_particles_t* p);
  virtual double one_particle_log_evaluate(all_particles_t * p,double r,int iP, int set_a);
  
  virtual double evaluate_derivative(all_particles_t* p);
  virtual double evaluate_derivative_second(all_particles_t* p,const double phi);
  virtual double getParameter(int i){return this->jastrowc.getParameter(i);};
  virtual void setParameter(double x,int i){this->jastrowc.setParameter(x,i);}
  
  virtual void laplacianMinusGradientSquared(const all_particles_t & p, grad_t & grad,value_t & e);
  
  void gradient(const all_particles_t & p, grad_t & grad);
  
  virtual bill_jastrow_wavefunction_two_body_asymmetric<jastrow_t,comp> * clone()
  {
    
    bill_jastrow_wavefunction_two_body_asymmetric<jastrow_t,comp> * wave2;
    
    wave2=new bill_jastrow_wavefunction_two_body_asymmetric(this->jastrowc);
    this->copyTo(wave2);
    return wave2;
    
  }
protected:
  
  bill_jastrow_wavefunction_two_body_asymmetric(jastrow_t & jastrowo) : bill_jastrow_wavefunction< jastrow_t,comp>(jastrowo)
  {} 
};

template<typename  jastrow_t,class comp>
class bill_jastrow_wavefunction_spinor_two_body_symmetric : public wavefunction< comp >
{
  typedef typename wavefunction<comp>::qmc_t qmc_t;
  typedef typename wavefunction<comp>::all_particles_t all_particles_t;
  typedef typename wavefunction<comp>::particles_t particles_t;
  typedef typename wavefunction<comp>::potential_t potential_t;
  typedef typename wavefunction<comp>::value_t value_t;
  
 public:
  
  bill_jastrow_wavefunction_spinor_two_body_symmetric(qmc_t* qmc_obj_,string filenameUp,string filenameUpDown,string fileParams) : wavefunction<comp>(qmc_obj_),jastrowc(filenameUp,filenameUpDown,fileParams) {}
  
  double log_evaluate(all_particles_t *p);
  virtual double potential(rabiCoupling *p,all_particles_t* state);
  virtual double potential(empty_t *e,all_particles_t* state){ return 0;};
  void print(int i)
  {
    
    jastrowc.print(string("jastrow")+ int_to_string(i)+ string(".dat"),string("jastrow_1d")+ int_to_string(i) + string(".dat"),string("jastrow_2d")+ int_to_string(i) + ".dat");
    
  };
protected:
  jastrow_t jastrowc;
};

// total wavefunction of the system
template<class tm>
// total wavefuncton: made up of all other wavefunctions

class total_wavefunction
{
public:
  typedef tm qmc_t;
  typedef typename qmc_t::all_particles_t all_particles_t;
  typedef typename qmc_t::value_t value_t;
  typedef typename qmc_t::particles_t particles_t;
  typedef allParticlesGradient1D grad_t;
  
  string label;
  qmc_t* qmc_obj;
  vector<wavefunction<qmc_t>*> waves;
  
  total_wavefunction(const total_wavefunction<tm> *waveTot2)
  {
    label=waveTot2->label;
    qmc_obj=waveTot2->qmc_obj;
    waves.resize(waveTot2->waves.size());
    
    //deep copy of individual wavefunctions
    for(int i=0;i<waves.size();i++)
      {
	waves[i]=waveTot2->waves[i]->clone();
      }
    e=waveTot2->e;
    wavesOptim=waveTot2->wavesOptim;
    phisWork=waveTot2->phisWork;
    phis2Work=waveTot2->phis2Work;
  }
  
  // returns the total energy and updates the drift force of the status

  /* 
    - computes the laplacian and the gradient of the total wavefunction
    - the gradient is then overwritten
   */
  
  void laplacianGradient(const all_particles_t &p,value_t &e,value_t & eF,grad_t & grad);

  /* 
     - Computes just the gradient
     - The gradient is reset and overwritten
   */
  virtual void gradient(const all_particles_t & p,grad_t & grad);
  
  template<class mV_t>
  void pair_correlation_asymm(all_particles_t *p, mV_t* g,double max,int set_a,int set_b);

  
  
  template<class mV_t>
  void pair_correlation_symm(all_particles_t *p, mV_t* g,double max,int set_a);
  
  void center_of_mass_no_pbc(all_particles_t *p,measure_scalar* g);
  // asymmetric structure factor
  void structure_factor_asymm(all_particles_t *p, measure_vector* s,vector<double> &qs,int setA,int setB);
  // symmetric structure factor
  void structure_factor_symm(all_particles_t *p, measure_vector* s,vector<double> &qs,int setA);
  template<class measure_vec_t>
  void structure_factor(all_particles_t *p, measure_vec_t* s,vector<double> &qs,int setA,vector<complex< double> > & work);
  template<class measure_vec_t>
  void structure_factor(all_particles_t *p, measure_vec_t* s,vector<double> &qs,int setA,int setB,vector<complex< double> > & work);
  
  template<class measure_vec_t>
  void structure_factorSpin(all_particles_t *p, measure_vec_t* s,vector<double> &qs,int setA,int setB,vector<complex< double> > & work);
  
  void density(all_particles_t* p,space_measure<measure_vector>* density_m);
  // used to compute the overlap with some other wavefunctions
  //void overlap(all_particles* p, overlap_measure* m);
  // allows to cause an overlap with some other means of computations
  
  void oneBodyDensityMatrixOffdiagonal(all_particles_t* p,space_measure<measure_vector_mult_index>* g,int nMCM,int Max);
  
  void oneBodyDensityMatrixOffdiagonalFutureWalkers(all_particles_t* p,space_measure<measure_vector_mult_index>* g,int nMCM,int Max,miniDMCDriver< qmc_t >* driver);
  
  int getWaveid(string label)
  {
    int j;
    for(j=0;j<waves.size();j++)
      {
	if (label==waves[j]->getLabel())
	  {
	    return j;
	  }
      }
  }
  
  double log_evaluate(all_particles_t *p);
  // evaluates the derivative for a certain wavefunction
  
  double evaluate_derivative(all_particles_t*p,vector<int> & wavesOptim)
  {
    int i,j;
    double w;
    w=0;
    
    for(i=0;i<waves.size();i++)
      {
	if (waves[j]->isOptimized()==true)
	  {
	    w+=waves[j]->evaluate_derivative(p);
	  }
      }
    return w;
    
  }

  int findWavesToOptimize()
  {
    int i;
    wavesOptim.resize(0);
    
    for(i=0;i<waves.size();i++)
      {
	if (waves[i]->isOptimized()==true)
	  {
	    wavesOptim.push_back(i);
	  }
      } 
  }
  
  // computes the first parameter derivative of the wavefunction
  void computeDerivatesFirst(vector<double> & phis,all_particles_t* p,double &phi_T)
  {
    int j,i;
    phi_T=0;
    phis.resize(wavesOptim.size());
    
    for(j=0;j<wavesOptim.size();j++)
    {
      i=wavesOptim[j];	 
      phis[j]=waves[i]->evaluate_derivative(p);
      phi_T+=phis[j];
    }
    
  
  }
  // computes first and second derivatives of the wavefunction
  
  void computeDerivatesFirstSecond(vector<double> & phis,vector<double> & phis2,all_particles_t* p,double &phi_T,double &phi2_T)
  {
    int i,j;

    phis.resize(wavesOptim.size());
    phis2.resize(wavesOptim.size());
    // first compute the first derivatives
    computeDerivatesFirst(phis,p,phi_T);
    
    
    for(j=0;j<wavesOptim.size();j++)
    {
      i=wavesOptim[j];	 
      phis2[j]=waves[i]->evaluate_derivative_second(p,phis[j]);
      
    }
    
    // compute the total second derivative with respect to the wavefunction
    
    phi2_T=0;
    for(i=0;i<wavesOptim.size();i++)
      {
	// second derivatives
	phi2_T+=phis2[i];
	// mixed product of first derivatives
	for(j=0;j<wavesOptim.size();j++)
	  {
	    if (j==i)
	      {
		continue;
	      }
	    
	    phi2_T+=phis[i]*phis[j];
	  }
	
      }
  }
  
  void setParameter(double x,int i,int j)
  {
    
    waves[j]->setParameter(x,i);
    
  }
  
  void setParameter(double x)
  {
    
    for(int i=0;i<wavesOptim.size();i++)
      {
	waves[wavesOptim[i]]->setParameter(x);
      }
  }
  
  double getParameter(int i,int j) const
  {
    return waves[j]->getParameter(i);
  }
  
  double getParameter()
  {
    return waves[wavesOptim[0]]->getParameter();
  }
  
  
  // computes the total wavefunction
  total_wavefunction(qmc_t* qmc_obj_);
  void print_jastrows();
  void set_drift_force(int set,vector<double> &drift_force);
  
  template<int>
    void energy_difference(all_particles_t* p,measure_energy_difference* med);

  void link_wavefunctions(xml_input* xml_wave,vector<wavefunction<qmc_t>* > wave_table,string label);
  // returns the center of mass without pbc of a certain set
  double center_of_mass_no_pbc(all_particles_t* p,int set);
  void set_drift_force(all_particles_t* p);
  void potential(empty_t * p,all_particles_t* state,double &ev){ev=0;};
  // void potential(rabiCoupling* p,all_particles_t* state,double &ev)
  // {
  //   int i=0;
  //   ev=0;
  //   for(i=0;i<waves.size();i++)
  //     {
  // 	ev+=waves[i]->potential(p,state);
  //     }
  //   ev=ev*p->get_omega();
  // }
  
  int get_sign();
  double get_phase();
  double e;
private:
  vector<double> wavesOptim;
  vector<double> phisWork;
  vector<double> phis2Work;
};

template<class J,class comp>
typename comp::swave_t* createBillJastrowTwoBodySymmetric(comp* qmc_obj,xml_input* xml_wave ,const string &filename)
  {
    
    return new bill_jastrow_wavefunction_two_body_symmetric<J,comp>(qmc_obj,xml_wave,filename);
  };

template<class J,class comp>
typename comp::swave_t* createBillJastrowTwoBodyAsymmetric(comp* qmc_obj,xml_input* xml_wave ,const string &filename)
  {
    return new bill_jastrow_wavefunction_two_body_asymmetric<J,comp>(qmc_obj,xml_wave,filename);
  };

template<class J,class comp>
typename comp::swave_t* createBillJastrowOneBody(comp* qmc_obj,xml_input* xml_wave ,const string &filename)
  {
    return new bill_jastrow_wavefunction_one_body<J,comp>(qmc_obj,xml_wave,filename);
  };

template<class J,class comp>
typename comp::swave_t* createBillJastrowSpinorTwoBodySymmetric(comp* qmc_obj,xml_input* xml_wave ,const string &filename1,const string &filename2)
{
  return new bill_jastrow_wavefunction_spinor_two_body_symmetric<J,comp>(qmc_obj,xml_wave,filename1,filename2);
};


string createWavefunctionId(xml_input* xml_wave);
string createJastrowId(xml_input* xml_wave);

template<class comp>
void load_wavefunctions(xml_input * xml_wave, vector< typename comp::swave_t* > &waves,comp* qmc_obj);

#include "wavefunction.hpp"
#include "measure_functions.hpp"

#include "wavefunction/billJastrowWaveFunctionOneBody.hpp"
#include "wavefunction/billJastrowWaveFunctionTwoBodySymm.hpp"
#include "wavefunction/billJastrowWaveFunctionTwoBodyASymm.hpp"
#include "wavefunction/totalWavefunction.hpp"

#endif
