#ifndef OPTIMIZATIONOBSERVABLESLINEARMETHOD_H
#define OPTIMIZATIONOBSERVABLESLINEARMETHOD_H

#include<utility>
#include "../measures.h"

template<class walker_t,class wave_t>
class optimizationObservablesLinearMethod : public measurement<walker_t,wave_t,measure_vector >
{
public:
  typedef typename wave_t::grad_t grad_t;
  
  optimizationObservablesLinearMethod(measure_vector* ms_,const wave_t* wave_,double delta_,vector<int> & ns) : measurement<walker_t,wave_t,measure_vector>(ms_)
  {
    wave2=new wave_t(wave_);
    wave=wave_;
    setDelta(delta_);
    optObs.resize(7);
    gradTmp.resize(ns);
    
  };
  
  void make_measurement(walker_t* w,wave_t* wave);
  
  double getEnergy();
  
  vector<double> getMean();
  
  void addParameter(int i,int j);
  
  double getParameter();
  void setParameter(double p,wave_t * waveTmp) const ;
  
  void init();
  
  void setDelta(double delta_)
  {
    delta=delta_;
  }
  
private:
  
  grad_t gradTmp;
  wave_t* wave2;
  const wave_t* wave;
  double delta;
  /*
    stores the optimization plan
   */
  
  vector<double> optObs;
  vector<pair<int,int> > parametersToOptimize;
};


template<class walker_t,class wave_t>
optimizationObservablesLinearMethod<walker_t,wave_t>* buildOptimizationObservablesLinearMethod(xml_input* main_input,wave_t * wave)
{
  vector<int> ns;
  int l1,l2;
  optimizationObservablesLinearMethod<walker_t,wave_t> *objPtr;
  double delta;
  
  // gets the total number of particles(required to know the total number of particles)
  
  main_input->reset()->get_child("system")->get_first_child();
  
  while(main_input->check() )
    {
      if (main_input->get_name()=="particles")
	{
	  main_input->get_attribute("n");
	  ns.push_back(main_input->get_int());
	}
      
      main_input->get_next();
    }
  
  cout << "Optimization plan" <<endl;
  cout << "Ns: "<<endl;
  print_vector(ns);
  
  // get delta
  main_input->reset()->get_child("optimizationPlan");
  delta=main_input->get_attribute("delta")->get_real();
  cout << "delta: "<< delta<<endl;
  
  // create object
  objPtr=new optimizationObservablesLinearMethod<walker_t,wave_t>(new measure_vector(7,"optimizationObs",false),wave,delta,ns);
  
  main_input->get_child("parameter")->get_first_child();
  while(main_input->check() )
    {
      if (main_input->get_name() != "text")
	{
	  main_input->get_attribute("wavefunction");
	  l1=main_input->get_int();
	  main_input->get_attribute("param");
	  l2=main_input->get_int();
	  objPtr->addParameter(l1,l2);	  
	  cout << "Par " << l1 << " " << l2 << endl;
	}
      
      main_input->get_next();
    }
  
  main_input->reset();
  
  //allocate mememory and initialize the auxiliary wavefunction
  objPtr->init();
  
  return objPtr;
}

#endif
