#ifndef MEASURESPIGS_H
#define MEASURESPIGS_H

#include "../measures.h"

template<class walker_t,class wave_t>
class energyTail : public measurement_scalar<walker_t,wave_t>
{
public:
  typedef typename wave_t::grad_t grad_t;
  
  energyTail(measure_scalar* mScal,string filename) :  measurement_scalar<walker_t,wave_t>::measurement_scalar(mScal){init(filename);};
  
  void init(string filename)
  {
    gradParticles=*buildAllParticlesGradient1D(filename);
  }
  
  void make_measurement(walker_t* w,wave_t* wave)
  {
    double e, dummy;
    wave->laplacianGradient((*w)[0],e,dummy,gradParticles);
    this->ms->add(e,0);
  }
  
private:
  grad_t gradParticles;
};



template<class measures_t>
void buildPIGSMeasures(string filename,measures_t & ms)
{
  xml_input* main_input;
  string label;
  int i=0;
  main_input=new xml_input;
  main_input->open("input.xml");
  main_input->reset()->get_child("measures")->get_first_child();
  while( main_input->check() )
    {
      i=i+1;
      
      if (main_input->get_name() == "energyTail")
	{

	  if (main_input->get_attribute("label") != NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="energyTail";
	    }
	  
	  ms.push_back(new energyTail<typename measures_t::walker_t,typename measures_t::wave_t >(build_measure_scalar(main_input,label),"input.xml"));
	}
      
      main_input->get_next();
    }
  
  delete main_input;
}

#endif
