template<class comp>
measurementInterface<typename comp::walker_t,typename comp::wave_t>*  build_structure_factor(xml_input* main_input,double l_box, unsigned int & id)
{
  int bins;
  int set_a;
  int set_b;
  double deltaQ;
  string label;
  bool complex,spin;
  bool futureWalker;
  double max;
  bool futureWalkers;
  int nFutureWalkers;
  
  typedef typename comp::wave_t wave_t;
  typedef typename comp::walker_t walker_t;
  typedef typename comp::qmcKind qmcKind;
  
  structure_future_single_creator<walker_t,wave_t,qmcKind> structureFutureObj;
  
  if (main_input->get_name() == "static_structure_factor")
    {
      
      if (main_input->get_attribute("label")!=NULL)
	{
	  label=main_input->get_string();
	}
      else
	{
	  label="static_structure_factor";
	}
      
      if (main_input->get_attribute("futureWalkers")!=NULL)
	{
	  futureWalkers=main_input->get_bool();
	}
      else
	{
	  futureWalkers=false;
	}
      
      
      if (main_input->get_attribute("setA") != NULL)
	{
	  set_a=main_input->get_int();
	}
      else
	{
	  set_a=0;
	}

      
      if (main_input->get_attribute("nFutureWalkers") != NULL)
	{
	  nFutureWalkers=main_input->get_int();
	}
      else
	{
	  nFutureWalkers=100;
	}
      

      if (main_input->get_attribute("setB") != NULL)
	{
	  set_b=main_input->get_int();
	}
      else
	{
	  set_b=0;
	}
	  
      if (main_input->get_attribute("bins") != NULL)
	{
	  bins=main_input->get_int();
	}
      else
	{
	  bins=100;
	}
      
      if (main_input->get_attribute("spin") != NULL)
	{
	  spin=main_input->get_bool();
	}
      else
	{
	  spin=false;
	}
	  
      if (main_input->get_attribute("max")!= NULL)
	{
	  max=main_input->get_real();
	}
      else
	{
	  max=bins*2*M_PI/l_box;
	}
      
      deltaQ=getDeltaQ(l_box,max,bins);
      
      
      
      if (set_a==set_b)
	{
	  
	  if (spin==false)
	    {
	      
	      if (futureWalkers==true)
		{
		   return new structure_factor_density_complexOrbitalsForwardWalking<walker_t,wave_t>(
										    build_measure_vector(main_input,label),bins,deltaQ,l_box,set_a,set_b,id);
		   
		}
	      else
		{
		  return 
		    new structure_factor_single_complex<walker_t,wave_t>(
									 build_measure_vector(main_input,label),bins,deltaQ,l_box,set_a
									 );
		}
	    }
	
	
	  else
	    {
	      
	      if (futureWalkers==true)
		{
		  
		  return new structure_factor_spin_complexOrbitalsForwardWalking<walker_t,wave_t>(
										    build_measure_vector(main_input,label),bins,deltaQ,l_box,set_a,set_b,id);
		}
	      else
		{
		  return new structure_factor_spin_complexOrbitals<walker_t,wave_t>(
							       build_measure_vector(main_input,label),bins,deltaQ,l_box,set_a,set_b							       );
		}
	    }
	}
      else
	{
	  if (futureWalkers==true)
	    {
	      cout << "Future walkers not yet supported";
	      exit(1);
	    }
	    
	  if (spin==false)
	    {
	      return 
		new structure_factor_double_complex<walker_t,wave_t>(
								     build_measure_vector(main_input,label),bins,deltaQ,l_box,set_a,set_b
								     );
	    }
      else
	{
	  return 
	    new structure_factor_spin_complex<walker_t,wave_t>(
							       build_measure_vector(main_input,label),bins,deltaQ,l_box,set_a,set_b
							       );
	  
	}
	  
	}
    }
}

template<class comp>
measure_dynamic<comp> * buildStructureFactorWalker(xml_input* main_input,double l_box)
{
  int set_a;
  int set_b;
  int bins;
  double max;
  bool spin;

  
  
  if (main_input->get_attribute("setA") != NULL)
    {
      set_a=main_input->get_int();
    }
  else
    {
      set_a=0;
    }
  if (main_input->get_attribute("setB") != NULL)
    {
      set_b=main_input->get_int();
    }
  else
    {
      set_b=set_a;
    }
  
  if (main_input->get_attribute("bins") != NULL)
    {
      bins=main_input->get_int();
    }
  else
    {
      bins=100;
    }
  if (main_input->get_attribute("max")!= NULL)
    {
      max=main_input->get_real();
    }
  else
    {
      max=bins*2*M_PI/l_box;
    }

  
  if (main_input->get_attribute("spin")!= NULL)
    {
      spin=main_input->get_bool();
    }
  else
    {
      spin=false;
    }
  
  if(set_a==set_b)
    {
      return new structureFactorWalker<comp,double>(build_q_vector(bins,l_box,max),set_a);
    }
  else
    {
      if (spin==false)
	{
	  
	  return new structureFactorDoubleWalker<comp,double>(build_q_vector(bins,l_box,max),set_a,set_b);
	  
	}
      else
	{
	  
	  return new structureFactorSpinWalker<comp,double>(build_q_vector(bins,l_box,max),set_a,set_b);
	}
    }
  
  
}


template<class comp>
measure_dynamic<comp> * buildPairCorrelationFuture(xml_input* main_input,double l_box)
{
  int set_a;
  int set_b;
  int bins;
  double max;
  bool spin;

  
  if (main_input->get_attribute("setA") != NULL)
    {
      set_a=main_input->get_int();
    }
  else
    {
      set_a=0;
    }
  if (main_input->get_attribute("setB") != NULL)
    {
      set_b=main_input->get_int();
    }
  else
    {
      set_b=set_a;
    }
  
  if (main_input->get_attribute("bins") != NULL)
    {
      bins=main_input->get_int();
    }
  else
    {
      bins=100;
    }
  if (main_input->get_attribute("max")!= NULL)
    {
      max=main_input->get_real();
    }
  else
    {
      max=bins*2*M_PI/l_box;
    }
  mesh grid(0,max,bins);
  
  if (set_a == set_b)
    {
      return new pairCorrelationSymmFuture<comp,double,mesh>(grid,set_a,bins,max);
    }
  else
    {
      return new pairCorrelationAsymmFuture<comp,double,mesh>(grid,set_a,set_b,bins,max);
    }
}

template<class comp>
measurementInterface<typename comp::walker_t,typename comp::wave_t>* build_pair_correlation(xml_input* main_input,double lBox,unsigned int & id)
{
  double max;
  string label;
  int set_a,set_b,bins,nFutureWalkers;
  bool futureWalkers;

  typedef typename comp::walker_t walker_t;
  typedef typename comp::wave_t wave_t;
  typedef typename comp::qmcKind qmcKind;

  structure_future_single_creator<walker_t,wave_t,qmcKind> structureFutureObj;

  
  if (main_input->get_attribute("max") != NULL)
            {
	      max=main_input->get_real();
	    }
	  else
	    {
	      max=lBox/2;
	    }
	  if (main_input->get_attribute("label")!=NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="pair_correlation";
	    }

	  if (main_input->get_attribute("setA") != NULL)
	    {
	      set_a=main_input->get_int();
	    }
	  else
	    {
	      set_a=0;
	    }
	  
	  if (main_input->get_attribute("setB") != NULL)
	    {
	      set_b=main_input->get_int();
	    }
	  else
	    {
	      set_b=0;
	    }

	  if (main_input->get_attribute("bins") != NULL)
	    {
	      bins=main_input->get_int();
	    }
	  else
	    {
	      bins=100;
	    }

	  if (main_input->get_attribute("nFutureWalkers") != NULL)
	    {
	      nFutureWalkers=main_input->get_int();
	    }
	  else
	    {
	      nFutureWalkers=100;
	    }

	  if (main_input->get_attribute("futureWalkers")!=NULL)
	    {
	      futureWalkers=main_input->get_bool();
	    }
	  else
	    {
	      futureWalkers=false;
	    }

	  if (futureWalkers==true)
	    {
	      id++;
	  
	      return structureFutureObj.create(main_input,label,id-1,nFutureWalkers,bins);
	    }
	  else
	    {
	      return new pair_correlation_m<walker_t,wave_t>(build_space_vector<measure_vector>(main_input,label,0,max),set_a,set_b);
	    }
	  
}

#include "observables/imaginaryTimeCorrelation.h"

template<class walker_t,class wave_t>
measurementInterface<walker_t,wave_t>*  buildWindingNumber(xml_input* main_input, unsigned int id)
{
  int bins;
  int set_a;
  string label;
  
  if (main_input->get_attribute("label")!=NULL)
    {
      label=main_input->get_string();
    }
  else
    {
      label="windingNumber";
    }
  
  if (main_input->get_attribute("setA") != NULL)
    {
      set_a=main_input->get_int();
    }
  else
    {
      set_a=0;
    }

    if (main_input->get_attribute("bins") != NULL)
    {
      bins=main_input->get_int();
    }
  else
    {
      bins=100;
    }
    
  return new windingNumber<walker_t,wave_t>(build_measure_vector(main_input,label),bins,set_a,id);
      
}

template<class walker_t,class wave_t>
measurementInterface<walker_t,wave_t>*  buildWindingNumberSpin(xml_input* main_input, unsigned int id)
{
  int bins;
  int set_a;
  string label;
  
  if (main_input->get_attribute("label")!=NULL)
    {
      label=main_input->get_string();
    }
  else
    {
      label="windingNumberSpin";
    }
  
  if (main_input->get_attribute("setA") != NULL)
    {
      set_a=main_input->get_int();
    }
  else
    {
      set_a=0;
    }

    if (main_input->get_attribute("bins") != NULL)
    {
      bins=main_input->get_int();
    }
  else
    {
      bins=100;
    }
    
  return new windingNumberSpin<walker_t,wave_t>(build_measure_vector(main_input,label),bins,set_a,id);
      
}
