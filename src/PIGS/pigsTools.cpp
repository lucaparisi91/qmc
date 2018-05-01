#include "pigsTraits.h"
#include "pigsTools.h"
#include <sstream>


bool logMetropolis(double log_ratio,system_t::rand_t* rand)
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

string getAttribute(xml_input * main_input,string name,string defaultValue="Unkown")
{
  string label;
  
  if (main_input->get_attribute(name) != NULL)
    {
      label=main_input->get_string();
    }
  else
    {
      label=defaultValue;
    }
  
  return label;
}

string toString(double value)
{
  ostringstream strs;
  strs << value;
  return strs.str();
  
}

towerSampling::towerSampling()
{
  pTot=0;
}
int towerSampling::sample()
{
  double c;
  c=randO->uniform()*pTot;
  
  for(int i=0;i<cumulative.size();i++)
    {
      if (cumulative[i]>=c)
	{
	  return i;
	}
    }
}

void towerSampling::add(double p)
{
  pTot+=p;
  cumulative.push_back(pTot);
}
