#include "optimizePlan.h"
#include "../xml-input.h"
#include <cassert>
// print the plan to file
void optimizePlan::print()
{
  registerType::iterator it;
  
  for(it = parameterRegister.begin(); it != parameterRegister.end(); it++)
    {
      printf("%s ",it->first.c_str());
      
      for(int i=0;i<it->second.size();i++)
	{
	  printf("(%i,%i) ",it->second[i].first,it->second[i].second);
	  
	}
      
      if (it->second.isBoundedBelow())
	{
	  printf("Min: %f ;",it->second.getMinParameter());
	}
      
      if (it->second.isBoundedAbove())
	{
	  printf("Max: %f ;",it->second.getMaxParameter());
	}
      
      printf("\n");
      
      printf("Delta: %f \n",getDelta());
    }
}


bool optimizePlan::checkBoundsParameters(const vector<double> &parameters) const
{
  registerType::const_iterator it;
  int k=0;
  bool pass=true;
  
  for(it = parameterRegister.begin(); it != parameterRegister.end(); it++)
    {
      assert(k<parameters.size());
      if (
	      
	    ( it->second.isBoundedBelow() ) and
	    ( parameters[k] < it->second.getMinParameter()  )
	  )
	{
	  pass=false;
	}

      if (
	  (it->second.isBoundedAbove()) and
	    (parameters[k] > it->second.getMaxParameter() )
	  )
	{
	  pass=false;
	}

      k++;
    }
      
  
  return pass;
}



void buildPlanOneParameter(xml_input &inputFile,optimizePlan & plan,string key)
{
  xmlNodePtr cur;
  
  int paramKey;
  int waveKey;
  
  cur=inputFile.cur;
  inputFile.get_child("p");
	  
  while( inputFile.check() )
    {
      //printf("name: %s\n",inputFile.get_name().c_str());
	  if(inputFile.get_name()=="p")
	    {
	      
	      if (inputFile.get_attribute("param")!=NULL)
		{
		  paramKey=inputFile.get_int();
		  
		}
	      
	      else
		{
		  printf("Param Key error\n");
		  exit(0);
		}


	      if (inputFile.get_attribute("wavefunction")!=NULL)
		{
		  waveKey=inputFile.get_int();
		}
	      
	      else
		{
		  printf("Wave key error\n");
		  exit(0);
		}
		    



		    
	      plan.add(key,waveKey,paramKey);
		    
	      //printf("%s %i %i\n",key.c_str(),waveKey,paramKey);
	    }
	  inputFile.get_next();
    }
	    
  inputFile.cur=cur;
      
    
}


optimizePlan buildOptimizePlan(string filename)
{
  xml_input inputFile;
  optimizePlan plan;
  string key;
  double delta;
  xmlNodePtr cur;
  double minParameter=0;
  double maxParameter=0;
  
  //printf("building optimization Plan...\n");
  //printf("Label Wave Param \n");
  inputFile.open(filename);
  
  inputFile.get_child("optimizationPlan");
  
   if (inputFile.get_attribute("delta")!=NULL)
     {
       delta=inputFile.get_real();
     }
   else
     {
       delta=1e-4;
     }
   
   plan.setDelta(delta);
   inputFile.get_child("parameter");
   
   if (!inputFile.check())
     {
       printf("No plan!");
       exit(1);
     }
   
   while( inputFile.check() )
     {// loop on individual wavefunction parameters

       if (inputFile.get_name()=="parameter")
	 {
	   if (inputFile.get_attribute("label")!=NULL)
	     {
	       key=inputFile.get_string();
	     }
      
	   else
	     {
	       printf("Key label error\n");
	       exit(1);
	     }
      
	   if (inputFile.get_attribute("minParameter")!=NULL)
	     {
	       minParameter=inputFile.get_real();
	       plan.setMinParameter(key,minParameter);
	     }
	
	   if (inputFile.get_attribute("maxParameter")!=NULL)
	     {
	       maxParameter=inputFile.get_real();
	       plan.setMaxParameter(key,maxParameter);
	     }
	   
	   buildPlanOneParameter(inputFile,plan,key);
	 }
       inputFile.get_next();
	
     }
    
  #ifdef VERBOSE
  
  plan.print();
  #endif
  return plan;
  
}
