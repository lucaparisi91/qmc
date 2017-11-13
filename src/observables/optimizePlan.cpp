#include "optimizePlan.h"
#include "../xml-input.h"

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
      printf("\n");
      printf("Delta: %f \n",getDelta());
    }
  
}

optimizePlan buildOptimizePlan(string filename)
{
  xml_input inputFile;
  optimizePlan plan;
  int paramKey;
  int waveKey;
  string key;
  double delta;
  xmlNodePtr cur;
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
   
  if(inputFile.check())
    {
      inputFile.get_child("parameter");
      
      while( inputFile.check() )
	{// loop on individual wavefunction parameters
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
		    
		    if (inputFile.get_attribute("label")!=NULL)
		      {
			key=inputFile.get_string();
		      }
		    
		    else
		      {
			printf("Key label error\n");
			exit(0);
		      }
		    
		    plan.add(key,waveKey,paramKey);
		    
		    //printf("%s %i %i\n",key.c_str(),waveKey,paramKey);
		  }
		inputFile.get_next();
	      }
	    
	    inputFile.cur=cur;
	    inputFile.get_next();
	}
    }
  else
    {
      printf("No optimization plan !");
    }
  #ifdef VERBOSE
  
  plan.print();
  #endif
  return plan;
  
}
