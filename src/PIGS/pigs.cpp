#include "pigs.h"
#include "logger.h"

bool running;

void pigsDriver::setWavefunction()
{
   vector<swave_t*> waves;
   load_wavefunctions<qmcDriver_t >(this->main_input,waves,this);
   wave=new wave_t(this);
   wave->link_wavefunctions(this->main_input,waves,"main");
   wave->print_jastrows();
}

#include "buildPIGS.h"

void pigsDriver::load()
  {
    pigsBuilder builder(this);
    
    printf("...loading configurations\n");
    builder.build(configurations);
    printf("...loading Moves\n");
    pigsMoverO=new pigsMover_t();
    builder.build(*pigsMoverO);
    
    printf("...setting measurements\n");
    builder.build(estimates);
    
  }

void pigsDriver::out()
{
  estimates.reduce();
  printf( "---------------------\n");
  printf("%s\n",pigsMoverO->info().c_str());
  estimates.out();
  estimates.clear();
  
}

void pigsDriver::run()
{
  printf("PIGS----------\n");
  printf("Warming...\n");
  for(int j=0;j<warmupBlocks;j++)
    {
      for(int i=0;i<stepsPerBlock;i++)
	{
	  step();
	}
      
    }
  
    printf("Measuring...\n");
    for(int j=0;j<nBlocks;j++)
    {
      for(int i=0;i<stepsPerBlock;i++)
	{
	  step();
	  
	  measure();
	}
      
      singletonLogger::getInstance()<<configurations;
      out();
      
    }
  printf("Finished!\n");
  
}

void pigsDriver::measure()
{
  estimates.increment();
  estimates.make_measurements(&configurations,wave);
}

void pigsDriver::step()
{
  pigsMoverO->move(configurations);
}

