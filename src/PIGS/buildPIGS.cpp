#include "buildPIGS.h"
#include "pigsTools.h"

void pigsBuilder::build(pigsBuilder::configurations_t & configurations)
{
    int nBeads=0;
    
    string kind="random";
    pos_t boxLength;
    
    boxLength=getLbox();

    nBeads=getNBeads();
    
    configurations.resize(nBeads);
    
    buildAllOrbitals("input.xml",configurations[0]);
    generateRandom(configurations[0],-boxLength/2.,boxLength/2.,getRandomGenerator());
    getGeometry()->all_particles_pbc(configurations[0]);
    
    for(int i=1;i<nBeads;i++)
      {
	configurations[i]=configurations[0];
      }
};

void pigsBuilder::getParameter(string name,double & value)
{
  
  value=qmcO->main_input->reset()->get_child("system")->get_child(name)->get_value()->get_real(); 
}

void pigsBuilder::build(pairParticleApproximationChain & prop)
{
  typedef pairParticleApproximationChain::pairParticlePropagator_t pairParticlePropagator_t;
  
  typedef pairParticleApproximationChain::freeParticlePropagator_t freeParticlePropagator_t;
  
  freeParticlePropagator_t * gf;
  pairParticlePropagator_t * gp;
  gf=new freeParticlePropagator_t();
  gp=new pairParticlePropagator_t();
  build(*gf);
  build(*gp);
  prop.setWave(getWave());
  prop.setFreePropagator(gf);
  prop.setPairPropagator(gp);
  prop.setGeometry(getGeometry());
  prop.setTimeStep(getTimeStep());
  
}

void pigsBuilder::build(freePropagator1D & gf)
{
  gf.setTimeStep(getTimeStep());
  gf.setGeometry(getGeometry());
}

void pigsBuilder::build(pairApproximationPropagator1DUnitary & gf)
{
  gf.setTimeStep(getTimeStep());
  gf.setGeometry(getGeometry());
  
}

void pigsBuilder::build(wiggleMove & mov)
{
  levyReconstructor *constructor=new levyReconstructor();
  build(*constructor);
  mov.setConstructor(constructor);
  mov.setRandomGenerator(getRandomGenerator());
}

void pigsBuilder::build(translateMove & mov)
{
  levyReconstructor *constructor=new levyReconstructor();
  build(*constructor);
  mov.setConstructor(constructor);
  mov.setRandomGenerator(getRandomGenerator());
  mov.setDisplacement(0.01);
  
}

void pigsBuilder::build(moveVariationalTails & mov)
{
  mov.setRandomGenerator(getRandomGenerator());
  mov.setGeometry(getGeometry());
  mov.setTimeStep(getTimeStep());
  
}

void pigsBuilder::build(pigsMover & mov)
{
  mov.setRandomGenerator(getRandomGenerator());
  propagator_t *prop;
  prop=new propagator_t();
  build(*prop);
  towerSampling * sampler=new towerSampling();
  build(*sampler);
  mov.setSampler(sampler);
  mov.setPropagator(prop);
  
  wiggleMove * move1=new wiggleMove();
  moveVariationalTails * move2=new moveVariationalTails();
  translateMove * move3=new translateMove();
  
  moveSingleBead * move4=new moveSingleBead();
  build(*move1);
  build(*move2);
  build(*move3);
  build(*move4);
  
  mov.addMove(*move1,0.6);
  mov.addMove(*move2,0.2);
  mov.addMove(*move3,0.2);
  //mov.addMove(*move4,1);
  
}

int pigsBuilder::getNBeads()
{
  return qmcO->main_input->reset()->get_child("method")->get_child("nBeads")->get_value()->get_int();
}

void pigsBuilder::build(levyReconstructor & constructor)
{
  configurations_t * configurationsBackup=new configurations_t();
  
  build(*configurationsBackup);
  constructor.setGeometry(getGeometry());
  constructor.setRandomGenerator(getRandomGenerator());
  constructor.setConfigurationsBackup(configurationsBackup);
  constructor.setTimeStep(getTimeStep());
  
}

void pigsBuilder::build(towerSampling & tower)
{
  tower.setRandomGenerator(getRandomGenerator());
}

void pigsBuilder::build(pairCorrelationPigs & m,xml_input * main_input)
{
  string label;
  double max;
  
  if (main_input->get_attribute("max") != NULL)
    {
      max=main_input->get_real();
    }
  else
    {
      max=getLbox()/2.;
    }
  
  if (main_input->get_attribute("label")!=NULL)
    {
      label=main_input->get_string();
    }
  else
    {
      label="pair_correlation";
    }
  
  space_measure<measure_vector> * ms;
  ms=build_space_vector<measure_vector>(main_input,label,0,max);
  m.setEstimator(ms);
  m.setTimeSliceBegin(getNBeads()/2);
  m.setTimeSliceEnd(getNBeads()/2.+1);
  //m.setTimeSliceBegin(0);
  //m.setTimeSliceEnd(getNBeads());
}

void pigsBuilder::build(energyThermodynamic & m,xml_input * main_input)
{
  string label;
  double max;
  
  
  
  if (main_input->get_attribute("label")!=NULL)
    {
      label=main_input->get_string();
    }
  else
    {
      label="energyThermodynamic";
    }
  
  measure_scalar * ms;
  ms=build_measure_scalar(main_input,label);
  m.setEstimator(ms);
  m.setPropagator(getPropagator());
  //m.setTimeSliceBegin(getNBeads()/2);
  //m.setTimeSliceEnd(getNBeads()/2.+1);
  //m.setTimeSliceBegin(0);
  //m.setTimeSliceEnd(getNBeads());
}

void pigsBuilder::build(moveSingleBead & move)
{
  
  move.setGeometry(getGeometry());
  move.setDisplacement(1.);
  move.setRandomGenerator(getRandomGenerator());
}

void pigsBuilder::build(measures<qmcDriver_t> & ms)
{
  xml_input* main_input=qmcO->main_input;
  
  string label;
  string name;
  int i=0;
  
  main_input->reset()->get_child("measures")->get_first_child();
  printf("Building measurements...\n");
  
  while( main_input->check() )
    {
      i=i+1;
      name=main_input->get_name();
      label=getAttribute(main_input,"label",name);
      if ( name== "energyTail")
	{
	  printf("\t energy...\n");
	  ms.push_back(new energyTail(build_measure_scalar(main_input,label),"input.xml"));
	}
      else if (name=="pair_correlation")
	{
	  printf("\tpair correlation..\n");
	  pairCorrelationPigs * m=new pairCorrelationPigs();
	  build(*m,main_input);
	  ms.push_back(m);
	  
	}
      else if (name=="energyThermodynamic")
	{
	  printf("\t energyThermodynamic.. \n");
	  energyThermodynamic * m=new energyThermodynamic();
	  build(*m,main_input);
	  ms.push_back(m);
	  
	}
      
      main_input->get_next();
    }
  
}

