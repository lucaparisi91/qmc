
template<class comp>
speciesHarmonicPotential<comp>* buildSpeciesHarmonicPotential(xml_input* main_input,typename speciesHarmonicPotential<comp>::geometry_t * geo)
{
  xmlNodePtr cur;
  vector<double> omegas;
  vector<double> centers;
  vector<int> sets;
  cur=main_input->cur;
  
  if (main_input->get_child("oneBodyPotential"))
    {
      do
	{
	 
	    
	      
	  omegas.push_back(main_input->get_attribute("omega")->get_real());
	  centers.push_back(main_input->get_attribute("center")->get_real());
	  sets.push_back(main_input->get_attribute("setA")->get_int());
	    
	  main_input->get_next("oneBodyPotential");
	}
      while( main_input->check() );
      
    }
  main_input->cur=cur;
  return new speciesHarmonicPotential<comp>(geo,omegas,centers,sets);
}




template<class comp>
speciesHarmonicPotential<comp>::speciesHarmonicPotential(speciesHarmonicPotential<comp>::geometry_t* geo_,vector<double> &omegas,vector<double> &x0,vector<int> &sets) : potential<comp>(geo_)
  {
    int i=0;
    if (omegas.size() != x0.size())
      {
	cout << "Frequencies and trap centers of different legth."<<endl;
	exit(1);
      }
    
    for(i=0;i<omegas.size();i++)
      {
	traps.push_back(new harmonicPotential(omegas[i],x0[i],sets[i]));
      }
    
  }

template<class comp>
double speciesHarmonicPotential<comp>::evaluate(speciesHarmonicPotential<comp>::all_particles_t* p)
  {
    double ev;
    int i;
    ev=0;
    
    for(i=0;i<traps.size();i++)
    {
      ev+=traps[i]->evaluate(this->geo,p);
    }
    return ev;
  }

template<class comp>
rabiCouplingPotential<comp>* buildRabiPotential(xml_input* main_input,typename comp::geometry_t * geo,typename comp::wave_t * wave)
{
  xmlNodePtr cur;
  double omega;
  int setA;
  cur=main_input->cur;
  
  if (main_input->get_child("oneBodyPotential"))
    { 
	  omega=main_input->get_attribute("omega")->get_real();
	  
	  setA=main_input->get_attribute("setA")->get_int();
	    
	  main_input->get_next("oneBodyPotential");
	
      
    }
  else
    {
      cout << "No one body potential"<<endl;
      exit(1);
    }
  
  main_input->cur=cur;
  
  return new rabiCouplingPotential<comp>(-omega,*geo,*wave);
}
