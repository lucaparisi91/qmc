#include "ptools.h"
#include "tools.h"
#include <cmath>

template<class qt>
center_of_mass_w<qt>* build_center_of_mass_w(xml_input* xml_m)
{
  int bins;
  int setA;
  string label;
  if (xml_m->get_attribute("bins") != NULL)
    {
      bins=xml_m->get_int();
    }
  else
    {
      bins=100;
    }
  if (xml_m->get_attribute("setA") != NULL)
    {
      setA=xml_m->get_int();
    }
  else
    {
      setA=0;
    }
  
  return new center_of_mass_w<qt>(bins,setA);
  
}

template<class mt>
void space_measure<mt>::add(double value,double x,double step)
{
  mt::add(value,grid->index(x),step);
}

template<class mt>
void space_measure<mt>::increment_value(double value,double x,double step)
{
  if (x > grid->a and x< grid->b)
    {
      mt::increment_value(value,grid->index(x),step);
    }
}
template<class mt>
space_measure<mt>::space_measure(double a,double b,int len_,string label_,bool history_) : mt(len_,label_,history_)
{
  // loads external and the number of bits
  // and create a new mesh
  
  grid=new mesh(a,b,len_);
}

template<class mt>
space_measure<mt>* build_space_vector(xml_input* xml_m,string label,double x1,double x2)
{
  int setA,setB,bins,jumps,skip;
  bool history;
  space_measure<mt>* m; // pointer to the space vector measurement
  
  if (xml_m->get_attribute("setA") != NULL)
    {
      setA=xml_m->get_int();
    }
  else
    {
      setA=0;
    }
  
  if (xml_m->get_attribute("setB") != NULL)
    {
      setB=xml_m->get_int();
    }
  else
    {
      setB=0;
    }
  if (xml_m->get_attribute("history") != NULL)
    {
      history=xml_m->get_bool();
    }
  else
    {
      history="false";
    }
  
  if (xml_m->get_attribute("bins") != NULL)
    {
      bins=xml_m->get_int();
    }
  else
    {
      bins=100;
    }

  if (xml_m->get_attribute("jumps") != NULL)
    {
      jumps=xml_m->get_int();
    }
  else
    {
      jumps=0;
    }
  if (xml_m->get_attribute("skip") != NULL)
    {
      skip=xml_m->get_int();
    }
  else
    {
      skip=0;
    }
  
  m=new space_measure<mt>(x1,x2,bins,label,history);
  m->set_a=setA;
  m->set_b=setB;
  m->set_history(history);
  m->set(jumps,skip);
  return m;
}


template<class comp>
void measure_dynamic_scalar<comp>::pack(packed_data* md_pack)
{
  md_pack->pack(last,1);
  md_pack->pack(bins,1);
  md_pack->pack(filled,1);
  md_pack->pack(ms.front(),ms.size());
}

template<class comp>
void measure_dynamic_scalar<comp>::unpack(packed_data* md_pack)
{
  md_pack->unpack(last,1);
  md_pack->unpack(bins,1);
  md_pack->unpack(filled,1);
  md_pack->unpack(ms.front(),ms.size());
}

template<class comp>
int measure_dynamic_scalar<comp>::get_pack_size()
{
  return pTools::get_pack_size(last) + pTools::get_pack_size(bins) + pTools::get_pack_size(filled) + pTools::get_pack_size(ms);
}

template<class qt>
void measure_dynamic_scalar<qt>::clear()
{
  int i;
  ofstream f;
  filled=0;
  last=0;
  ms[0]=ms[bins-1];
  //f.open("walker.hist",ios::app);
  //for(i=0;i<bins;i++)
  //  {
  //    f<<i<<" "<<ms[i]<<endl;
  //  }
  //f.close();
  //exit(0);
}

template <class T>
void measure_dynamic_scalar<T>::add(double m1)
{
  // checks wethever the vector was filled
  
  // adds the last bins to the systems
  
  last=(last+1)%bins;
  
  ms[last]=m1;
  
  if (last== (bins-1))
    {
      filled=1;
      
    }
  
}

template<class qt>
measure_dynamic_scalar<qt>::measure_dynamic_scalar(const int &bins_,const int &set) : measure_dynamic<qt>(set)
{
  bins=bins_;
  ms.resize(bins_);
  last=-1;
  filled=0;
}

template <class qt>
void measure_dynamic_scalar<qt>::print()
{
  int i;
  for (i=0;i<bins;i++)
    {
      cout<<i<<" "<<ms[i]<<endl;
    }
  cout<<endl;
}

template<class T>
measure_dynamic<T>& measure_dynamic_scalar<T>::operator=(measure_dynamic<T>& m1)
{
  int i;
  measure_dynamic_scalar<T>* m;
  m=static_cast< measure_dynamic_scalar<T> * >( &m1);
  //cout << "coping";
  
  for (i=0;i<bins;i++)
    {
      ms[i]=m->ms[i];
    }
  last=m->last;
  filled=m->filled;
}

template<class T>
void measure_dynamic_scalar<T>::reset()
{
  last=-1;
  filled=0;
}
template<class T>
void estimator_decorrelator_vector<T>::add(const vector<T> & m)
{
  int i;
  for(i=0;i<m.size();i++)
    {
      add(m[i],i);
    }
}

template<class T>
void estimator_decorrelator_vector<T>::add(const T& m,const int& j,const int& i)
{
  
  if (i>=depth[j])
    {
      
      depth[j]=depth[j]+1;
      sums[j].push_back(zero_el);
      sum_squares[j].push_back(zero_el);
      n[j].push_back(0);
      last[j].push_back(zero_el);
      waiting[j].push_back(0);
    }
  
  if(waiting[j][i]==1) // if i am waiting
    {
      
      sums[j][i]=sums[j][i]+m;// copy the added member
      sum_squares[j][i]=sum_squares[j][i] + m*m;
      n[j][i]=n[j][i]+1;
      waiting[j][i]=0;// stops waiting
      add((last[j][i] + m)/2,j,i+1);
      last[j][i]=zero_el;
    }
  else
    {
      last[j][i]=m;
      n[j][i]+=1;
      sums[j][i]+=m;
      sum_squares[j][i]+=  m*m;
      waiting[j][i]=1;
      
    }
}


template<class T>
void estimator_decorrelator<T>::add(const T &m,const int& i)
{
  
  if (i>=depth)
    {
      
      depth=depth+1;
      sums.push_back(zero_el);
      sum_squares.push_back(zero_el);
      n.push_back(0);
      last.push_back(zero_el);
      waiting.push_back(0);
    }
  
  if(waiting[i]==1) // if i am waiting
    {
      
      sums[i]=sums[i]+m;// copy the added member
      sum_squares[i]=sum_squares[i] + m*m;
      n[i]=n[i]+1;
      waiting[i]=0;// stops waiting
      add((last[i] + m)/2,i+1);
      last[i]=zero_el;
    }
  else
    {
      last[i]=m;
      n[i]=n[i]+1;
      sums[i]=sums[i]+m;
      sum_squares[i]=sum_squares[i] + m*m;
      waiting[i]=1;
      
    }
}

template<class T>
void estimator_decorrelator_vector<T>::add(const T &m,const int &j)
{
 
  
  add(m,j,0);
  
}
template<class T>
void estimator_decorrelator<T>::add(const T &m)
{
 
  
  add(m,0);
  
}

// compute the variances associated with a certain estimator
template<class T>
void estimator_decorrelator<T>::check_variances()
{
  int i;
  converged=false;

  if(vars.size()>5)
    {
      for(i=1;i<vars.size()-5;i++)
	{
	  if( vars[i]!=0 and (vars[i]-vars[i-1])/vars[i] < eps )
	    {
	      converged=true;
	  
	    }
	}
  
  
      error= vars[vars.size()-6];
    }
  
}

template<class T>
void estimator_decorrelator_vector<T>::check_variances()
{
  int i,j;
  
  for(j=0;j<depth.size();j++)
    {
      converged[j]=false;
      if(vars[j].size()>5)
	{
	  for(i=1;i<vars[j].size()-5;i++)
	    {
	      if( vars[j][i]!=0 and (vars[j][i]-vars[j][i-1])/vars[j][i] < eps )
		{
		  converged[j]=true;
	  
		}
	    }
  
	  
	  errors[j]= vars[j][vars[j].size()-6];
	}
    }
}



template<class T>
void estimator_decorrelator<T>::variances()
{
  int i;
  
  vars.resize(depth);
  
  //cout<<depth<<endl;
  for(i=0;i<depth;i++)
    {
      //cout << n[i]<<","<<sums[i]/n[i]<<"::";
      if (n[i]>1)
	{
	  vars[i]=sqrt(
	    abs(sum_squares[i]/n[i] - (sums[i]/n[i])*(sums[i]/n[i]))
	    /(n[i]-1)
		       );
	  
	  
	}
      else
	{
	  vars[i]=0;
	}
    }
  
  
  //cout<<endl;

}
template<class T>
void estimator_decorrelator_vector<T>::variances( )
{
  int i,j;
  vars.resize(depth.size());
  for (j=0;j<depth.size();j++)
    {
      vars[j].resize(depth[j]);
  //cout<<depth<<endl;
      for(i=0;i<depth[j];i++)
	{
	  //cout << n[i]<<","<<sums[i]/n[i]<<"::";
      if (n[j][i]>1)
	{
	  
	  vars[j][i]=sqrt(
	    abs(sum_squares[j][i]/n[j][i] - (sums[j][i]/n[j][i])*(sums[j][i]/n[j][i]))
	    /(n[j][i]-1)
		       );
	}	
      else
	{
	 
	  vars[j][i]=0;
	}
      
	}
    }
  

}

template<class T>
void estimator_decorrelator<T>::out(const string &filename)
{
  int i;
 
  ofstream f;
  
  variances(); // set the variances for the objects
  check_variances(); 
  //cout<<endl;
  if (vars.size() > 5)
    {
  f.open((filename + ".vars.dat").c_str());
  
  for(i=0;i<(vars.size()-5);i++)
    {
      f<< i << " "<< (vars[i]) << endl;
    }
  f.close();
    }
  
}
template<class T>
void estimator_decorrelator_vector<T>::out(const string &filename)
{
  int i,j,l;
  ofstream f;
  
  variances(); // set the variances for the objects
  check_variances();
  f.open((filename + ".vars.dat").c_str());
  l=depth[vec_max_int(depth)];
  //for(i=0;i<depth.size();i++)
  //  {
  //    cout<<depth[i]<<" ";
  //  }
  //cout<<endl;
  for(i=0;i<l;i++)
    {
      for(j=0;j<depth.size();j++)
	{
	  if (i>=depth[j] -5)
	    {
	      f<<0<<" ";
	    }
	  else
	    {
	      f<< vars[j][i] << " ";
	    }
	}
      f<<endl;
    }
  f.close();
  
  
}

// decorrelation of the estimators
template<class T>
estimator_decorrelator<T>::estimator_decorrelator()
{
  zero_element(zero_el);
  depth=0;
  converged=false;
  error=zero_el;
  eps=0.05;
}
// decorrelation of the estimators
template<class T>
estimator_decorrelator_vector<T>::estimator_decorrelator_vector(int n_)
{
  int i;
  zero_element(zero_el);
  for(i=0;i<n_;i++)
    {
      depth.push_back(0);
      errors.push_back(0);
      converged.push_back(false);
    }
  sums.resize(n_);
  n.resize(n_);
  last.resize(n_);
  sum_squares.resize(n_);
  waiting.resize(n_);
 
  eps=0.05;
  
  
}

// saves the correlator to a file
template<class T>
void estimator_decorrelator<T>::save(const string &filename)
{
  ofstream f;
  int i;
  
  f.open( (filename + string(".decorrelate.dat")).c_str());
  f<<depth<<endl;
  
  
  for(i=0;i<depth;i++)
    {
      
      f<<(n[i])<<" "<<(sums[i])<<" "<<(sum_squares[i])<<" "<<(last[i])<<" "<<( waiting[i])<<endl;
    }
  
  f.close();
  
}

template<class T>
void estimator_decorrelator_vector<T>::save(const string &filename)
{
  ofstream f;
  int i;
  int j;
  f.open( (filename + string(".decorrelate.dat")).c_str());
  //f<<depth<<endl;
  
  for(j=0;j<depth.size();j++)
    {
      f<<depth[j]<<endl;
      for(i=0;i<depth[j];i++)
	{
	  f<<n[j][i]<<" ";
      
	  f<<sums[j][i]<<" ";
     
	  f<<sum_squares[j][i]<<" ";
	
	  f<<last[j][i]<<" ";
	
          f<<waiting[j][i];
      
          f<<endl;
	}
    }
  f.close();

}

template<class T>
void estimator_decorrelator<T>::load(const string &filename)
{
  ifstream f;
  int i;
  
  f.open( (filename + string(".decorrelate.dat")).c_str());
  if (f.good())
    {
      //cout<<filename<<" ---"<<endl;
      f>>depth;
      //cout<<depth<<endl;
      //cout<<filename<<endl;
      //exit(0);
      for(i=0;i<depth;i++)
	{
	  n.push_back(0);
	  sums.push_back(zero_el);
	  sum_squares.push_back(zero_el);
	  last.push_back(zero_el);
	  waiting.push_back(0);
	  
	  f>>n[i];
	  f>>sums[i];
	  f>>sum_squares[i];
	  f>>last[i];
	  f>>waiting[i];
	  //cout << n[i]<< " "<<sums[i]<< " "<<sum_squares[i] << " " << last[i] << " " << waiting[i]<<endl;
	}
      //exit(0);
      
    }
  error=zero_el;
  converged=false;
  vars.resize(0);
  f.close();
}
template<class T>
void estimator_decorrelator_vector<T>::load(const string &filename)
{
  ifstream f;
  int i;
  int j;
  f.open( (filename + string(".decorrelate.dat")).c_str());
  if (f.good())
    {
      
      
      //cout << "----"<<endl;
      //cout << depth << endl;
      n.resize(depth.size());
      sums.resize(depth.size());
      sum_squares.resize(depth.size());
      last.resize(depth.size());
      waiting.resize(depth.size());
      
      for(j=0;j<depth.size();j++)
       {
         f>>depth[j];
	 n[j].resize(depth[j]);
	 sums[j].resize(depth[j]);
	 sum_squares[j].resize(depth[j]);
	 last[j].resize(depth[j]);
	 waiting[j].resize(depth[j]);
	 for(i=0;i<depth[j];i++)
	   {
	     f>>n[j][i];
      
	     f>>sums[j][i];
     
	     f>>sum_squares[j][i];
	
	     f>>last[j][i];
	
	     f>>waiting[j][i];
      
         
	   }
       }
    }
  f.close();
}

template<class walker_t,class wave_t>
void pair_correlation_m<walker_t,wave_t>::make_measurement(walker_t* w,wave_t* wave)
  {
    
    if (setA==setB)
      {
	wave->pair_correlation_symm(w->state,this->ms,this->ms->grid->b,setA);
      }
    else
      {
	wave->pair_correlation_asymm(w->state,this->ms,this->ms->grid->b,setA,setB);
      }
  };

// makes a Reptation Quantum Monte Carlo measurement
template<class wave_t,class m>
rqmcMeasurement<wave_t,m>::rqmcMeasurement(m* m1,int slices_)
{
  slices=slices_;
  int i=0;
  
  for(i=0;i<slices;i++)
    {
      ms.push_back(new m(*m1));
    }
  
}

template<class wave_t,class m>
void rqmcMeasurement<wave_t,m>::make_measurement(vector<rqmc_walker*> p,wave_t * w)
{
  // symmetric measurement
  int i;
  int j;
  int s;
  j=p.size()/2;
  for(i=0;i<ms.size();i++)
    {
      ms[i]->make_measurement(p[i+j],w);
      ms[i]->make_measurement(p[j-i],w);
    };
  
}

template<class walker_t,class wave_t,class qmc_kind>
class center_of_mass_difference_future_walker_creator
{
public:
  typedef measurementInterface<walker_t,wave_t> measurementInterface_t;
  
  bool isSupported(){return false;}
  
  center_of_mass_difference_future_walker<walker_t,wave_t>* create(xml_input* main_input,string label,int setA,int setB,int id,int nSteps)
  {
    throw notYetSsupported("No dynamic measurement for not DMC");
    
  }
  
  void append(vector<measurementInterface_t * > &ms,int id,xml_input* main_input,string label,int setA,int setB,int nSteps)
  {
    throw notYetSsupported("No dynamic measurement for not DMC");
  }
  
};

template<class walker_t,class wave_t,class qmc_kind>
class structure_future_single_creator
{
public:
  typedef measurementInterface<walker_t,wave_t> measurementInterface_t;
  
  bool isSupported(){return false;}
  
  structure_factor_single_f<walker_t,wave_t>* create(xml_input* main_input,string label,int id,int nFutureWalkers,int bins)
  {
    throw notYetSsupported("No dynamic measurement for not DMC");
    
  }
  
  void append(vector<measurementInterface_t * > &ms,int id,xml_input* main_input,string label,int nSteps,int bins)
  {
    throw notYetSsupported("No dynamic measurement for not DMC");
  }
  
};


template<class walker_t,class wave_t>
class structure_future_single_creator<walker_t,wave_t,dmc_t>
{
public:
  typedef measurementInterface<walker_t,wave_t> measurementInterface_t;
  
  bool isSupported(){return true;}
  
  structure_factor_single_f<walker_t,wave_t>* create(xml_input* main_input,string label,int id,int nFutureWalkers,int bins)
  {
    return new structure_factor_single_f<walker_t,wave_t>(build_measure_vector(main_input,label),id,nFutureWalkers,bins);
    
  }
  
  void append(vector<measurementInterface_t * > &ms,int id,xml_input* main_input,string label,int nSteps,int bins)
  {
    ms.append(create(main_input,label,id,nSteps,bins));
  }
  
};

template<class walker_t,class wave_t,class qmc_kind>
class center_of_mass_differenceSquared_future_walker_creator
{
public:
  typedef measurementInterface<walker_t,wave_t> measurementInterface_t;
  
  bool isSupported(){return false;}
  
  center_of_mass_differenceSquared_future_walker<walker_t,wave_t>* create(xml_input* main_input,string label,int setA,int setB,int id,int nSteps)
  {
    throw notYetSsupported("No dynamic measurement for not DMC");
    
  }
  
  void append(vector<measurementInterface_t * > &ms,int id,xml_input* main_input,string label,int setA,int setB,int nSteps)
  {
    throw notYetSsupported("No dynamic measurement for not DMC");
  }
  
};

template<class walker_t,class wave_t>
class center_of_mass_difference_future_walker_creator<walker_t,wave_t,dmc_t>
{
public:
  typedef measurementInterface<walker_t,wave_t> measurementInterface_t;
  
  bool isSupported(){return true;}
  center_of_mass_difference_future_walker<walker_t,wave_t>* create(xml_input* main_input,string label,int setA,int setB,int id,int nSteps)
  {
    return new center_of_mass_difference_future_walker<walker_t,wave_t>(build_measure_scalar(main_input,label),setA,setB,id,nSteps);
    
  }
  
  void append(vector<measurementInterface_t * > &ms,int id,xml_input* main_input,string label,int setA,int setB,int nSteps)
  {
    ms.push_back(create(main_input,label,setA,setB,id,nSteps));
  }
  
};

template<class walker_t,class wave_t>
class center_of_mass_differenceSquared_future_walker_creator<walker_t,wave_t,dmc_t>
{
public:
  typedef measurementInterface<walker_t,wave_t> measurementInterface_t;
  
  bool isSupported(){return true;}
  center_of_mass_differenceSquared_future_walker<walker_t,wave_t>* create(xml_input* main_input,string label,int setA,int setB,int id,int nSteps)
  {
    return new center_of_mass_differenceSquared_future_walker<walker_t,wave_t>(build_measure_scalar(main_input,label),setA,setB,id);
    
  }
  
  void append(vector<measurementInterface_t * > &ms,int id,xml_input* main_input,string label,int setA,int setB,int nSteps)
  {
    ms.push_back(create(main_input,label,setA,setB,id,nSteps));
  }
  
};


template<class walker_t,class wave_t>
class  winding_number_creator<walker_t,wave_t,dmc_t>
{
  typedef measurementInterface<walker_t,wave_t> measurementInterface_t;
public:
  bool isSupported(){return true;};
  
  measurementInterface_t* create(int id,xml_input* main_input,string label)
    {
      int set_a;
      

      if (main_input->get_attribute("setA")!= NULL)
	{
	  set_a=main_input->get_int();
	}
      else
	{
	  throw bad_xml();
	}
      
      if (main_input->get_attribute("op")!= NULL)
	{
	  int set_b;
	  string op;
	  op=main_input->get_string();
	  if (main_input->get_attribute("setB")!=NULL)
	    {
	      set_b=main_input->get_int();
	    }
	  else
	    {
	      throw bad_xml();
	    }

	  if (op=="sum")
	    {
	      return new winding_number_time_sum<walker_t,wave_t>(id,build_measure_vector_mult_index(main_input,label),set_a,set_b );
	    }
	  else
	    {
	      if(op=="diff")
		{
		  return new winding_number_time_diff<walker_t,wave_t>(id,build_measure_vector_mult_index(main_input,label),set_a,set_b );
		}
	      else
		{
		  throw bad_xml();
		}
	    }
	}
      else
	{
	  return new winding_number_time<walker_t,wave_t>(id,build_measure_vector_mult_index(main_input,label),set_a );
	}
    };
  
  void append(vector<measurementInterface_t * > &ms,int id,xml_input* main_input,string label)
  {
    ms.push_back(create(id,main_input,label));
  }
  
};
double getDeltaQ(double l_box,double qMax,int bins);
vector<double> build_q_vector(int bins,double l_box,double qMax);

#include "measuresBuilds.hpp"

template<class tm>
measures<tm>::measures(string filename,tm *qmc_obj)
{
  xml_input* main_input;
  double n_particles=0;
  bool history=false;
  unsigned int i,id=0;
  int bins=0;
  int slices=0;
  bool futureWalkers=false;
  int nFutureWalkers=0;
  string label="Unkown";
  int set_a,set_b=0;
  double max=0;
  string qmc_kind="";
  int nMCM=0;
  double lBox;
  center_of_mass_difference_future_walker_creator<typename tm::walker_t,typename tm::wave_t,typename tm::qmcKind> center_of_mass_dyn_creator;
  center_of_mass_differenceSquared_future_walker_creator<typename tm::walker_t,typename tm::wave_t,typename tm::qmcKind> center_of_mass_dynSquared_creator;
  
  typedef typename tm::walker_t walker_t;
  typedef typename tm::wave_t wave_t;
  
  main_input=new xml_input;
  main_input->open("input.xml");
  bins=0;
  id=0;
  n_particles=qmc_obj->geo->l_box;
  lBox=qmc_obj->geo->l_box;
  bool complex;
  
  ms.push_back(new get_energy_dmc<walker_t,wave_t>("energy"));
  ms.push_back(new get_f_energy_dmc<walker_t,wave_t>("force_energy"));
  
  bool optimize;
  bool spin;
  
  if (
      main_input->reset()->get_child("method")->get_attribute("kind")!= NULL)
    {
      qmc_kind=main_input->get_string();
    }
  else
    {
      cout << "Invalid input file."<<endl;
      cout << "The kind of QMC calculation was not found."<<endl;
      cout<<"Stopping."<<endl;
      exit(1);
    }
  // check whatever the function must be optimized
  
  if (main_input->get_attribute("optimize")!=NULL)
    {
      optimize=main_input->get_bool();
    }
  
  i=1;
  main_input->reset()->get_child("measures")->get_first_child();

  
  while( main_input->check() )
    {
      i=i+1;
      // density
      if (main_input->get_name() == "center_of_mass")
	{
	  if (main_input->get_attribute("label") != NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="center_of_mass";
	    }
	  
	  
	  ms.push_back(new center_of_mass<walker_t,wave_t>(build_measure_scalar(main_input,label)));
	  
	}
      if (main_input->get_name() == "magnetization")
	{
	  if (main_input->get_attribute("label") != NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="magnetization";
	    }
	  
	  ms.push_back(new magnetizationMeasurement<walker_t,wave_t>(build_measure_scalar(main_input,label)));
	  
	}
	    
      if (main_input->get_name() == "centerOfMassDifference")
	{
	  if (main_input->get_attribute("label") != NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="centerOfMassDifference";
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
	  
	  // collect info on future walkers
	  if (main_input->get_attribute("futureWalkers") != NULL)
	    {
	      futureWalkers=main_input->get_bool();
	    }
	  else
	    {
	      futureWalkers=false;
	    }
	  
	  if (futureWalkers==false)
	    {
	      ms.push_back(new centerOfMassDifference<walker_t,wave_t>(build_measure_scalar(main_input,label)));
	    }
	  else
	    {
	      
	      ms.push_back(new centerOfMassDifferenceForwardWalking<walker_t,wave_t>(build_measure_scalar(main_input,label),id));
	      id++;
	    }
	  
       
	}

      if (main_input->get_name() == "centerOfMassDifferenceSquared")
	{
	  if (main_input->get_attribute("label") != NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="centerOfMassDifferenceSquared";
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
	  
	  // collect info on future walkers
	  if (main_input->get_attribute("futureWalkers") != NULL)
	    {
	      futureWalkers=main_input->get_bool();
	    }
	  else
	    {
	      futureWalkers=false;
	    }
	  
	  if (futureWalkers==false)
	    {
	      ms.push_back(new centerOfMassDifferenceSquared<walker_t,wave_t>(build_measure_scalar(main_input,label)));
	    }
	  else
	    { 
	      ms.push_back(new centerOfMassDifferenceSquaredForwardWalking<walker_t,wave_t>(build_measure_scalar(main_input,label),id));
	      id++;
	      
	    }
	  
       
	}


      

      if (main_input->get_name() == "meanSquaresEstimator")
	{
	  if (main_input->get_attribute("label") != NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="meanSquaresEstimator";
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
	  
	  // collect info on future walkers
	  if (main_input->get_attribute("futureWalkers") != NULL)
	    {
	      futureWalkers=main_input->get_bool();
	    }
	  else
	    {
	      futureWalkers=false;
	    }
	  
	  if (futureWalkers==false)
	    {
	      ms.push_back(new meanSquaresEstimator<walker_t,wave_t>(build_measure_scalar(main_input,label)));
	    }
	  else
	    {
	      ms.push_back(new meanSquaresEstimatorForwardWalking<walker_t,wave_t>(build_measure_scalar(main_input,label),id));
	      id++;
	    }
	  
	}
            
      if (main_input->get_name() == "density")
	{
	  bool totalDensity;
	  string centering;
	  
  	  if (main_input->get_attribute("label") != NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="density";
	    }

	  if (main_input->get_attribute("bins") != NULL)
	    {
	      bins=main_input->get_int();
	    }
	  else
	    {
	      bins=1000;
	    }

	  if (main_input->get_attribute("futureWalkers") != NULL)
	    {
	      futureWalkers=main_input->get_bool();
	    }
	  else
	    {
	      futureWalkers=false;
	    }

	  
	  if (main_input->get_attribute("centering") != NULL)
	    {
	      centering=main_input->get_string();
	    }
	  else
	    {
	      centering="none";
	    }

	  if (centering=="none")
	    {
	      // add a new density measurement
	      ms.push_back(new density<walker_t,wave_t>(build_space_vector<measure_vector>(main_input,label,-n_particles/2,n_particles/2)));
	    }
	  else if (centering=="cm")
	    {
	      if(!futureWalkers)
		{
		  ms.push_back(new densityTotalCentered<walker_t,wave_t>(build_space_vector<measure_vector>(main_input,label,-n_particles/2,n_particles/2)));
		}
	      else
		
		{
		  ms.push_back(new densityTotalCenteredForwardWalking<walker_t,wave_t>(build_measure_vector(main_input,label),bins,lBox,id));
		  
		  id++;
		}
	    }
	}
      
      if (main_input->get_name()=="energyHistogram")
	{
	  if (main_input->get_attribute("label") != NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="energyHistogram";
	    }
	  
	  if (main_input->get_attribute("setA") != NULL)
	    {
	      set_a=main_input->get_int();
	    }
	  else
	    {
	      cout << "No impurity set specified"<<endl;
	    }
	  
	  ms.push_back(new energyHistogramImpurity<walker_t,wave_t>(build_space_vector<measure_vector_mult_index>(main_input,label,0,n_particles),set_a));
	  
	}
      // pair correlation
      if (main_input->get_name() == "pair_correlation")
	{
	  ms.push_back(build_pair_correlation<tm>(main_input,n_particles,id));
	  
	}
	  // add a new pair correlation to the system
      
      
      if (main_input->get_name() == "static_structure_factor")
	{
	  ms.push_back( build_structure_factor<tm>(main_input,n_particles,id) );
	  if (futureWalkers==true) {id++;}
	}
      
      
      
      if (main_input->get_name() == "oneBodyDensityMatrix")
	{
	  if (main_input->get_attribute("max") != NULL)
	    {
	      max=main_input->get_real();
	    }
	  else
	    {
	      max=n_particles/2;
	    }

	  if (main_input->get_attribute("futureWalkers")!= NULL)
	    {
	      futureWalkers=main_input->get_bool();
	    }
	  else
	    {
	      futureWalkers=false;
	    }
	  
	  if(futureWalkers)
	    {
	      if (main_input->get_attribute("nFutureWalkers")!= NULL)
		{
		  nFutureWalkers=main_input->get_int();
		}
	      else
		{
		  cout << "Please specify the number of future walkers for the 1BDM."<<endl;
		  exit(1);
		}
	    }
	  
	  if (main_input->get_attribute("label")!=NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="oneBodyDensityMatrix";
	    }
	  
	  if (main_input->get_attribute("nMCM")!=NULL)
	    {
	      nMCM=main_input->get_int();
	    }
	  else
	    {
	      nMCM=100;
	    }

	  
	  if (!futureWalkers)
	    {
	      // add a new pair correlation to the system
	      ms.push_back(new oneBodyDensityMatrixOffdiagonal<walker_t,wave_t>(build_space_vector<measure_vector_mult_index>(main_input,label,0,max),nMCM,max));
	    }
	  else
	    {
	      ms.push_back(new oneBodyDensityMatrixOffdiagonalFutureWalkers<tm>(build_space_vector<measure_vector_mult_index>(main_input,label,0,max),nMCM,max,nFutureWalkers,qmc_obj));
	    }
	  
	}
      
      // winding number
      
      if (main_input->get_name() == "winding_number" and (qmc_kind=="dmc" or qmc_kind=="svmc"))
	{
	  if (main_input->get_attribute("label")!=NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="winding_number";
	    }
	  
	  ms.push_back(buildWindingNumber<walker_t,wave_t>(main_input,id));
	  id++;
	  
	}

      if (main_input->get_name() == "winding_number_spin" and (qmc_kind=="dmc" or qmc_kind=="svmc"))
	{
	  if (main_input->get_attribute("label")!=NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="winding_number_spin";
	    }
	  
	  ms.push_back(buildWindingNumberSpin<walker_t,wave_t>(main_input,id));
	  id++;
	}

      if (main_input->get_name() == "centerOfMassDifferenceImaginaryTimeCorrelation" and (qmc_kind=="dmc" or qmc_kind=="svmc"))
	{
	  if (main_input->get_attribute("label")!=NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="centerOfMassDifferenceImaginaryTimeCorrelation";
	    }
	  
	  ms.push_back(buildCenterOfMassDifferenceImaginaryTimeCorrelation<walker_t,wave_t>(main_input,id));
	  id++;
	}
      
      if (main_input->get_name() == "centerOfMassSumSquaredImaginaryTimeCorrelation")
	{
	  if (main_input->get_attribute("label") != NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="centerOfMassSumSquaredImaginaryTimeCorrelation";
	    }
	  ms.push_back(buildCenterOfMassSumSquaredImaginaryTimeCorrelation<walker_t,wave_t>(main_input,id));
	  id++;
	  
	}
      
      if (main_input->get_name() == "centerOfMassSpinDifferenceSquared")
	{
	  if (main_input->get_attribute("label") != NULL)
	    {
	      label=main_input->get_string();
	    }
	  else
	    {
	      label="centerOfMassSpinDifferenceSquared";
	    }
	  
	  
	  
	  if (main_input->get_attribute("futureWalkers") != NULL)
	    {
	      futureWalkers=main_input->get_bool();
	    }
	  else
	    {
	      futureWalkers=false;
	    }
	  
	  
	  
	  if (futureWalkers==false)
	    {
	      ms.push_back(new centerOfMassSpinDifferenceSquaredMeasurement<walker_t,wave_t>(build_measure_scalar(main_input,label)));
	    }
	  else
	    {
	      ms.push_back(new centerOfMassSpinDifferenceSquaredMeasurementForwardWalking<walker_t,wave_t>(build_measure_scalar(main_input,label),id));
	      id++;
	    }
	}
      
      
      main_input->get_next();
      
    }
  
  
  
      
      // if (main_input->get_attribute("time_slices")!= NULL)
      // 	{
      // 	  slices=main_input->get_int();
      // 	  ms[ms.size()-1]->set_slices(slices);
      // 	}
      
      
      
	
    
  
  for(i=0;i<ms.size();i++)
    {
      ms[i]->load();
    }
  
}

template<class qmc_t>
void measures<qmc_t>::out()
{
  // outputs some form of measurements
  unsigned int i=0;
  
  for(i=0;i<ms.size();i++)
    {
      ms[i]->out();
    }
}

template<class qmc_t>
void measures<qmc_t>::clear()
{
  unsigned int i;
  for(i=0;i<ms.size();i++)
    {
      ms[i]->clear();
      ms[i]->clear_history();
    }
}

template<class tm>
void measures<tm>::make_measurements(measures<tm>::measure_obj_t* w,measures<tm>::wave_t* wave)
{
  
  //ms[0]->make_measurement(w,wave);
  //ms[1]->make_measurement(w,wave);
  unsigned int i;
  for(i=0;i<ms.size();i++)
    {
     
      if (ms[i]->check() == 2)
	{
	  ms[i]->make_measurement(w,wave);
	} 
    }
}

template<class qmc_t>
void measures<qmc_t>::save()
{
  unsigned int i;
  for(i=0;i<ms.size();i++)
    {
      ms[i]->save();
    }
}
template<class qmc_t>
void measures<qmc_t>::increment()
{
  unsigned int i;
  for(i=0;i<ms.size();i++)
    {
      ms[i]->increment();
    }
}

template<class qmc_t>
void measures<qmc_t>::record(double step_)
{
  unsigned int i;
  
  for(i=0;i<ms.size();i++)
    {
      ms[i]->record(step_);
    }
    
}

template<class tm>
void measures<tm>::reduce()
{

  unsigned int i;
  
  //ms[0]->reduce(0);
  //ms[1]->reduce(0);
  //ms[2]->reduce(0);
  //ms[3]->reduce(0);
  //ms[4]->reduce(0);
  for(i=0;i<ms.size();i++)
    {
      ms[i]->reduce(0);
    }

}
