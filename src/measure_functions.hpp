#include "dmc.h"
#include "jastrow.h"
#include "geometry.h"
#include "system.h"
#include "wavefunction.h"
#include "measures.h"
#include "mesh.h"
#include "random.h"
#include <iostream>
#include <complex>
#include <cassert>
#include <cmath>
#include <cstdlib>

/*
// computes the kinetic energy of the system
void bill_jastrow_wavefunction::kinetic_energy(particles* state,double &e,double &e_f)
{
  int i,j,sign;
  double tmp,tmp1,tmp2,x;
  
  assert(!one_body_jastrows_present or jastrow1b );
  e=0;
  e_f=0;
  for(i=0;i<qmc_obj->geo->n_particles;i++)
    {
      state->drift_force[i]=0;
    }
  
  for(i=0;i<qmc_obj->geo->n_particles;i++)
    {
      for(j=0;j<i;j++)
	{
	  x=qmc_obj->geo->distance_pbc(state->position[i],state->position[j]);
	  
	  if (x>=0)
	    {
	      sign=1;
	    }
	  else
	    {
	      sign=-1;
	      x=-x;
	    }
	  
	  tmp=jastrow2b->d0(x);
	  tmp1=jastrow2b->d1(x);
	  tmp2=jastrow2b->d2(x);
	  
	  state->drift_force[i]=state->drift_force[i]+ sign*tmp1/tmp;
	  state->drift_force[j]=state->drift_force[j]- sign*tmp1/tmp;
	  
	  e=e+2*(tmp2/tmp - pow((tmp1/tmp),2));
	}
      
  
  if (one_body_jastrows_present)
    {
      x=state->position[i];
      if (x>=0)
	     {
              sign= 1;
              }
	  else
	    {
	      sign=-1;
	      x=-x;
	    }
	  
      tmp=jastrow1b->d0(x);
	  tmp1=jastrow1b->d1(x);
	  tmp2=jastrow1b->d2(x);
	  
	  state->drift_force[i]=state->drift_force[i] + sign*tmp1/tmp;
	  e=e+(tmp2/tmp - pow((tmp1/tmp),2));
  
    }
          
    }
  
  for (i=0;i<qmc_obj->geo->n_particles;i++)
    {
      e_f=e_f + pow(state->drift_force[i],2);
    }
  // accumulates the energy and the force energy
  e= e + e_f;
  e=-e/2;
  e_f=e_f/2;
  
}
*/

// measure and record the pair correlation function

template<class comp>
template<class mV_t>
void total_wavefunction<comp>::pair_correlation_asymm(total_wavefunction<comp>::all_particles_t *p, mV_t* g,double max,int set_a,int set_b)
{
  // compute the pair correlation of the system
  int i,j;
  double dg,x;
  particles_t* p1,*p2;
  //cout << max<<endl;
  p1=p->particle_sets[set_a];
  p2=p->particle_sets[set_b];
  //cout << g->set_a << " "<<g->set_b<<endl;
  // loop on all particles and update the frequency bin
  dg=(qmc_obj->geo->l_box/g->get_step())/(p1->n*p2->n);
  
  for (i=0;i<p1->n;i++)  
    {
      for(j=0;j<p2->n;j++)
	{
	  x=abs(qmc_obj->geo->distance_pbc(p1->position[i],p2->position[j]));
	  if(x<=max)
	    {
	      
	      g->increment_value(dg,x,qmc_obj->current_step);
	    }
	}
    }
  g->increment_index();
  
}

template<class comp>
void total_wavefunction<comp>::structure_factor_asymm(total_wavefunction<comp>::all_particles_t *p, measure_vector* s,vector<double> &qs,int setA,int setB)
{
  // compute the pair correlation of the system
  int i,j,k;
  double dg,x;
  particles_t* p1,*p2;
  //cout << max<<endl;
  
  p1=p->particle_sets[setA];
  p2=p->particle_sets[setB];
  
  //cout << g->set_a << " "<<g->set_b<<endl;
  // loop on all particles and update the frequency bin
  for(k=0;i<qs.size();k++)
    {
      for (i=0;i<p1->n;i++)  
	{
	  for(j=0;j<p2->n;j++)
	    {
	      x=qmc_obj->geo->distance_pbc(p1->position[i],p2->position[j]);
	  
	      s->increment_value(2*cos(x*qs[k]),k,qmc_obj->current_step);
	    }
	}
    }
  
  s->increment_index();
  
}

template<class comp>
void total_wavefunction<comp>::structure_factor_symm(total_wavefunction<comp>::all_particles_t *p, measure_vector* s,vector<double> &qs,int setA)
{
  // compute the pair correlation of the system
  int i,j,k;
  double dg,x;
  particles_t* p1,*p2;
  //cout << max<<endl;
  
  p1=p->particle_sets[setA];
  
  //cout << g->set_a << " "<<g->set_b<<endl;
  // loop on all particles and update the frequency bin
  
  for(k=0;k<qs.size();k++)
    {
      for (i=0;i<p1->n;i++)  
	{
	  for(j=0;j<= i;j++)
	    {
	      x=abs(qmc_obj->geo->distance_pbc(p1->position[i],p1->position[j]));
	      
	      s->increment_value(2*cos(x*qs[k]),k,qmc_obj->current_step);
	    }
	  
	}
    }
  
  s->increment_index();
  
}

template<class comp>
template<class measure_vec_t>
void total_wavefunction<comp>::structure_factor(total_wavefunction<comp>::all_particles_t *p, measure_vec_t* s,vector<double> &qs,int setA,vector<complex<double> > & work)
{
  // compute the pair correlation of the system
  int i,k;
  double dg,x;
  particles_t* p1,*p2;
  complex<double> j;
  j=(0,1);
  //cout << max<<endl;
  
  p1=p->particle_sets[setA];
  
  //cout << g->set_a << " "<<g->set_b<<endl;
  // loop on all particles and update the frequency bin
  
  for(k=0;k<qs.size();k++)
    {
      work[k]=complex<double>(0,0);
      
      for (i=0;i<p1->n;i++)  
	{    
	  work[k]+=exp(complex<double>( 0,qs[k]*p1->position[i]))	;
	}
      
      s->increment_value(pow(abs(work[k]),2)/p1->n,k,qmc_obj->current_step);
    }
  
  s->increment_index();
  
}

template<class comp>
template<class measure_vec_t>
void total_wavefunction<comp>::structure_factor(total_wavefunction<comp>::all_particles_t *p, measure_vec_t* s,vector<double> &qs,int setA,int setB,vector<complex<double> > & work)
{
  // compute the pair correlation of the system
  int i,k;
  double n;
  particles_t* p1,*p2;
  complex<double> j;
  j=(0,1);
  //cout << max<<endl;
  
  p1=p->particle_sets[setA];
  p2=p->particle_sets[setB];
  n=p1->n+p2->n;
  //cout << g->set_a << " "<<g->set_b<<endl;
  // loop on all particles and update the frequency bin
  
  for(k=0;k<qs.size();k++)
    {
      work[k]=complex<double>(0,0);
      
      for (i=0;i<p1->n;i++)  
	{    
	  work[k]+=exp(complex<double>( 0,qs[k]*p1->position[i]))	;
	}
      
      for(i=0;i<p2->n;i++)
	{
	  work[k]+=exp(complex<double>( 0,qs[k]*p2->position[i]))	;
	}
      
      s->increment_value(pow(abs(work[k]),2)/n,k,qmc_obj->current_step);
    }
  
  s->increment_index();
  
}

template<class comp>
template<class measure_vec_t>
void total_wavefunction<comp>::structure_factorSpin(total_wavefunction<comp>::all_particles_t *p, measure_vec_t* s,vector<double> &qs,int setA,int setB,vector<complex<double> > & work)
{
  // compute the pair correlation of the system
  int i,k;
  double n;
  particles_t* p1,*p2;
  //cout << max<<endl;
  complex<double> j;
  j=(0,1);
  p1=p->particle_sets[setA];
  p2=p->particle_sets[setB];
  n=p1->n+p2->n;
  //cout << g->set_a << " "<<g->set_b<<endl;
  // loop on all particles and update the frequency bins
  
  for(k=0;k<qs.size();k++)
    {
      work[k]=complex<double>(0,0);
      
      for (i=0;i<p1->n;i++)  
	{    
	  work[k]+=exp(complex<double>( 0,qs[k]*p1->position[i]))	;
	}
      
      for(i=0;i<p2->n;i++)
	{
	  work[k]-=exp(complex<double>( 0,qs[k]*p2->position[i]))	;
	}
      
      s->increment_value(pow(abs(work[k]),2)/n,k,qmc_obj->current_step);
    }
  
  s->increment_index();
  
}
// compute the off diagonal one body density matrix
template<class comp>
void total_wavefunction<comp>::oneBodyDensityMatrixOffdiagonal(total_wavefunction<comp>::all_particles_t * p,space_measure<measure_vector_mult_index>* g,int nMCM,int Max)
{
  int i,j,k;
  particles_t * p1;
  double r,rnew;
  double swap;
  double l,dg;
  double logRatio,logRatioUnmodified,ratio;
  ratio=0;
  l=qmc_obj->geo->l_box;
  dg=1;
  
  p1=p->particle_sets[g->set_a];
  
  for(i=0;i<p1->n;i++)
    {
      // compute the unmodified ratio
      logRatioUnmodified=0;
      
      for(k=0;k<waves.size();k++)
	{
	  // comput the jastrow terms un modified
	  logRatioUnmodified+=waves[k]->one_particle_log_evaluate(p,p1->position[i],i,g->set_a);
	}
      
      // generate distances randomly
      for(j=0;j<nMCM;j++)
	{
	  // generates uniformely in [-l/2,l/2]
	  r=(qmc_obj->rand->uniform())*(g->grid->b - g->grid->a ) + g->grid->a ;
	  rnew=qmc_obj->geo->pbc( p1->position[i] + r);
	  // compute the jastrow terms modified
	  logRatio=0;
	  
	  for(k=0;k<waves.size();k++)
	    {
	      
	      logRatio+=waves[k]->one_particle_log_evaluate(p,rnew,i,g->set_a);
	    }
	  
	  ratio=exp(logRatio-logRatioUnmodified);
	  
	  // add to the accumulator
	  g->add(ratio,r,qmc_obj->current_step);
	  
	  
	}
      
    }
  
  
}

// compute the off diagonal one body density matrix
template<class comp>
void total_wavefunction<comp>::oneBodyDensityMatrixOffdiagonalFutureWalkers(total_wavefunction<comp>::all_particles_t* p,space_measure<measure_vector_mult_index>* g,int nMCM,int Max,miniDMCDriver< comp >* driver)
{
  int i,j,k;
  particles_t * p1;
  double r,rnew;
  double weight;
  double swap;
  double l,dg;
  double logRatio,logRatioUnmodified,ratio;
  ratio=0;
  l=qmc_obj->geo->l_box;
  dg=1;
  p1=p->particle_sets[g->set_a];
  
  for(i=0;i<p1->n;i++)
    {
      // compute the unmodified ratio
      logRatioUnmodified=0;
      
      for(k=0;k<waves.size();k++)
	{
	  // comput the jastrow terms un modified
	  logRatioUnmodified+=waves[k]->one_particle_log_evaluate(p,p1->position[i],i,g->set_a);
	}
      
      // generate distances randomly
      for(j=0;j<nMCM;j++)
	{
	  // generates uniformely in [-l/2,l/2]
	  r=(qmc_obj->rand->uniform())*(g->grid->b - g->grid->a ) + g->grid->a ;
	  rnew=qmc_obj->geo->pbc( p1->position[i] + r);
	  // compute the jastrow terms modified
	  logRatio=0;
	  
	  for(k=0;k<waves.size();k++)
	    {
	      
	      logRatio+=waves[k]->one_particle_log_evaluate(p,rnew,i,g->set_a);
	    }  
	  ratio=exp(logRatio-logRatioUnmodified);
	  // set the initial configurations
	  driver->set(p);
	  // makes the driver evolve over time
	  driver->evolve();
	  // returns the weight gathered along the evolution
	  weight=driver->get_weight();
	  // add to the accumulator
	  g->add(ratio,r,weight,qmc_obj->current_step);
	  
	} 
    }
}

template<class comp>
template<class mV_t>
void total_wavefunction<comp>::pair_correlation_symm(total_wavefunction<comp>::all_particles_t *p, mV_t* g,double max,int set_a)
{
  // compute the pair correlation of the system
  int i,j;
  double dg,x;
  particles_t* p1;
  
  p1=p->particle_sets[set_a];
  
  // loop on all particles and update the frequency bin
  dg=(qmc_obj->geo->l_box/g->get_step())/pow(qmc_obj->geo->l_box,2);
  
  for (i=0;i<p1->n;i++)
    {
      for(j=0;j<i;j++)
	{
	  x=abs(qmc_obj->geo->distance_pbc(p1->position[i],p1->position[j]));
	  if(x<=max)
	    {
	      g->increment_value(dg,x,qmc_obj->current_step);
	    }
	}
    }
  g->increment_index();
}

// measurement of the center of mass of the original system
template<class comp>
void total_wavefunction<comp>::center_of_mass_no_pbc(total_wavefunction<comp>::all_particles_t* p,measure_scalar* winding)
{
  int i;
  double mean_position;
  particles_t* state;
  state=p->particle_sets[winding->set_a];
  
  mean_position=0;
  // loop on all particle positions
  for(i=0;i<state->position_no_pbc.size();i++)
    {
      mean_position+=state->position_no_pbc[i];
    }
  
  mean_position=mean_position/state->position_no_pbc.size();
  
  winding->add(mean_position,qmc_obj->current_step);
  
}

// add a measurement of the density
template<class comp>
void total_wavefunction<comp>::density(total_wavefunction<comp>::all_particles_t *p,space_measure<measure_vector>* g)
{
  int i;
  double x;
  particles_t* state;
  double alpha;
  
  state=p->particle_sets[g->set_a];
  
  for(i=0;i<state->position.size();i++)
    {
      x=state->position[i];
      
      g->increment_value(1./(g->grid->step),x,qmc_obj->current_step);
    }
  g->increment_index();
}


// // overlap with the total wavefunction stored in o
// void total_wavefunction::overlap(all_particles* p,overlap_measure* o)
// {
//   // two different values for the wavefunction
//   double w1, w2;
//   w1=log_evaluate(p);
//   w2=o->wave_o->log_evaluate(p);
//   //cout<< w2<< " "<<w1<<" "<<exp(w2-w1)<<endl;
//   // add the overlap value
//   o->add(exp(w2-w1),qmc_obj->current_step);
// }

// compute the center of mass mithout periodic boundary conditions
template<class comp>
double total_wavefunction<comp>::center_of_mass_no_pbc(total_wavefunction<comp>::all_particles_t* p,const int set)
{
  int i;
  particles_t* p1;
  double r;
  
  p1=p->particle_sets[set];
  r=0;
  for(i=0;i<p1->n;i++)
    {
      r+=p1->position_no_pbc[i];
    }
  return r/p1->n;
}

// template class total_wavefunction<dmc<D1_t> >;
// template class total_wavefunction<dmc<spinor1D> >;
// template class bill_jastrow_wavefunction<jastrow_delta,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_two_body_symmetric<jastrow_delta,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_two_body_symmetric<jastrow_delta_bound_state,dmc<D1_t> >;

// template class bill_jastrow_wavefunction_one_body<jastrow_delta,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_one_body<jastrow_delta_bound_state,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_two_body_asymmetric<jastrow_delta,dmc<D1_t> >;
// template class bill_jastrow_wavefunction_two_body_asymmetric<jastrow_delta_bound_state,dmc<D1_t> >;
