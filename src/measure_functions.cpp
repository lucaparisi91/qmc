#include <vector>
#include "measures.h"
#include "particles/orbital.h"

// using namespace std;

void structureFactorSpin(orbitals<spinOrbital1D> & p1, vector<double> &res,const vector<double> &qs,vector<complex<double> > &work)  
{
  // compute the pair correlation of the system
  int i,k;
  double dg,x;
  
  typedef orbitals<spinOrbital1D>  particles_t;
  
  complex<double> j;
  j=(0,1);
  //cout << max<<endl;
  
  
  
  //cout << g->set_a << " "<<g->set_b<<endl;
  // loop on all particles and update the frequency bin
  
  for(k=0;k<qs.size();k++)
    {
      work[k]=complex<double>(0,0);
      
      for (i=0;i<p1.size();i++)  
	{    
	  work[k]+=complex<double>(p1[i].spin())*exp(complex<double>( 0,qs[k]*p1[i].position()))	;
	}
      
      res[k]=pow(abs(work[k]),2)/p1.size();
      
     }
  
}

void structureFactorDensity(orbitals<spinOrbital1D> & p1, vector<double> &res,const vector<double> &qs,vector<complex<double> > &work)  
{
  // compute the pair correlation of the system
  int i,k;
  double dg,x;
  
  typedef orbitals<spinOrbital1D>  particles_t;
  
  complex<double> j;
  j=(0,1);
  //cout << max<<endl;
  
  
  
  //cout << g->set_a << " "<<g->set_b<<endl;
  // loop on all particles and update the frequency bin
  
  for(k=0;k<qs.size();k++)
    {
      work[k]=complex<double>(0,0);
      
      for (i=0;i<p1.size();i++)  
	{    
	  work[k]+=exp(complex<double>( 0,qs[k]*p1[i].position()))	;
	}
      
      res[k]=pow(abs(work[k]),2)/p1.size();
      
     }
  
}

double centerOfMassNoBC(orbitals<spinOrbital1D> & p1)
{
  double centerOfMass=0;
  
  for(int i=0;i<p1.size();i++)
    {
      centerOfMass+=p1[i].positionNoBC();
    }
  
  return centerOfMass/p1.size();
}

double centerOfMassSpinNoBC(orbitals<spinOrbital1D> & p1)
{
  double centerOfMass=0;
  
  for(int i=0;i<p1.size();i++)
    {
      centerOfMass+=p1[i].positionNoBC()*p1[i].spin();
    }
  
  return centerOfMass/p1.size();
}
