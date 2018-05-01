#include "pigs.h"
#include "PIGSpropagator.h"

double pairParticleApproximationChain::logEvaluateSingleLinks(pairParticleApproximationChain::configurations_t & configurations)
{
  double logProb=0;
  
  for(int i=0;i<configurations.size()-1;i++)
    {
      for(int j=0;j<configurations[i][0].size();j++)
	{
	  
	  logProb+=freeParticlePropagatorFunction->logEvaluate(configurations[i][0][j],configurations[i+1][0][j]);
	}
    }
  
  return logProb;
  
}



double pairParticleApproximationChain::logEvaluatePairLinks(pairParticleApproximationChain::configurations_t & configurations)
{
  double logProb=0;
  // single particle 
  for(int i=0;i<configurations.size()-1;i++)
    {
      for(int j=0;j<configurations[i][0].size();j++)
	{
	  for(int k=0;k<j;k++)
	    {
	      logProb+=pairParticlePropagatorFunction->logEvaluate(configurations[i][0][j],configurations[i][0][k],configurations[i+1][0][j],configurations[i+1][0][k]);
	    }
	}
    }
  
  return logProb;
}

double pairParticleApproximationChain::logEvaluate(pairParticleApproximationChain::configurations_t & configurations)
{
  return logEvaluateSingleLinks(configurations) + logEvaluatePairLinks(configurations) + logEvaluateTail(configurations) + logEvaluateHead(configurations) ;
  
}

double pairParticleApproximationChain::logEvaluateSingleLinks(int iParticle,int iTimeSliceBegin,int iTimeSliceEnd,configurations_t & configurations)
{
  
  double logProb=0;
  
  for(int i=iTimeSliceBegin;i<iTimeSliceEnd;i++)
    {
      logProb+=freeParticlePropagatorFunction->logEvaluate(configurations[i][0][iParticle],configurations[i+1][0][iParticle]);
    }
  
  return logProb;
}

double pairParticleApproximationChain::logEvaluatePairLinks(int iParticle,int iTimeSliceBegin,int iTimeSliceEnd,configurations_t & configurations)
{
  double logProb=0;
  
  for(int i=iTimeSliceBegin;i<iTimeSliceEnd;i++)
    {
      for(int k=0;k<configurations[i][0].size();k++)
	{
	  if (iParticle!=k)
	    {
	      logProb+=pairParticlePropagatorFunction->logEvaluate(configurations[i][0][iParticle],configurations[i][0][k],configurations[i+1][0][iParticle],configurations[i+1][0][k]);
	    }
	}
    }
  
  return logProb;
}

double pairParticleApproximationChain::logEvaluatePairLinks(int iParticle,int jParticle,int iTimeSliceBegin,int iTimeSliceEnd,configurations_t & configurations)
{
  double logProb=0;
  
  for(int i=iTimeSliceBegin;i<iTimeSliceEnd;i++)
    {
      if (iParticle != jParticle)
	{
	  logProb+=pairParticlePropagatorFunction->logEvaluate(configurations[i][0][jParticle],configurations[i][0][iParticle],configurations[i+1][0][jParticle],configurations[i+1][0][iParticle]);
	}
      
      for(int k=0;k<configurations[i][0].size();k++)
	{
	  if (iParticle!=k & jParticle!=k)
	    {
	      logProb+=pairParticlePropagatorFunction->logEvaluate(configurations[i][0][iParticle],configurations[i][0][k],configurations[i+1][0][iParticle],configurations[i+1][0][k]);

	      logProb+=pairParticlePropagatorFunction->logEvaluate(configurations[i][0][jParticle],configurations[i][0][k],configurations[i+1][0][jParticle],configurations[i+1][0][k]);
	    }
	}
      
    }
  
  return logProb;
}

double pairParticleApproximationChain::logEvaluateSingleLinks(int iParticle,int jParticle,int iTimeSliceBegin,int iTimeSliceEnd,configurations_t & configurations)
{
  return logEvaluateSingleLinks(iParticle,iTimeSliceBegin,iTimeSliceEnd,configurations) + logEvaluateSingleLinks(jParticle,iTimeSliceBegin,iTimeSliceEnd,configurations) ;
}


double pairParticleApproximationChain::logEvaluateTimeDerivative(configurations_t & configurations)
{
  double a;
  double dis=0;
  double e=0;
  
  for(int i=0;i<configurations.size()-1;i++)
    {
      for(int j=0;j<configurations[i][0].size()-1;j++)
	{
	  
	  dis=geo->distance_pbc(configurations[i][0][j].position(),configurations[i+1][0][j].position());
	  
	  e+=freeParticlePropagatorFunction->logEvaluateTimeDerivative(dis);
	}
    }
  
  int N=configurations[0][0].size();
  
  double dis2=0;
  for(int i=0;i<N;i++)
    {
      
    for(int j=0;j<i;j++)
      {
	 for(int k=0;k<configurations.size()-1;k++)
	   {	      
	     e+=pairParticlePropagatorFunction->logEvaluateTimeDerivative(configurations[k][0][i],configurations[k][0][j],configurations[k+1][0][i],configurations[k+1][0][j]);  
	   }
      }
    }
  
  return e;
}
