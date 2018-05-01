#include "pigs.h"
#include "pigsMover.h"
#include "pigsTools.h"

pigsMover::pigsMover()
{
  
}

void pigsMover::move(configurations_t & configurations)
{
  int i=sampler->sample();
  movesToAttempt[i]->move(configurations,*propagator);
  
}

string pigsMover::info()
{
  string res="Move\tratio\n";
  for(int i=0;i<movesToAttempt.size();i++)
    {
      res+=movesToAttempt[i]->getName()+"\t" + toString(movesToAttempt[i]->getAcceptanceRatio()*100) + "%\n";
    }
  
  return res;
  
}

void pigsMover::addMove(pigsMove & moveToAdd,double weight)
{
  movesWeights.push_back(weight);
  movesToAttempt.push_back(&moveToAdd);
  sampler->add(weight);
  
}

void levyReconstructor::save(int iParticle,int iTimeSliceBegin,int iTimeSliceEnd,const configurations_t & configurations)
{
  for(int i=iTimeSliceBegin;i<=iTimeSliceEnd;i++)
    {
      
      (*configurationsBackup)[i][0][iParticle].positionNoBC()=configurations[i][0][iParticle].positionNoBC();
    }
}

void levyReconstructor::recover(int iParticle,int iTimeSliceBegin,int iTimeSliceEnd,configurations_t & configurations)
{
  for(int i=iTimeSliceBegin;i<=iTimeSliceEnd;i++)
    {
      configurations[i][0][iParticle].positionNoBC()=(*configurationsBackup)[i][0][iParticle].positionNoBC();
      configurations[i][0][iParticle].positionBC()=geo->pbc(configurations[i][0][iParticle].positionNoBC());
    }  
}

void levyReconstructor::construct(int iParticle,int iTimeSliceBegin,int iTimeSliceEnd,configurations_t & configurations)
{
  double a,b,sigma;
  pos_t newPosition;
  pos_t prevPosition;
  pos_t nextPosition;
  pos_t newDistance;
  int j,jLast;
  randO->gaussian(gaussianVariables,iTimeSliceEnd-iTimeSliceBegin);
  double length=iTimeSliceEnd-iTimeSliceBegin;
  
  for(int i=iTimeSliceBegin+1;i<iTimeSliceEnd;i++)
    {
      j=i-iTimeSliceBegin;
      jLast=iTimeSliceEnd-iTimeSliceBegin;
      a=1./(jLast - j + 1);
      
      b=(jLast-j)*1./(jLast-j+1);
      sigma=sqrt(timeStep*b);
      
      newPosition=configurations[i-1][0][iParticle].positionNoBC()*b + configurations[iTimeSliceEnd][0][iParticle].positionNoBC()* a;
      
      newPosition+=gaussianVariables[j-1]*sigma;
      configurations[i][0][iParticle].positionNoBC()=newPosition;
      configurations[i][0][iParticle].positionBC()=geo->pbc(configurations[i][0][iParticle].positionNoBC());
    }
}

void wiggleMove::move(configurations_t & configurations,qmcDriver_t::propagator_t & prop)
{
  bool accept;
  double logRatio=0;
  int tmp;
  
  //select particle and imaginary time slice
  int set=0;
  iParticle=int(getRandomGenerator()->uniform()*configurations[0][set].size());
  i1=int(getRandomGenerator()->uniform()*(configurations.size()));
  i2=int(getRandomGenerator()->uniform()*(configurations.size()));
  //i1=1;
  //i2=configurations.size()-1;
  
  if (i1>i2) {tmp=i1;i1=i2;i2=tmp;}
  
  reconstructor->save(iParticle,i1,i2,configurations);
  logRatio-=prop.logEvaluatePairLinks(iParticle,i1,i2,configurations);
  reconstructor->construct(iParticle,i1,i2,configurations);
  logRatio+=prop.logEvaluatePairLinks(iParticle,i1,i2,configurations);
  
  accept=logMetropolis(logRatio,getRandomGenerator());
  
  
  if (accept)
    {
      recordSuccessfullMove();
    }
  else
    {
      reconstructor->recover(iParticle,i1,i2,configurations);
      recordUnSuccessfullMove();
    }  
}

void moveVariationalTails::move(configurations_t & configurations,propagator_t & prop)
{
  bool accept;
  int set=0;
  int iParticle;
  bool head;
  double logRatio;
  int i1,i2;
  pos_t newPosition,oldPosition;
  iParticle=int(getRandomGenerator()->uniform()*configurations[0][set].size());
  
  head=getRandomGenerator()->uniform()>=0.5 ? true : false;
  getRandomGenerator()->gaussian(gaussianRandomNumbers,2);
  
  if (head)
    {
      i2=configurations.size()-1;
      i1=i2-1;
      logRatio-=prop.logEvaluateHead(configurations);
      logRatio-=prop.logEvaluatePairLinks(iParticle,i1,i2,configurations);
      oldPosition=configurations[i2][0][iParticle].positionNoBC();
      configurations[i2][0][iParticle].positionNoBC()=gaussianRandomNumbers[0]*sqrt(timeStep) + configurations[i1][0][iParticle].positionNoBC();
      configurations[i2][0][iParticle].positionBC()=geo->pbc(configurations[i2][0][iParticle].positionNoBC());
      
      logRatio+=prop.logEvaluateHead(configurations);
      logRatio+=prop.logEvaluatePairLinks(iParticle,i1,i2,configurations);
      
      accept=metropolis(logRatio,getRandomGenerator());
      
      if (accept)
	{
	  recordSuccessfullMove();
	}
      else
	{
	  recordUnSuccessfullMove();
	  configurations[i2][0][iParticle].positionNoBC()=oldPosition;
	  configurations[i2][0][iParticle].positionBC()=geo->pbc(configurations[i2][0][iParticle].positionNoBC());
	}
    }
  else
    {
      i1=0;
      i2=i1+1;
      
      logRatio-=prop.logEvaluateTail(configurations);
      logRatio-=prop.logEvaluatePairLinks(iParticle,i1,i2,configurations);
      oldPosition=configurations[i1][0][iParticle].positionNoBC();

      
      configurations[i1][0][iParticle].positionNoBC()=gaussianRandomNumbers[0]*sqrt(timeStep) + configurations[i2][0][iParticle].positionNoBC();
      configurations[i1][0][iParticle].positionBC()=geo->pbc(configurations[i1][0][iParticle].positionNoBC());
      
      logRatio+=prop.logEvaluateTail(configurations);
      logRatio+=prop.logEvaluatePairLinks(iParticle,i1,i2,configurations);
      
      accept=metropolis(logRatio,getRandomGenerator());
      
      if (accept)
	{
	  recordSuccessfullMove();
	}
      else
	{
	  recordUnSuccessfullMove();
	  
	  configurations[i1][0][iParticle].positionNoBC()=oldPosition;
	  configurations[i1][0][iParticle].positionBC()=geo->pbc(configurations[i1][0][iParticle].positionNoBC());
	}
    }
  
}

void translateMove::move(configurations_t & configurations,qmcDriver_t::propagator_t & prop)
{
  bool accept;
  double logRatio=0;
  int i1,i2;
  
  //select particle and imaginary time slice
  int set=0;
  int iParticle=int(getRandomGenerator()->uniform()*configurations[0][set].size());
  
  i1=0;
  i2=configurations.size()-1;

  logRatio-=prop.logEvaluateHead(configurations);
  logRatio-=prop.logEvaluateTail(configurations);
  
  for(iParticle=0;iParticle<configurations[0][0].size();iParticle++)
    {
      double dis=displacement*2*(getRandomGenerator()->uniform()-0.5);
      reconstructor->save(iParticle,i1,i2,configurations);
      logRatio-=prop.logEvaluatePairLinks(iParticle,i1,i2,configurations);
      
      for(int i=0;i<configurations.size();i++)
	{
	  configurations[i][0][iParticle].positionNoBC()+=dis;
	  configurations[i][0][iParticle].positionBC()=reconstructor->getGeometry()->pbc(configurations[i][0][iParticle].positionNoBC());
	}
      logRatio+=prop.logEvaluatePairLinks(iParticle,i1,i2,configurations);
    }
  
  logRatio+=prop.logEvaluateHead(configurations);
  logRatio+=prop.logEvaluateTail(configurations);
  
  
  accept=logMetropolis(logRatio,getRandomGenerator());
  
  if (accept)
    {
      recordSuccessfullMove();
    }
  else
    {
      
      for(iParticle=0;iParticle<configurations[0][0].size();iParticle++)
	{
	  reconstructor->recover(iParticle,i1,i2,configurations);	  
	}
      recordUnSuccessfullMove();
    }  
}

void moveSingleBead::move(configurations_t & configurations,propagator_t & prop)
{
  int iParticle;
  double delta;
  double logRatio=0;
  int iTimeSlice;
  int set=0;
  int i1,i2;
  
  iParticle=int(getRandomGenerator()->uniform()*configurations[0][set].size());
  iTimeSlice=int(getRandomGenerator()->uniform()*configurations.size());  
  
  if(iTimeSlice==0)
    {
      i1=0;
      i2=1;
      logRatio-=prop.logEvaluateTail(configurations);
    }
  else if (iTimeSlice==(configurations.size()-1))
    {
      i2=configurations.size()-1;
      i1=i2-1;
      logRatio-=prop.logEvaluateHead(configurations);
    }
  else
    {
      i1=iTimeSlice-1;
      i2=iTimeSlice+1; 
    }
  
  logRatio-=prop.logEvaluatePairLinks(iParticle,i1,i2,configurations);
  
  logRatio-=prop.logEvaluateSingleLinks(iParticle,i1,i2,configurations);

  
  delta=(displacement*getRandomGenerator()->uniform()-0.5)*displacement;
  
  configurations[iTimeSlice][0][iParticle].positionNoBC()+=delta;
  configurations[iTimeSlice][0][iParticle].positionBC()=geo->pbc(configurations[iTimeSlice][0][iParticle].positionNoBC());
  
  logRatio+=prop.logEvaluatePairLinks(iParticle,i1,i2,configurations);
  logRatio+=prop.logEvaluateSingleLinks(iParticle,i1,i2,configurations);
  
  if(i2==(configurations.size()-1))
    {
      logRatio+=prop.logEvaluateHead(configurations);
    }
  else if(i1==0)
    {
      logRatio+=prop.logEvaluateTail(configurations);
    }
  
  bool accept=logMetropolis(logRatio,getRandomGenerator());
  
  if(accept)
    {
      recordSuccessfullMove();
    }
  else
    {
      configurations[iTimeSlice][0][iParticle].positionNoBC()-=delta;
      configurations[iTimeSlice][0][iParticle].positionBC()=geo->pbc(configurations[iTimeSlice][0][iParticle].positionNoBC());
      recordUnSuccessfullMove();
    }
  
}
