#include "../measures.h"

void timeDifferenceSquaresAverage(const vector<double> & positions,vector<double> & timeCorrelations)
{
  
  int n=positions.size();
  assert(timeCorrelations.size()==positions.size());
  
  for(int i=0;i<n;++i)
    {
      timeCorrelations[i]=0;
      for(int j=i;j<n;j++)
	{
	  timeCorrelations[i]+=pow(positions[j]-positions[j-i],2);
	}
      timeCorrelations[i]/=(n-i);
    }
  
}

void timeProductAverage(const vector<double> & positions,vector<double> & timeCorrelations)
{
  
  int n=positions.size();
  assert(timeCorrelations.size()==positions.size());
  
  for(int i=0;i<n;++i)
    {
      timeCorrelations[i]=0;
      for(int j=i;j<n;j++)
	{
	  timeCorrelations[i]+=positions[j]*positions[j-i];
	}
      timeCorrelations[i]/=(n-i);
    }
  
}


