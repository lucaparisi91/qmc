#include "qProbability.h"
#include "../tools.h"

/*
double Q_probability(vector<complex<double> > &drift_force_new,vector<complex<double> > &drift_force_old,vector<double> &position,vector<double> &position2,double delta_tau)
{
  int i;
  double q;
  q=0;
  
  for (i=0;i<drift_force_old.size();i++)
    {
      q=q+(real_part(drift_force_old[i]) + real_part(drift_force_new[i])  )*(
									     ( real_part(drift_force_old[i]) - real_part(drift_force_new[i]))*(1/2.)*delta_tau - (position2[i] - position[i])
									     );
    }
  return q;
}

double Q_probability(vector<double> &drift_force_new,vector<double> &drift_force_old,vector<double> &position,vector<double> &position2,double delta_tau)
{
  int i;
  double q;
  q=0;
  
  for (i=0;i<drift_force_old.size();i++)
    {
      q=q+(drift_force_old[i] + drift_force_new[i])*(
	  (drift_force_old[i] - drift_force_new[i])*(1/2.)*delta_tau - (position2[i] - position[i])
						     );
    }
  return q;
}

*/
