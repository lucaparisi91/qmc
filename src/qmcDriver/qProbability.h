#include <vector>
#include <complex>
using namespace std;

template<class grad_t,class all_particles_t>
double Q_probability(grad_t & drift_force_new,grad_t & drift_force_old,all_particles_t & position,all_particles_t & position2,double delta_tau)
{
  
  double q=0;
  
  for(int i=0;i<drift_force_old.size();i++)
    {
      for(int j=0;j<drift_force_old[i].size();j++)
	{
	  q=q+(drift_force_old[i][j] + drift_force_new[i][j])*(
							       (drift_force_old[i][j] - drift_force_new[i][j])*(1/2.)*delta_tau - (position[i][j].positionNoBC() - position2[i][j].positionNoBC())
							       );
	}
    }
  
  return q;
}

double Q_probability(vector<complex<double> > &drift_force_new,vector<complex<double> > &drift_force_old,vector<double> &position,vector<double> &position2,double delta_tau);

double Q_probability(vector<double> &drift_force_new,vector<double> &drift_force_old,vector<double> &position,vector<double> &position2,double delta_tau);

