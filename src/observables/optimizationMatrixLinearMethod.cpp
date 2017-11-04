/*
  Buils the H matrix as a column major vector array.
 */

#include "../tools.h"
#include "optimizationMatrixLinearMethod.h"
#include "lapacke.h"


optimizationMethodLinear::optimizationMethodLinear()
 {
   N=2;
   H=new double[N*N];
   S=new double[N*N];
   work=new double[8*N];
   alphaR=new double[N];
   alphaI=new double[N];
   beta=new double [N];
   leftEigenVectors=new double[N*N];
   rightEigenVectors=new double[N*N];
}

void optimizationMethodLinear::buildHMatrix(vector<double> &obs)
{
  H[0]=obs[0];
  H[1]=obs[1]-obs[2]*obs[0];
  H[2]=H[1] + obs[3];
  H[3]=obs[4]-2*obs[2]*obs[1] + obs[2]*obs[2]*obs[0] + obs[5]-obs[2]*obs[3];
}

void optimizationMethodLinear::buildSMatrix(vector<double> & obs)
{
  S[0]=1;
  S[1]=0;
  S[2]=0;
  S[3]=obs[6]-obs[2]*obs[2];
  
}

/*
  - solves for the generalized eigen value problem to obtain the next parameter
  
*/

double optimizationMethodLinear::getParameterVariation(vector<double> &obs)
{
  print_vector(obs);
  buildHMatrix(obs);
  buildSMatrix(obs);
  cout << "S"<<endl;
  printColumnArray(S,2,2);
  cout << "H"<<endl;
  printColumnArray(H,2,2);
  return getParameterVariation();
  
}

double optimizationMethodLinear::getParameterVariation()
{
  // compute right eigen vectors and eigen-values
  lapack_int info;
  int iMin,i;
  double alpha,alpha2;

  // solve the generalized eigen value problem
  info=LAPACKE_dggev(LAPACK_COL_MAJOR,'N','V',N,H,N,S,N,alphaR,alphaI,beta,leftEigenVectors,N,rightEigenVectors,N);
  if(info != 0)
    {
      cout << "Error in computeing generalized eigenvectors"<<endl;
      exit(3);
    }

  // determine the minimum eigenvalue
  if (beta[0]==0)
    {
      alpha=0;
    }
  else
    {
      alpha=alphaR[0]/beta[0];
    }
  
  iMin=0;
  for(i=0;i<N;i++)
    {
      if (beta[i]==0)
	{
	  alpha2=0;
	}
      else
	{
	  alpha2=alphaR[i]/beta[i];
	}
	
      if(alpha>alpha2)
	{
	  alpha=alpha2;
	  iMin=i;
	}
      
    }
  
  cout << "alphaR"<<endl;
  printColumnArray(alphaR,1,2);
  
  cout << "alphaI"<<endl;
  printColumnArray(alphaI,1,2);

  cout << "beta"<<endl;
  printColumnArray(beta,1,2);
  
  cout << "rightEigenVector"<<endl;
  printColumnArray(rightEigenVectors,2,2);
  
  
  // avoid an infinite parameter variation
  if (rightEigenVectors[iMin*N]==0)
    {
      return 0;
    }
  else
    {
      // returns v_i[1]/v_i[0]
      return rightEigenVectors[iMin*N+1]/rightEigenVectors[iMin*N];
    }
  
}

