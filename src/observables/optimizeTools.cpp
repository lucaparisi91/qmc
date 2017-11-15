#include "optimizeTools.h"
#include "lapacke.h"
#include <cassert>
#include <cstdlib>

template<class T>
rawMatrix<T>::rawMatrix()
{
  n=0;m=0;
}

template<class T>
const void rawMatrix<T>::print() const
{
  for(int i=0;i<n;i++)
    {
      for(int j=0;j<m;j++)
	{
	  printf("%f ",(*this)(i,j));
	}
      printf("\n");
    }
}


template<class T>
void rawMatrix<T>::resize(unsigned int i,unsigned int j)
{
  M.resize(i*j,0);
  n=i;m=j;
}

template<class T>
const T& rawMatrix<T>::operator()(int i,int j) const
{
  return M[index(i,j)];
}

template<class T>
T& rawMatrix<T>::operator()(int i,int j)
{
  return M[index(i,j)];
}

int getMinRealPartEigenValue(vector<double> & alphaR,vector<double> & alphaI,vector<double> &beta,double tol=1e-4)
{
  int iMin;
  double alpha,alpha2;
  assert(alphaR.size()==alphaI.size());
  assert(alphaR.size()==beta.size());
  int N=alphaR.size();

  //print(alphaR);
  //print(alphaI);
  //print(beta);
  
  if (beta[0]<=tol)
    {
      alpha=0;
    }
  else
    {
      alpha=alphaR[0]/beta[0];
    }
  
  iMin=0;
  
  for(int i=0;i<N;i++)
    {
      if (beta[i]<=tol)
	{
	  alpha2=alpha;
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
  return iMin;
}

int linearMethodStepEstimator::getStep(vector<double> &parameters)
{
  lapack_int info;
  int n;
  int status;
  double tol=1e-4;

  /* return:
     0 : OK
     i > 0 : complex i_th eigenvalue
     -i : failed on i_th eigenvalue
   */

  status=0;
  parameters.resize(nP);
  assert(H.getN()==H.getM());
  assert(S.getN()==S.getM());
  assert(S.getN()==H.getN());
  n=H.getN();
  assert(n>0);
  
  vector<double> alphaR;
  vector<double> alphaI;
  vector<double> beta;
  
  rawMatrix<double> LEIGS;
  rawMatrix<double> REIGS;
  LEIGS.resize(n,n);
  REIGS.resize(n,n);
  beta.resize(n);
  alphaR.resize(n);
  alphaI.resize(n);
  
  info=LAPACKE_dggev(LAPACK_COL_MAJOR,'N','V',n,H.raw(),n,S.raw(),n,&alphaR[0],&alphaI[0],&beta[0],LEIGS.raw(),n,REIGS.raw(),n);
  
  if(info != 0)
    {
      printf( "Error in computing generalized eigenvectors\n");
      exit(3);
    }
  
  // check if eigenvalues are real or imaginary
  
  for(int i=0;i<alphaI.size();i++)
    {
      if (alphaI[i]>tol)
	{
	  status=i;
	} 
    }
  
  //set step to minimum eigenvalue
  int iMin=getMinRealPartEigenValue(alphaR,alphaI,beta);
  #ifdef VERBOSE
  printf("Right eigenvectors\n");
  REIGS.print();
  #endif
  
  for(int j=0;j<nP;j++)
    {
      if (REIGS(0,iMin)==0)
	{
	  parameters[j]=0;
	  status=-j;
	}
      else
	{
	  parameters[j]=REIGS(j+1,iMin)/REIGS(0,iMin);
	  
	}
    }
  return status;
}

void linearMethodStepEstimator::buildMatrix(const vector<double> &accumulatedData,int nP_)
{
  int k;
  double e;
  nP=nP_;
  e=accumulatedData[accumulatedData.size()-1];
  H.resize(nP+1,nP+1);
  S.resize(nP+1,nP+1);

  // set up overlap matrix S
  S(0,0)=1;
  
  for(int j=0;j<nP;j++)
    {
      S(0,j+1)=0;
      S(j+1,0)=0;
    }
  
  k=nP*3;
  for(int i=0;i<nP;i++)
    {
      for(int j=0;j<=i;j++)
	{
	  S(i+1,j+1)=accumulatedData[k]-accumulatedData[3*i]*accumulatedData[3*j];
	  S(j+1,i+1)=S(i+1,j+1);
	  k+=2;
	}
    }
  
  //set up the H matrix
  for(int j=0;j<nP;j++)
    {
      H(j+1,0)=accumulatedData[3*j+1]-accumulatedData[3*j]*e;
      H(0,j+1)=H(j+1,0)+ accumulatedData[3*j+2];
    }
  
  k=nP*3;
  for(int i=0;i<nP;i++)
    {
      for(int j=0;j<=i;j++)
	{
	  H(i+1,j+1)=accumulatedData[k+1]-accumulatedData[3*i]*accumulatedData[3*j+1] - accumulatedData[3*j]*accumulatedData[3*i+1] + accumulatedData[3*j]*accumulatedData[3*i]*e;

	  H(j+1,i+1)=H(i+1,j+1);
	  k+=2;
	}
    }
  
  for(int i=0;i<nP;i++)
    {
      for(int j=0;j<nP;j++)
	{	  
	  H(i+1,j+1)+=accumulatedData[k]-accumulatedData[3*i]*accumulatedData[3*j+2];
	  k++;
	}
    }
  H(0,0)=e;  
}

void linearMethodStepEstimator::addDiagonal(double element)
{
  
  for(int i=1;i<nP+1;i++)
    {
      H(i,i)+=element;
    }
}

template class rawMatrix<double>;
