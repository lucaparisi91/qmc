#include "accumulator.h"
#include <iostream>
#include "mpi.h"

int main(int argc,char** argv)
{
  int task;
  vectorAccumulatorVariance<double> testAccumulator(4);
  vector<double> vec;
  vec.resize(4);
  vec[0]=1.;
  vec[1]=2.;
  vec[2]=3.;
  vec[3]=4.;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&task);
  
  testAccumulator.accumulate(vec);
  testAccumulator.accumulate(vec);
  testAccumulator.accumulate(vec);
  
  testAccumulator.transfer(0);
  
  if (task==0)
    {
      testAccumulator.getMean(vec);
      
      for(int i=0;i<vec.size();i++)
	{
	  cout << vec[i]<< " ";
	}
      cout<<endl;
    }
  
  MPI_Finalize();
}

