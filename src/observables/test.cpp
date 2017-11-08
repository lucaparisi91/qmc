#include "accumulator.h"
#include <iostream>
int main(int argc,char** argv)
{
  vectorAccumulator<double,varianceAccumulatorVectorStrategy<double> > testAccumulator(4);
  vector<double> vec;
  vec.resize(4);
  vec[0]=1.;
  vec[1]=2.;
  vec[2]=3.;
  vec[3]=4.;

  testAccumulator.accumulate(vec);
  testAccumulator.accumulate(vec);
  testAccumulator.accumulate(vec);
  testAccumulator.getMeanSquares(vec);
  for(int i=0;i<vec.size();i++)
    {
      cout << vec[i]<< " "<<endl;
    }
  
}

