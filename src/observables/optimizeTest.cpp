#include "optimizeTools.h"


int main(int argc,char** argv)
{
  vector<double> testAccums;
  testAccums.resize(7);
  vector<double> step;
  step.resize(1);
  step[0]=0;
  testAccums[0]=-161.549;
  testAccums[1]=-32948.8;
  testAccums[2]=8.95492;
  testAccums[3]=26251.6;
  testAccums[4]=5.34871e+06;
  testAccums[5]=-1078.34;
  testAccums[6]=204.167;
  
  print(testAccums);
  linearMethodStepEstimator est;
  est.buildMatrix(testAccums,1);
  
  est.print();
  est.getStep(step);
  print(step);
}
