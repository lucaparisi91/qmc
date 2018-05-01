#include "pigs.h"

int main(int argc,char **argv)
{
  MPI_Init(&argc,&argv);
  
  pigsDriver driverO;
  
  driverO.load();
  driverO.run();
  
  MPI_Finalize();
}
