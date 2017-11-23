#include "spinOrbital.h"
#include <fstream>
int main(int argc,char** argv)
{
  
  ofstream file;
  ifstream fileI;
  spinOrbital orbital1;
  
  //--------- test spin orbital
  orbital1.position()=1.5678956;
  orbital1.spin()=1;
  orbital1.positionBC()=orbital1.position();
  orbital1.position()=3;
  file.open("save.dat");
  file.precision(17);
  file<< orbital1<<endl;
  file.close();
  orbital1.spin()=-1;
  orbital1.position()=-1.9;
  orbital1.positionBC()=-3.9;
  fileI.open("save.dat");
  fileI.precision(17);
  fileI>>orbital1;
  fileI.close();
  cout.precision(17);
  cout<<orbital1<<endl;
  
  // -- test orbitals -----
  
  orbitals<spinOrbital> configMC(10);
  configMC.setLabel("SetOne");
  for(int i=0;i<configMC.size();i++)
    {
      configMC[i]=orbital1;
      
    }
  
  file.open("save.dat");
  file << configMC;
  file.close();
  
  fileI.open("save.dat");
  fileI >> configMC;
  fileI.close();
  
  orbitals<orbitals<spinOrbital> > allConfigMC(1);
  allConfigMC.setLabel("All");
  allConfigMC[0]=configMC;

  file.open("save.dat");
  file<<allConfigMC;
  file.close();

  fileI.open("save.dat");
  fileI>>allConfigMC;
  file.close();
  
  cout<<allConfigMC[0]<<endl;
  cout<<allConfigMC.nParticles()<<endl;
}
