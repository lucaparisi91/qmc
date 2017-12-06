#include "orbital.h"
#include <fstream>

int main(int argc,char** argv)
{
  
  ofstream file;
  ifstream fileI;
  orbital1D orbital1;
  
  //--------- test spin orbital
  orbital1.position()=1.5678956;
  //orbital1.spin()=1;
  orbital1.positionBC()=orbital1.position();
  orbital1.position()=3;
  file.open("orbital.dat");
  file.precision(17);
  file<< orbital1<<endl;
  file.close();
  //orbital1.spin()=-1;
  orbital1.position()=-1.9;
  orbital1.positionBC()=-3.9;
  fileI.open("orbital.dat");
  fileI.precision(17);
  fileI>>orbital1;
  fileI.close();
  cout.precision(17);
  //cout<<orbital1<<endl;
  
  // -- test orbitals -----
  
  orbitals<orbital1D> configMC(10);
  configMC.setLabel("SetOne");
  for(int i=0;i<configMC.size();i++)
    {
      configMC[i]=orbital1;
      
    }
  
  //file.open("orbitals.dat");
  //file << configMC;
  //file.close();
  orbitals<orbital1D> configMC2;
  fileI.open("orbitals.dat");
  fileI >> configMC2;
  fileI.close();
  
  //cout<<configMC2;
  
  orbitals<orbitals<orbital1D> > allConfigMC(1);
  allConfigMC.setLabel("All");
  allConfigMC[0]=configMC;

  file.open("save.dat");
  file<<allConfigMC;
  file.close();
  orbitals<orbitals<orbital1D> > allConfigMC2;
  fileI.open("save.dat");
  fileI>>allConfigMC2;
  fileI.close();
  
  cout<<allConfigMC2<<endl;
  //cout<<allConfigMC.nParticles()<<endl;

  
}
