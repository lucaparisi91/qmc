#include "configurations.h"

ostream& operator<<(ostream& out,configurationsPIGS & configurations)
{
  out<<"Beads "<< configurations.size()<<endl;
  for(int i=0;i<configurations.configurations.size();i++)
    {
      
      out<<configurations[i];
    }
}
