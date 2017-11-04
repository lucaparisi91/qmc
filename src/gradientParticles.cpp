#include "gradientParticles.h"
#include "../ptools.h"

allParticlesGradient1D::allParticlesGradient1D(const vector<int> & ns)
{
  gradient.resize(ns.size());
  resize(ns);
}

void allParticlesGradient1D::clone( allParticlesGradient1D & g2)
{
  gradient=g2.gradient;
}

void allParticlesGradient1D::pack(packed_data* pack)
{
  for(int i=0;i<gradient.size();i++)
    {
      pack->pack(gradient[i]);
    }
  
}

void allParticlesGradient1D::unpack(packed_data* pack)
{
  for(int i=0;i<gradient.size();i++)
    {
      pack->unpack(gradient[i]);
    }
  
}

int allParticlesGradient1D::getPackSize()
{
  int sizePack;
  sizePack=0;
  
  for(int i=0;i<gradient.size();i++)
    {
      sizePack+=pTools::get_pack_size(gradient[i]);
    }
  return sizePack;
  
}

void allParticlesGradient1D::resize(const vector<int> &ns)
{
  for(int i=0;i<ns.size();i++)
    {
      gradient.resize(ns[i]);
    }  
}
