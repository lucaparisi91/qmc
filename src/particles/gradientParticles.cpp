#include "gradientParticles.h"
#include "../ptools.h"
#include <cassert>

allParticlesGradient1D::allParticlesGradient1D(const vector<int> & ns)
{  
  
  resize(ns);
}

allParticlesGradient1D::allParticlesGradient1D()
{
  
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
  gradient.resize(ns.size());
  
  for(int i=0;i<ns.size();i++)
    {
      gradient[i].resize(ns[i]);
    }
}

void allParticlesGradient1D::reset()
{
  
  for(int i=0;i<gradient.size();i++)
    {
      for(int j=0;j<gradient[i].size();j++)
	{
	  gradient[i][j]=0;
	}
    }
    
}

 allParticlesGradient1D& allParticlesGradient1D::operator+=(const allParticlesGradient1D & grad2)
{
  
  assert(grad2.gradient.size() == gradient.size() );
  
  for(int i=0;i<gradient.size();i++)
    {
      assert(grad2[i].size() == gradient[i].size() );
      
      for(int j=0;j<gradient[i].size();j++)
	{
	  gradient[i][j]+=grad2[i][j];
	}
    }
  return (*this);
}

void allParticlesGradient1D::print()
{
  
  for(int i=0;i<gradient.size();i++)
    {
      for(int j=0;j<gradient[i].size();j++)
	{
	  cout<<gradient[i][j]<<",";
	}
      cout <<endl;
    }
    
}  
