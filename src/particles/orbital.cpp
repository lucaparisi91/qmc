#include "orbital.h"
#include "../ptools.h"

int spinOrbital::get_pack_size() const
{
  return pTools::get_pack_size(spinValue) + pTools::get_pack_size(x) + pTools::get_pack_size(xBC);
}

int orbital1D::get_pack_size() const
{
  return pTools::get_pack_size(x) + pTools::get_pack_size(xBC);
}

void spinOrbital::pack(packed_data* packO)
{
  packO->pack(spinValue);
  packO->pack(x);
  packO->pack(xBC);
}

void orbital1D::pack(packed_data* packO)
{
  
  packO->pack(x);
  packO->pack(xBC);
}

void spinOrbital::unpack(packed_data* packO)
{
  packO->unpack(spinValue);
  packO->unpack(x);
  packO->unpack(xBC);
}

void orbital1D::unpack(packed_data* packO)
{
  packO->unpack(x);
  packO->unpack(xBC);
}

spinOrbital& spinOrbital::operator=( const spinOrbital &s)
{
  spinValue=s.spinValue;
  x=s.x;
  xBC=s.xBC;
}

orbital1D& orbital1D::operator=( const orbital1D &s)
{
  x=s.x;
  xBC=s.xBC;
}

ostream& operator<<(ostream& output ,const spinOrbital &orbital2)
{
  
  output << orbital2.spinValue<< " ";
  output << orbital2.x << " ";
  output << orbital2.xBC;
  return output;
}

ostream& operator<<(ostream& output ,const orbital1D &orbital2)
{
  output << orbital2.x << " ";
  output << orbital2.xBC;
  return output;
}

istream& operator>>(istream& input ,spinOrbital &orbital2)
{
  
  input >> orbital2.spinValue;
  input >> orbital2.x;
  input >> orbital2.xBC;
  return input;
}


istream& operator>>(istream& input ,orbital1D &orbital2)
{
  input >> orbital2.x;
  input >> orbital2.xBC;
  return input;
}

template<class T>
ostream& operator<<(ostream& output,const orbitals<T> & orbitals2)
{
  int nOrbitals;
  
  output<<"Orbitals: "<<"'"<<orbitals2.getLabel()<<"'"<<endl;
  output<<"N orbitals: "<<orbitals2.size()<<endl;
  
  for(int i=0;i<orbitals2.size();i++)
    {
      output << orbitals2.orbitalsStorage[i];
      output << endl;
    }
  return output;
  }

  
template<class T>
 istream& operator>>(istream& input,orbitals<T> & orbitals2)
  {
    string dummy;
    int N;
    input>>dummy;
    input>>dummy;
    if (dummy=="") dummy="'Unkown'";
    dummy=trim(dummy);
    orbitals2.label=dummy.substr(1,dummy.size()-2);
    input>>dummy;
    input>>dummy;
    input>>N;
    orbitals2.resize(N);
    for(int i=0;i<orbitals2.size();i++)
      {
	input >> orbitals2.orbitalsStorage[i];
      }
    return input;
  }


void buildAllOrbitals(string filename,orbitals<orbitals<orbital1D> > & allOrbitals)
{
  xml_input xml_main_input;
  orbitals<orbital1D > orbitals;
  int n;
  
  xml_main_input.open(filename)->reset();

  xml_main_input.reset()->get_child("system")->get_child("particles");
  
  do
    {
      xml_main_input.get_attribute("n");
      n=xml_main_input.get_int();
      orbitals.resize(n);
      allOrbitals.push_back(orbitals);
      
      xml_main_input.get_next("particles");
    }
  
  while(xml_main_input.check());
  
  xml_main_input.reset();
  
}

template class orbitals<orbital1D>;
template class orbitals<orbitals<orbital1D> >;
template istream& operator>>(istream& input,orbitals<orbital1D> & orbitals2);
template istream& operator>>(istream& input,orbitals<orbitals<orbital1D> > & orbitals2);
template ostream& operator<<(ostream& output,const orbitals<orbital1D> & orbitals2);
template ostream& operator<<(ostream& output,const orbitals<orbitals<orbital1D> > & orbitals2);

void generateUniform(orbitals<orbital1D > & orbitals,double xMin,double xMax)
{
  double step;
  step=(xMax-xMin)/orbitals.size();
  for(int i=0;i<orbitals.size();i++)
    {
      orbitals[i].positionNoBC()=xMin+i*step;
      
    }
}

void generateUniform(orbitals<orbitals<orbital1D> > & orbitals,double xMin,double xMax)
{
  
  for(int i=0;i<orbitals.size();i++)
    {
      generateUniform(orbitals[i],xMin,xMax);
    }
  
}



template<class randomGenerator_t>
void generateRandom(orbitals<orbital1D> & orbitals,double xMin,double xMax,randomGenerator_t* randO)
{
  
  for(int i=0;i<orbitals.size();i++)
    {
      orbitals[i].positionNoBC()=randO->rand()*(xMax-xMin) + xMin;
    }
}
