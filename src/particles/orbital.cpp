#include "orbital.h"
#include "../ptools.h"

int spinOrbital1D::get_pack_size() const
{
  return pTools::get_pack_size(spinValue) + orbital1D::get_pack_size();
  
}

int orbital1D::get_pack_size() const
{
  return pTools::get_pack_size(x) + pTools::get_pack_size(xBC);
}

void spinOrbital1D::pack(packed_data* packO)
{
  orbital1D::pack(packO);
  packO->pack(spinValue);
}

void orbital1D::pack(packed_data* packO)
{
  
  packO->pack(x);
  packO->pack(xBC);
}

void spinOrbital1D::unpack(packed_data* packO)
{
  orbital1D::unpack(packO);
  packO->unpack(spinValue);
}

void orbital1D::unpack(packed_data* packO)
{
  packO->unpack(x);
  packO->unpack(xBC);
}
ostream& operator<<(ostream& out    ,const orbital1D &orbital2)
{ 
  return orbital2.output(out);
}

ostream& operator<<(ostream& out    ,const spinOrbital1D &orbital2)
{ 
  return orbital2.output(out);
}

istream& operator>>(istream& in ,orbital1D &orbital2)
{
  return orbital2.input(in);
}

istream& operator<<(istream& in ,spinOrbital1D &orbital2)
{
  return orbital2.input(in);
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

template class orbitals<orbital1D>;
template class orbitals<orbitals<orbital1D> >;
template istream& operator>>(istream& input,orbitals<orbital1D> & orbitals2);
template istream& operator>>(istream& input,orbitals<orbitals<orbital1D> > & orbitals2);
template ostream& operator<<(ostream& output,const orbitals<orbital1D> & orbitals2);
template ostream& operator<<(ostream& output,const orbitals<orbitals<orbital1D> > & orbitals2);

template class orbitals<spinOrbital1D>;
template class orbitals<orbitals<spinOrbital1D> >;
template istream& operator>>(istream& input,orbitals<spinOrbital1D> & orbitals2);
template istream& operator>>(istream& input,orbitals<orbitals<spinOrbital1D> > & orbitals2);
template ostream& operator<<(ostream& output,const orbitals<spinOrbital1D> & orbitals2);
template ostream& operator<<(ostream& output,const orbitals<orbitals<spinOrbital1D> > & orbitals2);


void setNDown(orbitals<spinOrbital1D> &orbitals,int n)
{
  assert(n<=orbitals.size());
  
  for(int i=0;i<n;i++)
    {
      orbitals[i].spin()=-1;
    } 
}

void setNDown(orbitals<orbitals<spinOrbital1D> > &orbitals,int n)
{
  
  for(int i=0;i<orbitals.size();i++)
    {
      setNDown(orbitals[i],n);
    }
}



template<>
void buildAllOrbitals(string filename,orbitals<orbitals<spinOrbital1D> > & allOrbitals)
{
  xml_input xml_main_input;
  orbitals<spinOrbital1D > orbitals;
  int n;
  int nDown;
  xml_main_input.open(filename)->reset();

  xml_main_input.reset()->get_child("system")->get_child("particles");
  
  do
    {
      xml_main_input.get_attribute("n");
      n=xml_main_input.get_int();
      xml_main_input.get_attribute("nUp");
      nDown=n-xml_main_input.get_int();
      orbitals.resize(n);
      setNDown(orbitals,nDown);
      allOrbitals.push_back(orbitals);
      
      xml_main_input.get_next("particles");
    }
  
  while(xml_main_input.check());
  
  xml_main_input.reset();
  
}

int getMagnetization(orbitals<spinOrbital1D> &p)
{
  int magnetization;
  magnetization=0;
  for(int i=0;i<p.size();i++)
    {
      magnetization+=p[i].spin();
    }
  return magnetization;
  
}

int getMagnetization(orbitals<orbital1D> &p)
{
  cout << "Magnetization of a non spin orbit. Makes no sense."<<endl;
  exit(0);
  
}
