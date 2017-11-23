#include "spinOrbital.h"
#include "../ptools.h"

int spinOrbital::getPackSize() const
{
  return pTools::get_pack_size(spinValue) + pTools::get_pack_size(x) + pTools::get_pack_size(xBC);
}

void spinOrbital::pack(packed_data* packO)
{
  packO->pack(spinValue);
  packO->pack(x);
  packO->pack(xBC);
}

void spinOrbital::unpack(packed_data* packO)
{
  packO->unpack(spinValue);
  packO->unpack(x);
  packO->unpack(xBC);
}

spinOrbital& spinOrbital::operator=( const spinOrbital &s)
{
  spinValue=s.spinValue;
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

istream& operator>>(istream& input ,spinOrbital &orbital2)
{
  
  input >> orbital2.spinValue;
  input >> orbital2.x;
  input >> orbital2.xBC;
  return input;
}

template<class T>
ostream& operator<<(ostream& output,const orbitals<T> & orbitals2)
  {
    output<<"Orbitals: "<<"'"<<orbitals2.getLabel()<<"'"<<endl;
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
    input>>dummy;
    input>>dummy;
    if (dummy=="") dummy="'Unkown'";
    dummy=trim(dummy);
    orbitals2.label=dummy.substr(1,dummy.size()-2);
    for(int i=0;i<orbitals2.size();i++)
      {
	input >> orbitals2.orbitalsStorage[i];
      }
    return input;
  }

template class orbitals<spinOrbital>;
template class orbitals<orbitals<spinOrbital> >;
template istream& operator>>(istream& input,orbitals<spinOrbital> & orbitals2);
template istream& operator>>(istream& input,orbitals<orbitals<spinOrbital> > & orbitals2);
template ostream& operator<<(ostream& output,const orbitals<spinOrbital> & orbitals2);
template ostream& operator<<(ostream& output,const orbitals<orbitals<spinOrbital> > & orbitals2);

