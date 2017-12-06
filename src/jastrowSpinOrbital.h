#ifndef JASTROWSPINORBITAL_H
#define JASTROWSPINORBITAL_H

#include <string>

template<class jastrow1_t,class jastrow2_t>
class jastrowSpinOrbitalSymmetricalDelta
{
public:
  jastrowSpinOrbitalSymmetricalDelta(string filenameSame,string filenameDiff) : jastrowSame(filenameSame),jastrowDifferent(filenameDiff){};
  
  double d0(const double x,int i,int j)
  {
    return (1+i*j)/2 * jastrowSame.d0(x) + (1-i*j)/2 * jastrowDifferent.d0(x);
  }
  
  double d1(const double x,int i,int j)
  {
    return (1+i*j)/2 * jastrowSame.d1(x) + (1-i*j)/2 * jastrowDifferent.d1(x);
  }

  double d2(const double x,int i,int j)
  {
    return (1+i*j)/2 * jastrowSame.d2(x) + (1-i*j)/2 * jastrowDifferent.d2(x);
  }

  void print(string name_0d,string name_1d,string name_2d)
  {
    jastrowSame.print(name_0d + "Same.dat",name_1d + "Same.dat",name_2d+"Same.dat");
    jastrowDifferent.print(name_0d + "Different.dat",name_1d + "Differnt.dat",name_2d + "Different.dat");
  }

  void setParameter(double x,int i)
  {
    if (i==0) jastrowSame.setParameter(0);
    if (i==1) jastrowDifferent.setParameter(0);    
  }

  double  getParameter(int i)
  {
    if (i==0) return jastrowSame.getParameter(0);
    if (i==1) return jastrowDifferent.getParameter(1);
    
  }
private:
jastrow1_t jastrowSame;
  jastrow2_t jastrowDifferent;

}
  
#endif
