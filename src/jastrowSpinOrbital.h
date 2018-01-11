#ifndef JASTROWSPINORBITAL_H
#define JASTROWSPINORBITAL_H

#include <string>

template<class jastrow1_t,class jastrow2_t>
class jastrowSpinOrbitalTwoBody
{
public:
  jastrowSpinOrbitalTwoBody(string filenameSame,string filenameDiff) : jastrowSame(filenameSame),jastrowDifferent(filenameDiff){};
  
  jastrowSpinOrbitalTwoBody(string file) :jastrowSame(file),jastrowDifferent(file){};
  
  inline double d0(const double x,int i,int j)
  {
    return i==j ? jastrowSame.d0(x) : jastrowDifferent.d0(x);
  }
  
  inline double d1(const double x,int i,int j)
  {
    return i==j ? jastrowSame.d1(x) : jastrowDifferent.d1(x);
  }

  inline double d2(const double x,int i,int j)
  {
    return i==j ? jastrowSame.d2(x) : jastrowDifferent.d2(x);
  }
  
  void print(string name_0d,string name_1d,string name_2d)
  {
    jastrowSame.print(name_0d + "Same.dat",name_1d + "Same.dat",name_2d+"Same.dat");
    jastrowDifferent.print(name_0d + "Different.dat",name_1d + "Differnt.dat",name_2d + "Different.dat");
  }

  void setParameter(double x,int i)
  {
    if (i==0) jastrowSame.setParameter(x,0);
    if (i==1) jastrowDifferent.setParameter(x,0);    
  }
  
  double  getParameter(int i) const
  {
    if (i==0) return jastrowSame.getParameter(0);
    if (i==1) return jastrowDifferent.getParameter(1);
    
  }
private:
  jastrow1_t jastrowSame;
  jastrow2_t jastrowDifferent;
  
};


template<class jastrow1_t,class jastrow2_t>
class jastrowSpinOrbitalOneBody
{
public:
  
  jastrowSpinOrbitalOneBody(string filenameUp,string filenameDown) : jastrowUp(filenameUp),jastrowDown(filenameDown){center[0]=jastrowUp.center;center[1]=jastrowDown.center;};
  
  jastrowSpinOrbitalOneBody(string file) :jastrowUp(file),jastrowDown(file){center[0]=jastrowUp.center;center[1]=jastrowDown.center;};  
  
  inline double d0(const double x,int i)
  {
     return i==1 ? jastrowUp.d0(x) : jastrowDown.d0(x);
  }
  
  inline double d1(const double x,int i)
  {
    return i==1 ? jastrowUp.d1(x) : jastrowDown.d1(x);
  }
  
  double d1d0(const double x,int i)
  {
    return i==1 ? jastrowUp.d1d0(x) : jastrowDown.d1d0(x);
  }

  double d2d0(const double x,int i)
  {
    return i==1 ? jastrowUp.d2d0(x) : jastrowDown.d2d0(x);
  }
  
  inline double d2(const double x,int i)
  {
    return i==1 ? jastrowUp.d2(x) : jastrowDown.d2(x);
  }
  
  

  void print(string name_0d,string name_1d,string name_2d)
  {
    jastrowUp.print(name_0d + "Up.dat",name_1d + "Up.dat",name_2d+"Up.dat");
    jastrowDown.print(name_0d + "Down.dat",name_1d + "Down.dat",name_2d + "Down.dat");
  }
  
  void setParameter(double x,int i)
  {
    if (i==0) jastrowUp.setParameter(x,0);
    if (i==1) jastrowDown.setParameter(x,0);    
  }
  
  double  getParameter(int i) const
  {
    if (i==0) return jastrowUp.getParameter(0);
    if (i==1) return jastrowDown.getParameter(1);
    
  }
  
  double getCenter(int s) const
  {
    return center[(1-s)/2];
  }
  
private:
  
  jastrow1_t jastrowUp;
  jastrow2_t jastrowDown;
  double center[2];
  
};



#endif
