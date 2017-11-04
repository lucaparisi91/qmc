#include <string>
#include<vector>
using namespace std;

// a general function with parameters to be defined
template<class T>
class func
{
 public:
  vector<double> parameters;
  double a;
  double b;
  double d0(double x)
  {
    return static_cast<T*>(this)->d0(x);
  }
  
  double d1(double x)
  {
    return static_cast<T*>(this)->d1(x);
  }
  
  double d2(double x)
  {
    return static_cast<T*>(this)->d2(x);
  }
  void print(string filename);
};
// a function built from two other functions
template<class T>
  class double_func : func<T> 
{
 public:
  func<T>* f1;
  func<T>* f2;
  
  double_func(func<T>* f1_,func<T>* f2_);
  double d0(double x);
  double d1(double x);
  double d2(double x);
  /*
  double d1(double x)
  {
    return T1->d1(x)*T2->d0(x) + T1->d0(x)*T2->d1(x);
  }
  
  double d2(double x)
  {
    return T1->d2(x)*T2->d0(x) + T1->d1(x)*(T2->d1(x)*T2->d2(x)) + T1->d0(x)*T2->d2(x);
  }
  */
  
};
  

class delta_wavefunction : func<delta_wavefunction>
{
  delta_wavefunction(string filename);
  double d0(double distance);//evaluate the function
  double d1(double distance);//evaluates the first derivative
  double d2(double distance);//evaluates the second derivative
  void load_parameters(string filename);//load parameters from a file
};
