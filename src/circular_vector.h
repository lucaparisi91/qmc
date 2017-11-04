#include <vector>
using namespace std;
class circular_vector
{
 public:
  vector<double> vec;
  int len; //actual length of the vector
  int head; // position of the head
  int index(int i);
  double get(int i);// get an element of array with pbc
  void set(int i,double value);
  void increment(int i);
  void stack(double value);// add an element at the end of the vector and increse the head position by one
  void increment_head();
  
  circular_vector(int len_);
  
};
