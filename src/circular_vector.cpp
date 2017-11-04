#include "circular_vector.h"
#include "cassert"

circular_vector::circular_vector(int len_)
{
  
  len=len_;
  head=0;
  vec.resize(len);
}
int circular_vector::index(int i)
{
  i=(i+ head)%len;
  if (i<0)
    {
      return i+len;
    }
  else
    {
      return i;
    }
  assert(i<len);
  assert(i>=0);
}
double circular_vector::get(int i)
{
  return vec[index(i)];
}

void circular_vector::set(int i,double value)
{
  vec[index(i)]=value;
}

void circular_vector::increment(int i)
{
  vec[index(i)]=vec[index(i)]+1;
}

void circular_vector::stack(double value)
{
  increment_head();
  set(0,value);
}
void circular_vector::increment_head()
{
  head=(head+1)%len;
}
