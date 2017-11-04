#ifndef TYPELIST_H
#define TYPELIST_H

#define TYPELIST_1(T1) typeList<T1,TL::empty_t>
#define TYPELIST_2(T1,T2) typeList< T1, TYPELIST_1(T2) >

template<class T,class U>
class typeList
  {
    typedef T head;
    typedef U tail;
  };

template<int i>
class intToType
{
  enum { value_t = i};
};

namespace TL
{
  class empty_t {};
  
  template<class Tlist,int i>
  class getType
  {
    typedef typename getType<typename Tlist::tail,i-1>::value_t value_t;
  };
  
  template<class Tlist>
  class getType<Tlist,0>
  {
    typedef typename Tlist::head value_t;
  };
}


#endif
