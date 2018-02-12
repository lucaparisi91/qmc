#ifndef PIGS_H
#define PIGS_H

#include "../random.h"


template<class comp>
class pigsDriver : public qmc<comp>
{
public:
  typedef typename comp::rand_t rand_t;
  
  
  pigsDriver(){};
  void run(){printf("PIGS run not implemented yet");};
  
private:
  
};



#endif

