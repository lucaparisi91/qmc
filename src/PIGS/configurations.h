#ifndef PIGS_CONFIGURATIONS_H
#define PIGS_CONFIGURATIONS_H


template<class configuration_t>
class configurations_t
{
public:
  inline configuration_t & operator[](int i){return configurations[i];}
  inline const configuration_t &  operator[](int i) const {return configurations[i];} 
  
private:
  vector<configuration_t> configurations;
};

#endif
