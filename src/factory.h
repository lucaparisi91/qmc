#ifndef ABSTRACT_FACTORY_H
#define ABSTRACT_FACTORY_H

#include <map>
#include <iostream>
#include <stdlib.h>
#include <string>
#include "xml-input.h"
#include "exception_classes.h"

using namespace std;

// the wavefactory class

template<class comp>
class waveFactory
{
 public:
  
  typedef typename comp::swave_t object;
  typedef object* (*funcCreatorType) ( comp*,xml_input *,const string & );
  typedef string idType;
  typedef map<idType,funcCreatorType> creatorMap;
  
  
  // create the wavefunction
  object* create(const idType &id,comp* qmc_obj,xml_input* xI, const string &filename)
  {
    typename creatorMap::const_iterator i;
    i=creators.find(id);
    if (i!= creators.end())
      {
	return (i->second)(qmc_obj,xI,filename);
      }
    else
      {
	throw factoryIdNotRecorded(id);
      }
  };
  
  bool registerType(const idType &id, funcCreatorType creator )
  {
    return creators.insert(typename creatorMap::value_type(id,creator) ).second;
  };
  
 private:
  
  creatorMap creators;
};


#endif
