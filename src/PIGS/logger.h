#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
using namespace std;

class singletonLogger
{
public:
  template<class T>
  singletonLogger& operator<<(T & data){logStream<<data;logStream<<endl;};
  static  singletonLogger& getInstance()
  {
    static singletonLogger singletonIstance;
    return singletonIstance;
  }
  
private:
  static ofstream logStream;

  ~singletonLogger(){singletonLogger::logStream.close();}
  singletonLogger(){singletonLogger::logStream.open("pigs.log");};


  
};

#endif
