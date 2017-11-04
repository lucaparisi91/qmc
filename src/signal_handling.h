#include "globals.h"
#include <iostream>
using namespace std;
void signal_handler( int signal)
{
  cout <<"Interruption caught."<<endl;
  // a global variable to catch that the program is not running
  running=false;
  
}
