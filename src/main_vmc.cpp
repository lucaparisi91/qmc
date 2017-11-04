#include<iostream>
#include "vmc.h"
#include "signal_handling.h"
#include "csignal"

using namespace std;

int main()
{
  signal(SIGINT, signal_handler);
  signal(SIGTERM, signal_handler);
  vmc* sim;
  
  sim=new vmc();
  sim->run();
}



  
