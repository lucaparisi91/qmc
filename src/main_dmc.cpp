#include<iostream>
#include "vmc.h"
#include "signal_handling.h"
#include "csignal"

using namespace std;

int main()
{
  signal(SIGINT, signal_handler);
  signal(SIGTERM, signal_handler);
  dmc* sim;
  
  sim=new vmc();
  sim->run();
}
