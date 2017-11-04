#ifndef TIMER_H
#define TIMER_H
#include<string>

using namespace std;

class timer
{
 public:
  timer(string label_=" ");
  double total_time;
  double start_time;
  double stop_time;
  void start();
  void stop();
  void reset();
  void setLabel(string label);
  string getLabel();
  double get_time();
  double get_total_time();
private:
  string label;
};

#endif
