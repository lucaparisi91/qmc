#include "timer.h"
#include "mpi.h"
void timer::start()
{
  start_time=MPI_Wtime();
}
void timer::stop()
{
  stop_time=MPI_Wtime();
  total_time=total_time + (stop_time - start_time);
}

void timer::reset()
{
  start_time=0;
  stop_time=0;
  total_time=0;
}
double timer::get_time()
{
  return stop_time -start_time;
}

double timer::get_total_time()
{
  return total_time;
}

timer::timer(string label_)
{
  start_time=0;
  stop_time=0;
  total_time=0;
  label=label_;
}

void timer::setLabel(string label_)
{
  label=label_;
}

string timer::getLabel()
{
  return label;
}



