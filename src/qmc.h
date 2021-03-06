#ifndef QMC
#define QMC
extern bool running;
#include <vector>
#include "mpi.h"
#include "traits.h"
#include "xml-input.h"
#include "potential.h"

template<class T> class geometry;
class random1;
class counter;
class xml_input;
class timer;
using namespace std;

template<class comp>
class qmc : comp
{
public:
  typedef typename comp::geometry_t geometry_t;
  typedef typename comp::all_particles_t all_particles_t;
  typedef potential<comp> potential_t;
  
public :
  qmc();

  
  counter *c;
  
  int mpi_tasks;
  int mpi_task;
  double t_init;
  double max_steps;
  double t_end;
  
  xml_input* xml_save;
  xml_input* main_input;
  xml_input* xml_load;
  vector<MPI_Request> walkers_send_requests;
  vector<MPI_Request> walkers_recv_requests;
  double current_step;
  vector<int> core_populations;
 
  geometry_t * geo;// describe the geometry of the system
  rand_t *rand;
  double jumps;
  double delta_tau;
  double skip;
  int warmupBlocks;
  int stepsPerBlock;
  double n_metropolis;
  double success_metropolis;
  int nBlocks;
  
  vector<timer*> timers;
  void saveGeneralQmc();
  void loadGeneralQmc();
  void startTimers();
  void stopTimers();
  void printTimers();
  void setInputFile(string inputFileName_){inputFileName=inputFileName_;}
  string getInputFileName() const {return inputFileName;}
  
  //virtual void run()=0;
  //virtual void save()=0;// saves the state of the system
  //virtual void load()=0;// loads the state of the system
  //virtual void step()=0;//performs a single state of the Quantum Monte Carlo algorithm
  //virtual void load_wavefunctions(xml_input* xml_wave);
private:
  string inputFileName;
  
};


#include "qmc.hpp"


#endif
