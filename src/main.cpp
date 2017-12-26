#include "dmc.h"
#include "traits.h"
#include "signal_handling.h"
#include "csignal"
#include "input.h"
#include "measures.h"
#include <cstdlib>
#include <sys/types.h>
#include <unistd.h>
#include "xml-input.h"
#include "omp.h"
#include "mpi.h"
#include <fenv.h>
#include "vmc.h"
using namespace std;

int main(int argc,char** argv)
{
  
  //signal(SIGINT, signal_handler);
  //signal(SIGTERM, signal_handler);
 
  //rqmc* rqmc_o;
  //gatherer* g;
  ofstream log_file;
  //g=new gatherer();
  xml_input* main_input;
  string calculation;
  
  string bcInput;

  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  try
    {
      bool isSpinor;
  main_input=new xml_input;
  
  main_input->open("input.xml");

  bcInput=get_bc(main_input);
  calculation=main_input->get_child("method")->get_attribute("kind")->get_string();
  
  if (main_input->get_attribute("spinor") )
    {
      isSpinor=main_input->get_bool();
    }
  else
    {
      isSpinor=false;
    }
  cout << "--"<<calculation<<"--"<<endl;
  log_file.open("qmc.log");
  log_file << "PID="<<getpid()<<endl;
  log_file.close();
  MPI_Init(&argc,&argv);
  
  if (calculation == "dmc" or calculation=="svmc")
    {
      // if (isSpinor)
      // 	{
      // 	  dmc<spinor1D>* dmc_o;
      // 	  dmc_o=new dmc<spinor1D>;
      // 	  dmc_o->run();
      // 	}
      // else
      // 	{

      //cout << bcInput << endl;
      if (bcInput == "periodic")
	{
	dmc< D1_t<pbc1d> >* dmc_o;
        dmc_o=new dmc<D1_t<pbc1d> >;
        dmc_o->run();
	}
      
      if (bcInput=="none")
	{
	 dmc< D1_t<noPbcD1> >* dmc_o;
         dmc_o=new dmc<D1_t<noPbcD1> >;
         dmc_o->run();
	}
	// }
    }
  else
    {
      if(calculation=="vmc")
   // 	 {
   // 	   if (isSpinor)
   // 	     {
   // 	       vmc<spinor1D>* vmc_o;
   // 	       vmc_o=new vmc<spinor1D>;
   // 	       vmc_o->run();
   // 	     }
   // 	   else
   // 	     {
	  if (bcInput == "periodic")
	{
	vmc< D1_t<pbc1d> >* vmc_o;
        vmc_o=new vmc<D1_t<pbc1d> >;
        vmc_o->run();
	}
      
      if (bcInput=="none")
	{
	 vmc< D1_t<noPbcD1> >* vmc_o;
         vmc_o=new vmc<D1_t<noPbcD1> >;
         vmc_o->run();
	}
   // 	     }
   // 	 }
      }
     
      MPI_Finalize();
    }
  
  catch(my_exception& feat)
    {
      
      feat.what();
    }
}
  


