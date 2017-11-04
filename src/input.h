#include <string>
using namespace std;

// load from an input file a variable by name
string load_parameter(string filename,string param_name); //loads a parameter by reading through a list
bool check_file_exists(string filename);
bool load_integer_paramater(string filename,string name,int &value);
  bool load_real_paramater(string filename,string name,double &value);  

  
