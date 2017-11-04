#include<iostream>
#include <string>
#include<fstream>
#include"input.h"
#include"tools.h"
#include <cstdlib>
using namespace std;

string load_parameter(string filename,string param_name)
{
  ifstream input_file;
  string line="";
  string value="";
  int pos;
  int in_line;
  
  param_name=trim(param_name);
  input_file.open(filename.c_str());
   if (! input_file.good()){
     cout << "Cannot open the file.";
   }
   
  
    while(getline(input_file,line))
      {
	//cout << line<<endl;
	// avoid line beginning with a '#'
	if (line[0] == '#')
	  {
	    continue;
	  }
	
	pos=line.find("=");
       
      
	if (pos !=string::npos)
	  {
	    if (in_line==1)
	      {
		return value;
	      }
	    if (line.substr(0,pos).compare(param_name)==0 )
	      {
	  
		// found the parameter
		if ((pos+1) < line.length() && line[pos+1] == '>')
		  {
		    in_line=1;	   
		    }
		else
		  {
		    return line.substr(pos+1,line.length()-pos-1);
		  }
	
	      }
	  }
	else
	  {
	    if (in_line==1)
	  {
	    value=value+line+"\n";
	  }
	    
	  }
      }
    
    
    if (in_line ==1)
      {
	return value;
      }
    else
      {
    return "";
      }
    
}

bool check_file_exists(string filename)
{
  // checks whatever a file exist
  ifstream f(filename.c_str());
  return f.good();
}

bool load_real_paramater(string filename,string name,double &value)
{
  string tmp;
  tmp=load_parameter(filename,name);
  if (tmp.compare(string("")) == 0)
    {
      return false;
    }
  else
    {
      value=atof(tmp.c_str());
      return true;
      
    }
}
bool load_integer_paramater(string filename,string name,int &value)
{
  string tmp;
  tmp=load_parameter(filename,name);
  if (tmp.compare(string("")) == 0)
    {
      return false;
    }
  else
    {
      value=atoi(tmp.c_str());
      return true;
      
    }
}
