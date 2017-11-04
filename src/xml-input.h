#ifndef XML_INPUT_H
#define XML_INPUT_H

#include <libxml/parser.h>
#include <libxml/xmlIO.h>
#include <libxml/xinclude.h>
#include <libxml/tree.h>
#include <string>
#include <complex>
using namespace std;

class xml_input
{
public:
  xmlDocPtr doc;
  xmlNodePtr cur;
  xmlChar* value;
  // move to the child
  xml_input* open(string filename);
  xml_input* get_child(string child_name); // find the next child with that name
  xml_input* get_child(string child_name,string attr, string value);
  xml_input* get_value();
  xml_input* new_doc(string filename);
  xml_input* get_attribute(string filename);
  xml_input* get_next();
  string toString();
  xml_input* get_next(string name);
  xml_input* get_next(string name,string attribute,string attribute_value);
  int get_n_child(string filename);
  bool check();
  string get_string(); // treats the field as a string
  int get_int(); // treats the input field as an integer number
  double get_real(); // treats the input field as a real number
  bool get_bool();
  string get_name();
  complex<double> get_complex();
  xml_input* save(string filename);
  xml_input* add_child(string name,string value);
  xml_input* open();
  xml_input* get_first_child();
  
  xml_input* reset();
  xml_input* close();
  xml_input* set(string value);
  xml_input* set(string attr,string value);
  
  ~xml_input();

  
};

#endif
