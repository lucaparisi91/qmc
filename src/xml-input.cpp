#include <iostream>
#include <cstdlib>
#include "xml-input.h"
#include <string>

using namespace std;

xml_input* xml_input::open(string filename)
{
  // it parses an empty file in a system
  doc=xmlParseFile(filename.c_str());
  
  if (doc==NULL)
    {
      return NULL;
    }
  
  // reset the wavefunction of the system
  reset();
  
  return this;

  
}
// creates an empty xml document
xml_input* xml_input::new_doc(string filename="root")
{
  doc = xmlNewDoc(BAD_CAST "1.0");
  cur = xmlNewNode(NULL, BAD_CAST filename.c_str());
  xmlDocSetRootElement(doc, cur);
  return this;
}

xml_input* xml_input::save(string filename)
{
  xmlSaveFormatFile (  filename.c_str(), doc, 1);
  return this;
}
// adds a child to the document
xml_input* xml_input::add_child(string name,string value)
{
  cur=xmlNewTextChild (cur, NULL, (const xmlChar* )name.c_str(), (const xmlChar* )value.c_str());
  return this;
}
// reset the pointer to the initial value
xml_input* xml_input::reset()
{
  cur = xmlDocGetRootElement(doc);
  return this;
}
// get the next element in the system

xml_input* xml_input::get_next()
{
  cur=cur->next;
  return this;
}
xml_input* xml_input::get_next(string filename)
{
  xmlNodePtr cur2;
  cur2=cur;
  cur2=cur2->next;
  while(cur2!=NULL)
    {
      
          if ( 
	  xmlStrcmp(cur2->name, (const xmlChar *)(filename.c_str())) == 0
	  )
	    {
	      cur=cur2;
	      return this;
	    }
	  cur2=cur2->next;
	  
  
    }
  cur=NULL;
  return NULL;
}
xml_input* xml_input::get_next(string filename,string attribute,string attribute_value)
{
  xmlNodePtr cur2;
  cur2=cur;
  
  while (get_next(filename))
    {
      
      if (get_attribute(attribute) !=NULL)
	{
	  if (get_string() == attribute_value)
	    {
	      return this;
	    }
	}
    }
  cur=cur2;
  return NULL;
}

// get the number of children with filename filename
int xml_input::get_n_child(string filename)
{
  xmlNodePtr cur2;
  int i;
  cur2=cur;
  i=0;
  while( cur2!=NULL)
    {
      cur2=cur2->next;
      if (
	  xmlStrcmp(cur->name, (const xmlChar *)(filename.c_str())) == 0
	  )
	{
	  i=i+1;
	}
    }
  
  return i;
  
}
bool xml_input::check()
{
  if (cur == NULL)
    {
      return false;
    }
  else
    {
      return true;
    }
}
xml_input* xml_input::get_first_child()
{
  if (cur == NULL)
    {
      cout<<"No node";
      exit(2);
    }
  cur=cur->xmlChildrenNode;
  return this;
}
// returns the name of the current node
string xml_input::get_name()
{
  const char *name;
  name=(const char*)cur->name;
  return string(name);
}
xml_input* xml_input::get_child(string child_name)
{
  
  if (cur == NULL)
    {
      cout<<"No node";
      exit(2);
    }
  cur=cur->xmlChildrenNode;
  
  while (cur != NULL)
    {
      if ( xmlStrcmp(cur->name, (const xmlChar *)child_name.c_str()) == 0 )
      {
	return this;
      }
		 
    cur = cur->next;
    }
  cur=NULL;
  return this;
}

xml_input* xml_input::get_child(string child_name,string label,string value)
{
  
  if (cur == NULL)
    {
      cout<<"No node";
      exit(2);
    }
  cur=cur->xmlChildrenNode;
  
  while (cur != NULL)
    {
      if ( xmlStrcmp(cur->name, (const xmlChar *)child_name.c_str()) == 0 )
      {
	if (get_attribute("label")->get_string() == value)
	  {
	    return this;
	  }
      }
		 
    cur = cur->next;
    }
  cur=NULL;
  return this;
}
// returns the value for the string
xml_input* xml_input::get_value()
{
  value=xmlNodeListGetString(doc, cur->xmlChildrenNode, 1); 
  return this;
}
xml_input* xml_input::get_attribute(string name)
{
  value=xmlGetProp(cur,(const xmlChar* )name.c_str());
  if (value==NULL)
    {
      return NULL;
    }
  return this;
}

bool xml_input::get_bool()
{
  if (get_string()=="true")
    {
      return true;
    }
  else
    {
      if(get_string()=="false")
	{
	  return false;
	}
      else
	{
	  cout<<"Is not a bool"<<endl;
	  exit(2);
	}
    }
}

string xml_input::get_string()
{
  if (value==NULL)
    {
      return string("");
    }
  else
    {
      return string( reinterpret_cast<const char*>(value));
    }
}
int xml_input::get_int()
{
  return atoi(get_string().c_str());
}
double xml_input::get_real()
{
  return atof(get_string().c_str());
}

complex<double> xml_input::get_complex()
{
  complex<double> a;
  stringstream ss;
  ss << get_string();
  ss >> a;
  return a;
}

xml_input* xml_input::close()
{
  xmlFreeDoc(doc);
  return this;
}

xml_input* xml_input::set(string value)
{
  
  if (cur == NULL)
    {
      cout<<"Nothing to set."<<endl;
      exit(2);
    }
  else
    {
      xmlNodeSetContent(cur, (xmlChar*)value.c_str());
    }
  return this;
}
xml_input* xml_input::set(string attr,string value)
{
  
  if (cur == NULL)
    {
      cout<<"Nothing to set."<<endl;
      exit(2);
    }
  else
    {
      xmlSetProp (cur, 
		  (xmlChar *)attr.c_str(), 
                  (xmlChar *) value.c_str());
    }
  return this;
}

xml_input::~xml_input()
{
  //xmlFreeNode(cur);
  xmlFreeDoc(doc); 
}

string xml_input::toString()
{
  xmlChar *s;
  int size;
  string out;
  xmlBufferPtr buffer;
  buffer= xmlBufferCreate();
  size = xmlNodeDump(buffer, doc, cur, 0, 1);

  out=(char *)buffer->content;
  xmlBufferFree(buffer);
  return out;
}
