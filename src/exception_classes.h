#ifndef EXCEPTION_CLASSES_H
#define EXCEPTION_CLASSES_H

#include <exception>

class my_exception
{
public:
virtual void what(){};
};
  
class bad_dimension : public my_exception
{
public:
virtual void what() {cout<< "Bad dimension number."<<endl;};
};

class bad_xml : public my_exception
{
public:
virtual void what() {cout<< "Bad xml file."<<endl;};
};

class get_value_null_xml_element :public my_exception
{
public:
virtual void what() {cout<< "Tring to get a value from a null xml element."<<endl;};
};

class missing_input : public my_exception
{
public:
  missing_input(string missed_ ) : missed(missed_) {}
  virtual void what(){cout<< "Missing " <<missed<< " from input file"<<endl;};
private:
  string missed;
};

class notYetSsupported : public my_exception
{
public:
notYetSsupported(string feature_): feature(feature_){};
  virtual void what(){cout<< "Feature not supported: "<<feature<<endl;}
private:
  string feature;
};

class factoryIdNotRecorded : public my_exception
{
public:
factoryIdNotRecorded(string obj_):obj(obj_){};
virtual void what(){cout << "Factory object '"<<obj <<"' not found."<<endl;}
private:
string obj;

};
class noOptimization : public my_exception
{
public:
virtual void what(){cout<<"Wrong wavefunction to optimize."<<endl;}

};

class jastrowNotFound : public my_exception
{
public:
jastrowNotFound(string obj_):obj(obj_){};
virtual void what(){cout << "Jastrow '"<<obj <<"' not found."<<endl;}
private:
string obj;

};


class unkown_parameter
{
public:
  unkown_parameter(string label_) : label(label_){}
  virtual void wwhat(){cout<<"Unkown parameter: "<<label<<endl;}
private:
  string label;
};

#endif
