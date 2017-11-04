// one dimensional mesh
#ifndef MESH_H
#define MESH_H

class mesh
{
 public:
  // extrems of the mesh
  double a;
  double b;
  double step;
  // number of bins in a mesh
  int bins;
  // remeber the index given a real value
  int index(double x);
  mesh(double a_,double b_,int bins_);
  
  double get_step(){return step;}
};

#endif
