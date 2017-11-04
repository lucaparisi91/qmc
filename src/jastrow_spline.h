#include "spline.h"
#include "tools.h"

class jastrow_spline : public jastrow<jastrow_spline>
{
public:
  typedef double position_t;
  typedef double value_t;
  inline double dP1(const double &x){return 0;};
  inline double dP2(const double &x){return 0;};
  jastrow_spline(vector<double> &X,vector<double> &Y,double center=0)
  {
    set_spline(X,Y,center=center);
  }
  void set_spline(vector<double> &X,vector<double> &Y,double center=0)
  {
    assert(X.size() > 0);
    s.set_points(X,Y);
    
    this->a=X[0]; // lower value of X
    this->b=X[X.size()-1]; // higher value of X
    this->center=center;
  }
  inline double d0(const double &x) const { return s(x) ;} ;
  inline double d1(const double &x) const {return s.d1(x);};
  inline double d2(const double &x) const {return s.d2(x);};
  
  jastrow_spline(string filename )
  {
    
    vector<vector <double> > M;
    M=readMatrixFromFile(filename);
    //for(int i=0;i<M[0].size();i++)
    //  {
    //	cout << M[0][i] << " " <<M[1][i]<<endl;
    //  }
    
    set_spline(M[0],M[1]);
    
  }
  
  template<class T> jastrow_spline(T const & jastrowOr, int bins)
  {
    vector<double> X;
    vector<double> Y;
    double step;
    double x;
    int i;
    assert(bins > 1);

    this->center=jastrowOr.center;
    step= (b - a)/(bins-1);
    
    for(i=0;i<bins;i++)
      {
	x=i*step + jastrowOr->a;
	X.push_back(x);
	Y.push_back(jastrowOr->d0(x) );
      }
    
    jastrow_spline(X,Y,center=jastrowOr.center);
    
  }
private:
  
  tk::spline s;

};

template<class J>
class jastrow_spline_from_jastrow : jastrow_spline
{
public:
  jastrow_spline_from_jastrow(string filename) : jastrowInit(filename){jastrow_spline::jastrow_spline(jastrowInit);};
private:
  
  J jastrowInit;
};
