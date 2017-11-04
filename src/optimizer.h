
#include <complex>
template<class T>
class optimizer : T
{
  
};

// a linear policy for the hamiltonian
template<int nP>
class linearOptimizer
{
public:
  
  linearOptimizer()
  {
    // allocate the H and S matrix
    int i;
    int d;
    d=nP+1;
    ipiv.resize(d);
    H.resize(d);
    S.resize(d);
    for(i=0;i<d;i++)
      {
	H[i].resize(d);
	S[i].resize(d);
      }
    
  }
  
  // estimates the step to follow
  double estimateStep(vector<double> &HSV)
  {
    int n ;
    n=2;
    // builds the H matrix
    H[0][0]=HSV[0];
    H[0][1]=HSV[1]-HSV[2]*HSV[0];
    H[1][0]=HSV[1]-HSV[2]*HSV[0] + HSV[3];
    H[1][1]=HSV[4]-2*HSV[2]*HSV[1]
      + HSV[2]*HSV[2]*HSV[0]
      + HSV[5]-HSV[2]*HSV[3];
    // builds the S matrix
    S[0][0]=1;
    S[0][1]=0;
    S[1][0]=0;
    S[1][1]=HSV[6]
      - HSV[2]*HSV[2];
    
    
    //print_matrix(H);
    cout << "S" << S[1][1]<<endl;
    // builds S^-1H
    H[0][0]=H[0][0]*S[1][1];
    H[0][1]=H[0][1]*S[1][1];
    print_matrix(H);
    
    // finds the eigenvalues

    double b;
    double c;
    double landa;
    double delta;
    b=-(H[0][0]+H[1][1]);
    c=-H[0][1]*H[1][0];
    // lowest eigen value
    delta=b*b-4*c;
    if (delta<0)
      {
	delta=0;
      }
    
    landa=(-b - sqrt(delta))/(2.);
    cout << ((landa-H[0][0])/(H[0][1]))<<endl;
    return (landa-H[0][0])/(H[0][1]);
    
  }
  
private:
  
  vector<int> ipiv;
  vector<vector<double> > H;
  vector<vector<double> > S;
  
  
};

// implements the newton method for optimization
template<int np>
class NewtonOptimizer
{
public:
  
  double estimateStep(vector<double> &HSV)
  {
    double a,b,c,d,g;
    double step;
    // gradient
    g=2*(HSV[1]-HSV[0]*HSV[2]);
    
    a=2*(HSV[8]-HSV[7]*HSV[0] - HSV[4]+HSV[6]*HSV[0]);
    
    b=4*(HSV[4]-HSV[6]*HSV[0]-HSV[2]*g);
    d=2*(HSV[5]-HSV[2]*HSV[3]);
    cout << "d: "<<d << endl;
    cout << "g: "<<g << endl;
    cout << "a:"<<a<<endl;
    
    //cout << a << endl;
    //cout << b << endl;
    //cout << d << endl;
    //cout << g << endl;
    //cout << g << endl;
    
    step=-g/(a+b+d);
    
    //cout << step<<endl;
    //cout << a+b+d<<endl;
    return step;
  }
  
};
