#include <cstdlib>
// use the quadratic brente algorithm to find the root of a continous function between a and b

template<class T>
double findRootBrente(T f,double a,double b,double eps)
{
  
  if (f(a)*f(b)>=0)
    {
      cout << "Not alternating ends!";
      exit(0);
      return a;
    }
  
  double c,s,tmp,d;
  bool flag;
  if (abs(f(a)) < abs(f(b)) )
    {
      // swap a and b
      tmp=a;
      a=b;
      b=tmp;
    }
  
  c=a;
  s=a;
  d=a;
  
  flag=true;
  
  while( ! (f(b)==0 or f(s)==0 or abs(b-a)<=eps))
    {
      if (  (f(a) != f(b)) and (f(b)!=f(c)))
	{
	  // use quadratic interpolation to find the closest root
	  s=a*f(b)*f(c)/( (f(a) - f(b))*(f(a)-f(c)) ) + b*f(a)*f(c)/((f(b)-f(a))*( f(b)-f(c) )) + c*f(a)*f(b)/( (f(c)-f(a))*( f(c)-f(b) ) );
	}
      else
	{
	  // metodo della secante
	  s=b-f(b)*(b-a)/(f(b)-f(a));
	}
      if (
	  (s < (3*a+b)/4 or s>b) or
	  (flag and (abs(s-b) >= abs(b-c)/2)) or
	  (! flag and (abs(s-b)>=abs(c-d)/2)  ) or
	  (  flag and (abs(b-c) < eps  ) ) or
	  ( !flag and (abs(c-d) < eps)  ) 
	  )
	{
	  s=(a+b)/2.;
	  flag=true;
	}
      else
	{
	  flag=false;
	}
      d=c;
      c=b;
      if (f(a)*f(s)<0)
	{
	  b=s;
	}
      else
	{
	  a=s;
	}

      if (abs(f(a))<abs(f(b)))
	{
	  // swap a and b
	  tmp=a;
	  a=b;
	  b=tmp;
	}
	
    }
  
  if (f(s)==0)
    {
      return s;
    }
  else
    {
      return b;
    }
}

template<class T>
void findFirstAlternating(T f,double & a,double & b,int bins=10)
{
  double step;
  double amin,bmax;
  double fa,fb;

  amin=a;
  bmax=b;
  

  while(true)
    {
      step=(bmax-amin)/bins;
      b=a;
  
      while ( b< bmax)
	{
	  fa=f(a);
	  fb=f(b);
      
	  if ( fa*fb<=0 )
	    {
	      return;
	    }
	  b=b + step;
	}
      bins=bins*2;
    }
}
