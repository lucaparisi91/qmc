import pandas as pd
import scipy as sp
import numpy as np
import matplotlib.pylab as plt
import scipy.interpolate as interp
from math import *
import exceptions


class model1d:
    def __init__(self,x,y,delta=None,minX=None,maxX=None):
        
        self.minX=min(x)
        self.maxX=max(x)

        if minX==None:
            self.minX=min(x)
        else:
            self.minX=minX

        if maxX==None:
            self.maxX=max(x)
        else:
            self.maxX=self.maxX
            

        # ordering and slicing of the input data
        data=pd.DataFrame({"x":np.array(x),"y":np.array(y)})
        
        if delta is not None:
            data["delta"]=np.array(delta)
        data=data.sort_values(by="x")
        
        data=data[ (data["x"] >=  self.minX) &( data["x"]<=self.maxX)]
        
        self.x=np.array(data["x"])
        self.y=np.array(data["y"])
        
        if delta is not None:
            self.delta=np.array(delta)
        else:
            self.delta=delta
        
        
    def plot(self,bins=1000,der=0,log=False,error=False,color=None):
        if log==True:
            x=np.logspace(log10(self.minX)+1e-7,log10(self.maxX)-1e-7,num=bins)
            plt.xscale("log")
        else:
            x=np.linspace(self.minX,self.maxX,num=bins)

        y=self.__call__(x,der=der)
        plt.plot(x,y,color=color)
        
        if(error==True):
            ymin=self.__call__(x,der=der,value="top")
            ymax=self.__call__(x,der=der,value="bottom")
            
            plt.plot(x,ymax,color=color)
            plt.plot(x,ymin,color=color)
            plt.fill_between(x, ymin,ymax,alpha=0.5,color=color)
            
    def plotTraining(self):
        if self.delta is None:
            plt.plot(self.x,self.y,"o")
        else:
            print(len(self.x))
            plt.errorbar(self.x,self.y,self.delta,fmt="o")
        
    
class splineModel(model1d):

    def __init__(self,x,y,s=0.):
        
        model1d.__init__(self,x,y)
        self.spline=interp.splrep(self.x,self.y,s=s)
        
    def __call__(self,x,der=0):
        
        if ( (min(x)<self.minX) or (max(x) > self.maxX) ):
            raise exceptions.outOfRange()
        return interp.splev(x,self.spline,der=der)
    
# model with a polynomial of a certain degree
class poly1dModel(model1d):
    def __init__(self,x,y,n,minX=None,maxX=None):
        self.n=n
        model1d.__init__(self,x,y,minX=minX,maxX=maxX)
        self.pol=np.poly1d(np.polyfit(self.x,self.y,n))
        
    # calls the values of a certain wavefunction
    def __call__(self,x,der=0):
        return np.polyder(self.pol,m=der)(x)
    
class fitModel(model1d):
    
    def __init__(self,x,y,f,p0=None,minX=None,maxX=None,f1=None,f2=None,delta=None):
        model1d.__init__(self,x,y,minX=minX,maxX=maxX,delta=delta)
        #
        self.f=f
        self.f1=f1
        self.f2=f2
        
        self.fit=sp.optimize.curve_fit(f,self.x,self.y,sigma=self.delta,p0=p0)
        
        
    def __call__(self,x,der=0,value="mean"):
        
        # get the required parameters
        
        params=self.fit[0]
        paramErrors=np.sqrt(np.diag(self.fit[1]))

        # check if computer the mean value or the top and bottom values for the parameters
        if value=="top":
            params=params + paramErrors
        else:
            if value=="bottom":
                params=params-paramErrors
            else:
                if value=="mean":
                    pass
                else:
                    raise exceptions.invalideParameter()
        
        # check which derivative i have to compute
        if der==0:
            return self.f(x,*params)
        if der==1:
            if self.f1==None:
                raise exceptions.missing()
            else:
                return self.f1(x,*params)
        if der==2:
            if self.f2==None:
                raise exceptions.missing()
            else:
                return self.f2(x,*params)

        raise exceptions.invalideParameter()
    
    
    # with with the model  f(x)= a*x/(x+b)
class growthModel1(fitModel):
    
    def __init__(self,x,y,p0=None,minX=None,maxX=None,delta=None):
        def f(x,a,b):
            return a*x/(x+b)
        def f1(x,a,b):
            return a*(1/(x+b) - x/(x+b)**2)
        def f2(x,a,b):
            return a*(-2/(x+b)**2 + 2*x/(x+b)**3)
        
        fitModel.__init__(self,x,y,f,p0=None,f1=f1,f2=f2,minX=minX,maxX=maxX,delta=delta)
        
class growthModel2(fitModel):
            
    def __init__(self,x,y,p0=None,minX=None,maxX=None,delta=None):
        def f(x,a,b,c,k):
            return a*x/(x+b) + c*np.tanh(k*x)
        def f1(x,a,b,c,k):
            return ( a*(1/(x+b) - x/(x+b)**2) ) + ( c*k/np.cosh(k*x)**2   ) 
        def f2(x,a,b,c,k):
            return a*(-2/(x+b)**2 + 2*x/(x+b)**3) + c*( -2/(np.cosh(k*x)**3) *k**2 )

        fitModel.__init__(self,x,y,f,p0=None,f1=f1,f2=f2,minX=minX,maxX=maxX,delta=delta)
