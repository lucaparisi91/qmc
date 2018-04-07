#!/home/luca/software/anaconda/bin/python
from math import *
import scipy as split
import numpy as np
from scipy import optimize
import tools
import  xml.etree.ElementTree as ET
import matplotlib.pylab as plt

class jastrow:

    def __init__(self):
        # jastrow parameters
        self.parameters={}
        self.bins=1000
        self.a=0
        self.b=50
        self.output_parameters={}
    def plot_f_root(self,a=None,b=None,bins=None):
        if a==None:
            a=self.a
        if b==None:
            b=self.b
        if bins==None:
            bins=self.bins
            
        a=float(a)
        b=float(b)
        
        bins=float(self.bins)
        x=np.arange(a,b,(b-a)/bins)
        g=np.vectorize(self.f_root)
        y=g(x)
        plt.scatter(x,y,marker="o")
    
    def plot_jastrow(self,a=None,b=None,bins=None):
        if a==None:
            a=self.a
        if b==None:
            b=self.b
        if bins==None:
            bins=self.bins
            
        a=float(a)
        b=float(b)
        
        bins=float(self.bins)
        x=np.arange(a,b,(b-a)/bins)
        g=np.vectorize(self.jastrow)
        y=g(x)
        plt.plot(x,y)
    
    def plot_jastrow_1d(self,a=None,b=None,bins=None):
        if a==None:
            a=self.a
        if b==None:
            b=self.b
        if bins==None:
            bins=self.bins
            
        a=float(a)
        b=float(b)
        
        bins=float(self.bins)
        x=np.arange(a,b,(b-a)/bins)
        g=np.vectorize(self.jastrow_1d)
        y=g(x)
        plt.plot(x,y)
        
    def plot_kinetic_eigen(self,a=None,b=None,bins=None):
        if a==None:
            a=self.a
        if b==None:
            b=self.b
        if bins==None:
            bins=self.bins
            
        a=float(a)
        b=float(b)
        
        bins=float(self.bins)
        x=np.arange(a,b,(b-a)/bins)
        g0=np.vectorize(self.jastrow)
        g2=np.vectorize(self.jastrow_2d)
        
        y=g2(x)/g0(x)
        plt.plot(x,y)
        
    def plot_jastrow_2d(self,a=None,b=None,bins=None):
        if a==None:
            a=self.a
        if b==None:
            b=self.b
        if bins==None:
            bins=self.bins
            
        a=float(a)
        b=float(b)
        
        bins=float(self.bins)
        x=np.arange(a,b,(b-a)/bins)
        g=np.vectorize(self.jastrow_2d)
        y=g(x)
        plt.plot(x,y)
        
    def print_jastrow(self,filename,a,b,bins):
        # open jastrow file
        f=open(filename,"w+")
        dx=(b-a)*1./(bins)
        f.write(str(a)+"\n")
        f.write(str(b)+"\n")
        f.write(str(bins)+"\n")
        for i in range(0,bins):
            value=self.jastrow((i+1)*dx+a)
            #if (value < 1E-7):
            #    value=1E-7
            
            f.write(str(value)+"\n")
        f.close()
    def print_jastrow_1d(self,filename,a,b,bins):
        # open jastrow file
        f=open(filename,"w+")
        dx=(b-a)*1./(bins)
        f.write(str(a)+"\n")
        f.write(str(b)+"\n")
        f.write(str(bins)+"\n")
        for i in range(0,bins):
            value=self.jastrow_1d((i+1)*dx+a)
            #if (value < 1E-7):
            #    value=1E-7
            
            f.write(str(value)+"\n")
        f.close()
    def print_jastrow_2d(self,filename,a,b,bins):
        f=open(filename,"w+")
        dx=(b-a)*1./(bins)
        f.write(str(a)+"\n")
        f.write(str(b)+"\n")
        f.write(str(bins)+"\n")
        for i in range(0,bins):
            value=self.jastrow_2d((i+1)*dx+a)
            #if (value < 1E-7):
            #    value=1E-7
            
            f.write(str(value)+"\n")
        f.close()
        
    def import_data_dict(self,d):
        for key in d:
            self.parameters[key]=float(d[key])
        return self
    def import_data_file(self,filename):
        
        d=tools.read_settings(filename)
        for key in d:
            self.parameters[key]=float(d[key])
            
    def print_jastrows(self,dirname=".",appendix="",begin=0,end=0,steps=10000):
        if (end==0 ):
            end=float(self.parameters["l_box"])
        self.print_jastrow(dirname+"/jastrow"+appendix+".dat",begin,end,steps)
        self.print_jastrow_1d(dirname+"/jastrow_1d"+appendix+".dat",begin,end,steps)
        self.print_jastrow_2d(dirname+"/jastrow_2d"+appendix+".dat",begin,end,steps)
    def print_jastrow_parameters(self,filename):
        root=ET.Element("jastrow")
        tools.dictionary_to_xml(self.parameters,root)
        
            
        tree=ET.ElementTree(root)
        
        tree.write(filename)
        
    #  makes a spline out of a certain number of points
    def make_spline(self,filename,bins=1000):
        filename=str(filename)
        l=float(self.parameters["l_box"])/2.;
        # define mesh
        X=np.arange(0,l,l/(bins-1));
        X=np.append(X,l)
        self.parameters["spline"]=filename;
        f=open(filename,"w+")
        for x in X:
            f.write( str(x) + " " + str(self.jastrow(x)) + "\n")
        f.close();
        
        
# a class to find the root of a general continuos function
    def find_root(self,step_root=10**-4,a_root=0.,b_root=10**6,eps=1e-4):
        
        lower_root=a_root
        upper_root=a_root + 2*step_root
        limit_root=b_root
                
        while (self.f_root(lower_root)*self.f_root(upper_root)>=0):
            
            if (abs(self.f_root(lower_root)) < eps) and (abs(self.f_root(upper_root)) < eps):
                return (lower_root + upper_root)/2.
            lower_root=lower_root + step_root
            upper_root=upper_root + step_root
            if upper_root > limit_root:
                print "error: Cannot find a positive root"
                exit()
       
        return optimize.brentq(self.f_root,lower_root,upper_root)


    
    def print_parameters(self):
        
        for key in self.parameters :
            print str(key) + ": " + str(self.parameters[key]) + "\n"
        
#######################delta potential##########################
class jastrow_delta(jastrow):
    
    # computes the required parameters
    def f_root(self,x) :
        return self.parameters["g"]*tan(x) - (pi/2 - x)*1./self.parameters["cut_off"]
    
    def process(self):        
        self.parameters["delta"]=self.find_root()
        
        self.parameters["k"]=self.parameters["g"]*tan(self.parameters["delta"])
        
        return self

    def jastrow(self,x):
        x=abs(x)
        if x<=self.parameters["cut_off"]:
            return sin(self.parameters["k"]*x+self.parameters["delta"])
        else:
            return self.jastrow(self.parameters["cut_off"])
    def jastrow_1d(self,x):
        x=abs(x)
        if x<=self.parameters["cut_off"]:
            return self.parameters["k"]*cos(self.parameters["k"]*x+self.parameters["delta"])
        else:
            return 0
    def jastrow_2d(self,x):
        if x<=self.parameters["cut_off"]:
            return -(self.parameters["k"]**2)*sin(self.parameters["k"]*x+self.parameters["delta"])
        else:
            return 0
        
###################### jastrow delta with gaussians#######################

class jastrow_gaussian(jastrow):
    def __init__(self,l,alpha=0.5,position=0):
        jastrow.__init__(self)
        self.parameters["alpha"]=str(alpha)
        self.parameters["position"]=str(position)
        self.parameters["l_box"]=str(l)


class jastrow_orbital(jastrow):
    def __init__(self,l,c=1,position=0):
        jastrow.__init__(self)
        self.parameters["c"]=str(c)
        self.parameters["position"]=str(position)
        self.parameters["l_box"]=str(l)

#########################3 jastrow smooth step

class jastrowSmoothStep(jastrow):
    def __init__(self,r0,p,position=0,l=100.,delta=1e-5):
        jastrow.__init__(self)
        self.parameters["r0"]=str(r0)
        self.parameters["p"]=str(p)
        self.parameters["position"]=str(position)
        self.parameters["l_box"]=str(l)
        self.parameters["delta"]=delta
        
############################ jastrow delta with phonons ########
# a delta jastrow with phononic contributions
class jastrow_delta_phonons(jastrow):
    
    def __init__(self,z,l,g):
        jastrow.__init__(self)
        self.a=0
        self.b=float(l)/2.
        self.parameters["g"]=float(g)
        self.parameters["z"]=float(z)
        self.parameters["l_box"]=float(l)
    def f_root(self,x):
        z=self.parameters["z"]
        g=self.parameters["g"]
        l_box=self.parameters["l_box"]
        # finds the root of a certain system
        delta=atan(x/g)
        beta=(x*tan((pi/l_box)*z))/( (pi/l_box)  * tan(x*z+delta))
        return (sin(x*z+delta) - sin((pi/l_box)*z)**beta)
    
    def process(self,step_root=10**-4,a_root=10**-4,b_root=100):
        self.parameters["k"]=self.find_root(step_root=step_root,a_root=a_root,b_root=10**6)
        x=self.parameters["k"]
        g=self.parameters["g"]
        z=self.parameters["z"]
        l_box=self.parameters["l_box"]
        # set the delta parameter for the simulation
        self.parameters["delta"]=atan(x/g )
        
        
        self.parameters["beta"]=(x*tan(pi/l_box*z))/(pi/l_box  * tan(x*z+self.parameters["delta"]))
        
    def jastrow(self,x):
        k=self.parameters["k"]
        z=self.parameters["z"]
        delta=self.parameters["delta"]
        l_box=self.parameters["l_box"]
        beta=self.parameters["beta"]
        y=abs(x)
        if y < z:
            return sin(k*y+delta)
        else:
            return sin(pi/l_box*y)**beta
    def jastrow_1d(self,x):
        k=self.parameters["k"]
        z=self.parameters["z"]
        delta=self.parameters["delta"]
        l_box=self.parameters["l_box"]
        beta=self.parameters["beta"]
        y=abs(x)
        if y < z:
            return cos(k*y+delta)*k
        else:
            return beta*sin(pi/l_box*y)**(beta-1)*cos(pi/l_box*y)*pi/l_box
    def jastrow_2d(self,x):
        k=self.parameters["k"]
        z=self.parameters["z"]
        delta=self.parameters["delta"]
        l_box=self.parameters["l_box"]
        beta=self.parameters["beta"]
        y=abs(x)
        if y < z:
            return -sin(k*y+delta)*(k**2)
        else:
            return beta*(beta - 1)*sin(pi/l_box*y)**(beta-2)*(cos(pi/l_box*y)*pi/l_box)**2 - beta*(pi/l_box)**2*sin(pi/l_box*y)**beta
    
        
#############################barrier########################
class jastrow_barrier(jastrow):
    
    def f_root(self,x):
        result=0.
        n=0
        k0=0.
        k1=0.
        '''
        while( k0**2/2. <= self.parameters["v0"] ):
            
            n=n+1
            k0=(n*pi/2 - x)*1./self.parameters["l_box"]      
        
        self.parameters["k0"]=k0
        k1=sqrt(k0**2-self.parameters["v0"]*2)
        self.parameters["k1"]=k1
        result=x + atan(k0/(k1*tan(k1*self.parameters["l_barrier"]))) - k0*self.parameters["l_barrier"]
        '''
        k0=(pi/2-x)/self.parameters["l_box"]
        k1=sqrt(self.parameters["v0"]*2 - k0)
        
        result= x - atan(k0/(k1*tanh(k1*self.parameters["l_barrier"]))) + k0*self.parameters["l_barrier"]
        return result
    
    def process(self):
        #self.parameters["k0"]=sqrt(self.parameters["e"]*2)
        #k0=self.parameters["k0"]
        #self.parameters["k1"]=sqrt(k0**2-self.parameters["v0"]*2)
        #k1=self.parameters["k1"]
        #self.parameters["delta"]=atan(k0/(k1*tanh(k1*self.parameters["l_barrier"]))) - k0*self.parameters["l_barrier"]
        self.parameters["delta"]=self.find_root(step_root=1e-5,a_root=-pi,b_root=pi)
        self.parameters["k0"]=(pi/2 - self.parameters["delta"])/self.parameters["l_box"]
        self.parameters["k1"]=sqrt(self.parameters["v0"]*2 -self.parameters["k0"] )
        
        
        self.parameters["A"]=sin(self.parameters["k0"]*self.parameters["l_barrier"] + self.parameters["delta"])/cosh(self.parameters["k1"]*self.parameters["l_barrier"])
        
        return self
    
    def jastrow(self,x):
        x=abs(x);
        if (x<=self.parameters["l_barrier"]):
            return self.parameters["A"]*cosh(self.parameters["k1"]*x)
        else:
            return sin(self.parameters["k0"]*x + self.parameters["delta"])
        
##################### jastrow delta in trap #########################################

class jastrow_delta_in_trap(jastrow):
    
    def __init__(self,lbox,g,m=1./2,position=0):
        jastrow.__init__(self)
        
        if g=="inf":
            c=0
        else:
            c=1/(m*float(g))
        
        self.parameters["c"]=c
        self.parameters["l_box"]=float(lbox)
        self.parameters["position"]=float(position)
        
    def jastrow(self,x):
        return abs(x) + parameters["g"]


##################### jastrow delta in trap #########################################

class jastrow_delta_bound_state_no_pbc(jastrow):
    
    def __init__(self,lbox,g,beta,alpha,m=1./2,position=0,delta=1e-5):
        jastrow.__init__(self)
        
        self.parameters["k"]=g*m
        self.parameters["l_box"]=float(lbox)
        self.parameters["position"]=float(position)
        self.parameters["beta"]=float(beta)
        self.parameters["alpha"]=float(alpha)
        self.parameters["delta"]=float(delta)
        
    def jastrow(self,x):
        return exp(-self.parameters["k"]*x)

class jastrow_delta_bound_state_no_pbc2(jastrow):
    
    def __init__(self,lbox,g,k2,xI,m=1./2,position=0):
        jastrow.__init__(self)
        
        k=g*m*1.
        
        if k2>k:
            print "Invalid k2"
            return None

        if xI< 1./(k-k2):
            b=1./(1/(k-k2)-xI)
        else:
            print "invalid xI"
            return None
        
        a=exp(-(k-k2)*xI)*(1+b*xI)
        
        #b=k*exp(-(k-k2)*x)/(1/x**2+k2/x)
        #b=1/(1/(k-k2) - x)
        
        self.parameters["k"]=g*m
        self.parameters["b"]=b
        self.parameters["a"]=a
        
        self.parameters["xI"]=xI
        self.parameters["l_box"]=float(lbox)
        self.parameters["position"]=float(position)
        self.parameters["k2"]=k2
        
    def __call__(self,xTot):
            sol=np.zeros(len(xTot))
            x=xTot[xTot<=self.parameters["xI"] ]
            
            sol[xTot<=self.parameters["xI"] ]=np.exp(-self.parameters["k"]*x)
            
            x=xTot[xTot>self.parameters["xI"] ]
            
            #sol[xTot>self.parameters["xI"] ]=self.parameters["b"]*np.exp(-self.parameters["k2"]*x)/x

            sol[xTot>self.parameters["xI"] ]=self.parameters["a"]*np.exp(-self.parameters["k2"]*x)/(1+self.parameters["b"]*x)

            

            return sol
    def d1(self,xTot):
        sol=np.zeros(len(xTot))
        x=xTot[xTot<=self.parameters["xI"] ]
            
        sol[xTot<=self.parameters["xI"] ]= -self.parameters["k"]*np.exp(-self.parameters["k"]*x)
            
        x=xTot[xTot>self.parameters["xI"] ]

        sol[xTot>self.parameters["xI"] ]=-self.parameters["a"]*np.exp(-self.parameters["k2"]*x)*(self.parameters["k2"]/(1+self.parameters["b"]*x) + self.parameters["b"]/(1+self.parameters["b"]*x)**2 )

        return sol


class jastrow_delta_bound_state_no_pbc3(jastrow):
    
    def __init__(self,lbox,g,r1,p=1,m=1./2,position=0):
        jastrow.__init__(self)
        
        self.parameters["r0"]=1/(g*m)
        self.parameters["k"]=g*m
        self.parameters["r1"]=r1
        self.parameters["p"]=p
        self.parameters["delta"]=1e-5
        
        
        self.parameters["l_box"]=float(lbox)
        self.parameters["position"]=float(position)

        
        
    def __call__(self,x):
        alpha=self.parameters["r0"]/self.parameters["r1"]
        p=self.parameters["p"]
        return np.exp(-x/self.parameters["r0"] *(alpha*x+p)/(x+p))
    
    
################### jastrow delta bound state  ###################

class jastrow_delta_bound_state(jastrow):
    # returns the root of a certain function
    # returns the root of a certain function
    def f_root(self,x):
        beta=exp(-x*self.parameters["cut_off"]*2)
        return x - self.parameters["m"]*self.parameters["g"]*(1+beta)/(1-beta)
    
    def process(self,step=1e-5,a=None,b=20*pi,log=True):
        if a is None:
            a=step
        self.parameters["k"]=self.find_root(step_root=step,a_root=a,b_root=b)
        self.parameters["beta"]=exp(-2*self.parameters["cut_off"]*self.parameters["k"])
        if log:
            print self.parameters["k"]
        
        #print self.parameters["beta"]
    def jastrow(self,x):
        x=abs(x)
        if (x<=self.parameters["cut_off"]):
            return exp(-self.parameters["k"]*x) + self.parameters["beta"]*exp(self.parameters["k"]*x)
        
def delta_bound_state(g,l_box,cut_off,m=1./2,step=1e-5):
    j=jastrow_delta_bound_state();
    j.parameters["g"]=abs(g)
    j.parameters["l_box"]=l_box
    j.parameters["cut_off"]=cut_off
    j.parameters["m"]=m
    j.process(step=step)
    return j

def TG_jastrow(n_particles,cut_off=None,l_box=None,position=0):
    j=jastrow_delta()
    if l_box == None:
        l_box=n_particles
    if cut_off == None:
        cut_off=n_particles/2.
    
    j.parameters["k"]=pi/(2*cut_off)
    j.parameters["delta"]=0
    j.parameters["l_box"]=l_box
    j.parameters["cut_off"]=cut_off
    j.parameters["g"]="inf"
    j.parameters["position"]=position;
    print "k=" + str(j.parameters["k"])
    print "delta=" + str(j.parameters["delta"])
    return j

def jastrow_general(n_particles,g,cut_off=None,l_box=None,position=0):
    if l_box == None:
        l_box=n_particles
    if cut_off == None:
        cut_off=n_particles*1./2
        
    j=jastrow_delta()
    j.parameters["g"]=g
    j.parameters["l_box"]=l_box
    j.parameters["cut_off"]=cut_off
    j.parameters["position"]=position
    j.process()
    #print "k=" + str(j.parameters["k"])
    #print "delta=" + str(j.parameters["delta"])
    return j
