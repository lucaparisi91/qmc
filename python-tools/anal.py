import scipy as sp
import scipy.stats as stats
import numpy as np
from math import *
import tools
import os
import re
import  xml.etree.ElementTree as ET
from scipy.optimize import curve_fit
import statsmodels.tsa.stattools as tsa
import matplotlib.pylab as plt
import pandas as pd
import string
# error classes
class bound_error:
    pass

class too_many_bins:
    pass
class not_found:
    pass

class not_enough_data:
    def warn():
        print "The block as zero size"

class bad_line_format():
    pass
class zero_block_size:
    def warn():
        print "The block as zero size"

class error_estimate(object) :
    def __init__(self,values):
        self.values=values
        self.errors=[]
        self.mean=np.mean(self.values)
        self.mean2=np.mean(self.values**2)
        self.error=0
        self.indices_autocorr=np.array([])
        self.indices=np.array([])
        self.autocorr=np.array([])
    def estimate(self):        
        print "Estimating not defined"
        return self
    def autocorrelate_value(self, offset):
        a=0
        n=len(self.values)-offset
        for i in range(0,n):
            a=a+self.values[i]*self.values[i+offset]
        a=a/n
        
            
        return (a - self.mean**2)/(self.mean2-self.mean**2)
    
    def autocorrelate(self,bins=100):
        
        self.tmax=len(self.values)/2
        self.tmin=0
        self.autocorr=np.zeros(bins)
        self.indices_autocorr=np.zeros(bins)
        step=(self.tmax-self.tmin)/bins
        for i in range(0,bins):
            self.autocorr[i]=self.autocorrelate_value(i*step+self.tmin)
            self.indices_autocorr[i]=i*step+self.tmin
        
        return self
    
    def plot_autocorrelate(self):
        if len(self.autocorr) <= 1:
            print "Nothing to plot"
            return self
        indices=[x for x in range(0,len(self.autocorr))]
        plt.plot(self.indices_autocorr,self.autocorr,linestyle='solid')
        plt.show()
        return self
    def get_mean(self):
        return self.mean
    def get_error(self):
        return self.errors[-1]
    def get_errors(self):
        return (self.errors,self.indices)
    
    def plot_error(self):
        if len(self.errors) <= 1:
            print "Nothing to plot"
            return self
        
        plt.plot(self.indices,self.errors,marker='o',linestyle='dashed')
        plt.show()
        return self
    def save(self,dirname="."):
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        m=np.zeros((len(self.indices),2))
        m[:,0]=self.indices
        m[:,1]=self.errors
        tools.write_matrix_in_file(m,dirname+"/errors.out")
        #tools.write_matrix_in_file([self.indices_autocorr,self.autocorr],dirname+"/autocorr.out")
        d={}
        d["mean"]=self.mean
        if (len(self.errors) > 0):
            d["error_mean"]=self.errors[-1]
        tools.write_dictionary_in_file(d,dirname+"/status.log")
            
class binning(error_estimate) :
    
    def binning(self,vec):
        #./print vec
        if (len(vec) <= 1):

            return vec[0]
        self.errors.append(vec.std()/sqrt(len(vec)))
        
        #print vec.mean()
        for i in range(0,len(vec)/2):
            vec[i]=(vec[2*i]+vec[2*i+1])/2.
    
        return self.binning(vec[0:len(vec)/2])
    def estimate(self):
        
        values_to_bin=np.copy(self.values)
        mean=self.binning(self.values)
        self.indices=[x for x in range(self.errors)]
        return self
    
class reblock(error_estimate):
    def __init__(self,values):
         super(reblock, self).__init__(values)
         
         self.min_blocks=10
         self.max_blocks=len(self.values)-1
         self.bins=10
    def reblock(self,block_size):
        if (block_size==0):
            raise zero_block_size
        n_blocks=len(self.values)/block_size
        
        mean_block=0.
        mean_blocks=0.
        mean2_blocks=0.
        for i in range(0,n_blocks):
            mean_block=0
            for j in range(0,block_size):
                mean_block=mean_block+self.values[i*block_size + j]
            mean_block=mean_block/block_size
            mean_blocks=mean_blocks+mean_block
            mean2_blocks=mean2_blocks+mean_block*mean_block
        mean_blocks=mean_blocks/n_blocks
        mean2_blocks=mean2_blocks/n_blocks

        return [sqrt((mean2_blocks - mean_blocks**2)/n_blocks),sqrt((mean2_blocks - mean_blocks**2)) ]
    
    def estimate(self):
        if (len(self.values) < self.max_blocks):
            print "Not enogh values. Decrease max_blocks or increase statistics"
            return self
        
        block_size_min=len(self.values)/self.max_blocks
        block_size_max=len(self.values)/self.min_blocks
        step=(block_size_max-block_size_min)/self.bins
        self.errors=np.zeros(self.bins)
        self.variances=np.zeros(self.bins)
        self.indices=np.zeros(self.bins)
        if block_size_min==0:
            raise too_many_bins
        for i in range(0,self.bins):
            blocksize=i*step +block_size_min
            self.errors[i]= self.reblock(blocksize)[0]
            self.variances[i]=self.reblock(blocksize)[1]
            self.indices[i]=blocksize
        return self
    
class observable(object):
    
    def __init__(self):
        self.values=np.array([])
        
    def import_data(self,values):
        self.values=values
        
    def make_binning(self):
        self.error=binning(self.values)
        self.error.estimate()
        self.error.plot_error()
        return self.error
    def reblock(self):
        self.error=reblock(self.values)
        self.error.estimate()
        return self.error
    def plot(self):
        k=[ x for x in range(0,len(self.values))]
        plt.plot(k,self.values,marker='o',linestyle='dashed')
    
def anal_energy_dmc(dirname="."):
    
    e=observable()
    e.values=tools.read_array_from_file(dirname+"/e_dmc.dat")
    e.values=e.values[jumps::]
    #e.reblock().autocorrelate().plot_autocorrelate()
    e.reblock().estimate().autocorrelate().save(dirname+"/energy_dmc")
    return e

def anal_g(dirname="."):
    data=tools.read_matrix_from_file(dirname+"/g.dat")
    values=[]
    errors=[]
    
    n=data.shape[0]
    i=0
    for column in data:
        
        e=observable()
        e.import_data(column)
        r=e.reblock()
        values.append(r.get_mean())
        errors.append(r.get_error())
        i=i+1
        print "completed " + str(i) +" of " +str(n)
        
    
    tools.write_matrix_in_file([values,errors],"g.out")

def plot_linear((a,b),(x1,x2),bins=100):
    x=np.array([k for k in np.arange(x1,x2,(x2-x1)/bins)])
    plt.plot(x,a*x+ b)

        
def plot_single(filename,begin=0,alpha=1,end=0):
    values=tools.read_array_from_file(filename,begin=begin,end=end)
    
    end=len(values)
    #print len(values)
    k=[ x for x in range(begin,end+begin)]
    #print len(k)
    plt.plot(k,values,marker='o',linestyle='none',alpha=alpha)
def plot_double(filename):
    values=tools.read_matrix_from_file(filename)
    plt.plot(values[0],values[1],marker='o',linestyle='none')
def plot_pair_correlation(filename):
    values=tools.read_matrix_from_file(filename)
    x=[k for k in range(0,len(values[0]))]
    plt.plot(x,values[0],marker='o',linestyle='none')

def extrapolate_error(values,makePlot=False):
    index=np.arange(1,len(values)+1)
    index2=1./np.sqrt(index)
    res=stats.linregress(index2,values)
    if makePlot:
        plt.scatter(index,values)
        plt.plot(index,res[1] + res[0]/np.sqrt(index))
    p=(res[0]/res[1])**2/(len(values)+1.)
    return (res,p)

def anal_winding(filename,filename2):
    
    data=mean_vector(filename);
    f=open(filename2,"w+")
    for i in range(0,len(data[0])):
        f.write(str(i) + " " + str(data[0][i]) + " " + str(data[1][i]) + "\n")
        
    f.close()



def anal_energy(dirname=".",jumps=0,nMax=None,bins=None,filename="e.dat",sub_directory="energy",e_min=None,e_max=None,makePlot=False):
    
    values=np.array(tools.read_matrix_from_file(dirname+"/"+filename,nMax=nMax)[1])
    values=values[jumps::]
    
    return anal_energy_values(values,bins=bins,dirname=dirname,sub_directory=sub_directory,makePlot=makePlot)


def anal_energy_values(values,bins=None,dirname=".",sub_directory="energy",makePlot=False):
    e=observable()
    e.values=values
    if bins==None:
        bins=max(  (10,int(sqrt(len(e.values))))  )
    #e.import_data(tools.read_array_from_file(dirname+"/"+filename))
    #e.reblock().autocorrelate().plot_autocorrelate()
    r=e.reblock()
    r.bins=bins
    r.estimate().save(dirname+"/"+sub_directory)
    res=extrapolate_error(r.errors,makePlot=makePlot)
    #res2=extrapolate_error(r.variances,makePlot=True)
    return [r.mean,res[0][1],res[1],np.var(e.values)]

def analVectorHistory(filename,directory="analVectorData",bins=None,makePlot=False,jumps=0):
    df=pd.read_csv(filename,sep=" ",index_col=False,header=None,names=["i","value","conv"])
    xs=list(set(df["i"]))
    res=[]
    for x in xs:
        df1=df[df["i"]==x]
        row=[x]
        row=row + anal_energy_values(np.array(df1["value"][jumps::]),dirname=".",sub_directory=directory,bins=bins,makePlot=makePlot)
        res.append(row)
        
    return np.array(res)
    
def anal_pair_correlation(dirname=".",jumps=0,bins=10,filename="g.dat"):
    matrix=tools.read_matrix_from_file(dirname+"/"+filename)
    values=np.zeros((matrix.shape[0],2))
    
    for j in range(0,matrix.shape[0]):
        e=observable()
        e.values=matrix[j,jumps::]
        #e.reblock().autocorrelate().plot_autocorrelate()
        r=e.reblock()
        r.bins=bins
        r.estimate().save(dirname+"/pair_correlation/"+str(j))
        values[j,0]=r.mean
        values[j,1]=r.errors[-1]
        tools.write_matrix_in_file(values,dirname+"/pair_correlation/g.dat")
    
    return None

def fit_single_file(filename,e_min=None,e_max=None,begin=None,end=None):
    values=tools.read_array_from_file(filename,e_min=e_min,e_max=e_max)
    x=[k for k in range(1,len(values)+1)]
    if begin == None:
        begin=0
    if end== None:
        end=len(x)
        
    res=stats.linregress(x[begin:end],values[begin:end])
    plot_linear((res[0],res[1]),(begin,end))
    return res

def winding_number(filename,window,bins=0):
    
    data=tools.read_matrix_from_file(filename);
    print "file read."
    print len(data[1])

    if (bins ==0):
        bins=window
    
    times=np.zeros(bins)
    values=np.zeros(bins)
    ns=np.zeros(bins)
    step=window/bins
    
    for i in range(0,bins):
        times[i]=i*step
    
    for k in range(0,len(data[1])-window,window):
    
        for j in range(k,window+k,step):
        
            for i in range(0,bins):
                
                values[i]=values[i]+(data[1][j]-data[1][j - i*step])**2
                ns[i]=ns[i]+1
    
    for i in range(0,bins):
        if (ns[i] != 0):
            values[i]=values[i]/ns[i]
        else:
            values[i]=0
    plt.scatter(times,values)
    return [times,values]
            
#anal_energy_vmc("run_1")
#plot_double("run_1/energy/errors.out")
#anal_g("run_1")
#plt.show()

def mean_vector(filename):
    i=0;
    maxI=0;
    value=0;
    vec=np.array([]);
    vec2=np.array([]);
    
    ns=np.array([])
    f=open(filename)
    # check the size of the vector
    for line in f:
        i=int(line.split(" ")[0])
        
        if (i>=maxI):
            maxI=i
            
    f.seek(0)
    
    vec=np.resize(vec,maxI+1)
    vec2=np.resize(vec2,maxI+1)
    ns=np.resize(ns,maxI+1)
    vec.fill(0)
    vec2.fill(0)
    ns.fill(0)
    
    for line in f:
        line=line.split(" ");
        i=int(line[0])
        value=float(line[1])
        vec[i]=vec[i] + value
        vec2[i]=vec2[i] + value**2 
        ns[i]=ns[i]+1;

    for i in range(0,maxI+1):
        vec[i]=vec[i]/ns[i];
        vec2[i]=sqrt((vec2[i]/ns[i] - vec[i]**2)/ns[i]) 
    f.close()

    return [vec,vec2]
    
def plot_file(filename,error=True,label="",reset_index=False,jumps=0,linear_increment=False):
    data=np.array(tools.read_matrix_from_file(filename))
    if len(data)<3:
        error=False
    
    if reset_index:
        data[0]=data[0]-data[0][0]
        
    if linear_increment:
        data[0]=range(1,len(data[0])+1)
    data=data[:,jumps:]
    
    if error :
        return plt.errorbar(data[0],data[1], yerr=data[2],label=label,marker='o',linestyle='none')
    else:
        return plt.plot(data[0],data[1],label=label,marker='o')[0]
# use this function in the fit to extract the effective mass

def effective_mass_exp_curve(x,a,b,c):
    return a + b/x*(1-np.exp(-c*x))
def effective_mass_curve_1(x,a,b):
    return a + b/x
def linear(x,a,b):
    return a*x+b
def hyperbole(x,a,b):
    return a+b/x
def effective_mass_curve_2(x,a,b,c):
    return a*x + b*(1-np.exp(-c*x))
def plot_effective_mass(data_file,settings_file="input.xml",label="",param_file="effective_mass_t.dat",factor=1,fit_lines=True):
    root=ET.parse(settings_file).getroot()
    skip=int(root.find("measures").find("winding_number").attrib["skip"])
    time_step=float(root.find("method").find("delta_tau").text)
    fit_mean=[]
    fit_low=[]
    fit_up=[]
    
    data=tools.read_matrix_from_file(data_file)
    data[1]=data[1]*factor
    data[2]=data[2]*factor
    #plt.ylim(0,1)
    if (os.path.isfile(param_file) and fit_lines==True):
        
        with open(param_file,"r") as f:
            cut_low=int(f.readline())
            cut_heigh=int(f.readline())
            params=f.readline().split();
            errors=f.readline().split();
            a=(time_step*(skip+1))
            for x in (data[0][cut_low:cut_heigh]*a):

                #values.append(effective_mass_exp_curve(float(x),float(params[0]),float(params[1]),float(params[2])))
                fit_mean.append(hyperbole(float(x),float(params[0]),float(params[1])*a))
                fit_low.append(hyperbole(float(x),float(params[0])-float(errors[0]),float(params[1])*a - float(errors[1])*a))
                fit_up.append(hyperbole(float(x),float(params[0])+float(errors[0]),float(params[1])*a + float(errors[1])*a))
            plt.plot(data[0][cut_low:cut_heigh]*a,fit_mean,linewidth=4)
            plt.fill_between(data[0][cut_low:cut_heigh]*a,fit_low,fit_up,alpha=0.4)
            
    return plt.errorbar(data[0]*time_step*(skip+1),data[1]/(data[0]*time_step*(skip+1)), yerr=data[2]/(data[0]*time_step*(skip+1)),label=label,fmt='o',ms=1)

# analize the effective mass for the system
def anal_effective_mass(data_file,settings_file="input.xml",label="",out_file="effective_mass_fit.dat",cut_low=500,cut_heigh=None,n_cuts=10,factor=1):
    root=ET.parse(settings_file).getroot()
    
    skip=int(root.find("measures").find("winding_number").attrib["skip"])
    time_step=float(root.find("method").find("delta_tau").text)
    #print "anal" + str(data_file) + "\n"
    data=tools.read_matrix_from_file(data_file)
    if (factor != 1):
        for i in range(len(data[1])):
            data[1][i]=factor*data[1][i]
            
    p0=[1,1,0.001]
    errors=[]
    #params,covs=curve_fit(effective_mass_exp_curve,data[0][700:],(data[1]/(data[0]*time_step*(skip+1)))[700:],p0=p0)
    #params,covs=curve_fit(effective_mass_curve_2,data[0][cut_low:],data[1][cut_low:]/(time_step*(skip+1)),p0=p0,maxfev=100000)
    #errors.append(covs[0][0])
    #errors.append(covs[1][1])
    #errors.append(covs[2][2])
    
    # convert time bounds to array indexes
    cut_low= int(cut_low/(time_step*(skip+1)) )
    
    if cut_heigh==None:
        cut_heigh=len(data[0])-1
    else:
        cut_heigh=int(cut_heigh/(time_step*(skip+1)) )
    # check bounds
    if cut_low<0 or cut_heigh<0:
        raise bound_error()
    else:
        if ( cut_low >= len(data[0]) or cut_heigh >= len(data[1]) ):
            
            raise bound_error()
        
    # size of the cut
    size_of_cut=(cut_heigh - cut_low)/n_cuts
    
    #mean_params=np.zeros((3,n_cuts))
    #mean2_params=np.zeros((3,n_cuts))
    #error_params=np.zeros((3,n_cuts))

    # mean values and related errors
    means=np.zeros(3)
    means2=np.zeros(3)
    errors=np.zeros(3)
    
    # in the range of the number of cuts
    for i in range(0,n_cuts):
        cut_low_current=cut_low + i*size_of_cut
        cut_heigh_current=cut_low + (i+1)*size_of_cut
        if (cut_heigh_current - cut_low_current<= 2):
            raise not_enough_data
        params,covs=curve_fit(linear,data[0][cut_low_current:cut_heigh_current],data[1][cut_low_current:cut_heigh_current]/(time_step*(skip+1)),sigma=data[2][cut_low_current:cut_heigh_current]/(time_step*(skip+1)),maxfev=100000)
        means[0]=means[0]+ params[0]
        means[1]=means[1] + params[1]
        # accumulate the mean of the parameter
        means2[0]=means2[0] + params[0]**2
        means2[1]=means2[1] + params[1]**2
        # accumulates the square of the errors
        errors[0]=errors[0] + covs[0][0]
        errors[1]=errors[0] + covs[1][1]
    means=means/n_cuts
    means2=means2/n_cuts
    errors=np.sqrt(abs((means2-means**2) + errors/n_cuts))      
    
    with open(out_file,"w") as f:
        #f.write(str(params[0])+" "+ str(params[1]) + " " + str(params[2])+"\n")
        f.write(str(cut_low)+ "\n")
        f.write(str(cut_heigh)+"\n")
        
        f.write(str(means[0])+" "+ str(means[1]) + " " + "0"+"\n")
        #f.write(str(covs[0][0]) + " " + str(covs[1][1]) + " "+ str(covs[2][2]) +"\n")
        f.write(str(errors[0]) + " " + str(errors[1]) + " "+ "0" +"\n")
        
    return (params,covs)


def plot_all(dirname=".",target="g_tilde",jumps=0,reset_index=False):
    plot_all_energy(target)
    plt.show()
    plot_all_density(dirname,target)
    plt.show()
    plot_all_effective_mass(dirname,target)
    plt.show()
    plot_all_pair_correlation(dirname,target)
    plt.show()

def get_effective_mass(filename):
        with open(filename,"r") as f:
            cut_low=int(f.readline())
            cut_heigh=int(f.readline())
            
            params=f.readline().split();
            errors=f.readline().split();
        return [params[0],errors[0]]

# check if stationarity is reached
def is_stationary(filename,bins=10,eps=0.01):
    # the values
    values=np.array(tools.read_matrix_from_file(filename)[1])
    # number of steps
    steps=np.array(tools.read_matrix_from_file(filename)[0])
    #steps=np.arange(0,len(values))
    #tab=np.zeros((len(steps),2))
    #tab[:,0]=steps 
    #tab[:,1]=values
    #print tab
    
    converged=False
    step=(len(values)-1)/bins
    
    for i in range(0,bins):
        low=i*step
        heigh=(i+1)*step
        tsa_r=tsa.adfuller(values[low:heigh],regression="c")
        if (tsa_r[0] < tsa_r[4]["5%"]):
            converged=True
        
    return converged

def interpolatePoint(x1,x2,x):
    return x1[1] + (x2[1] - x1[1])*1./(x2[0]-x1[0])*(x-x1[0])
def interpolatePoints(points,x):
    # find the range of x
    for i in range(0,len(points)-1):
        if x>=points[i][0] and x<points[i+1][0]:
            return interpolatePoint(points[i],points[i+1],x)
    if x>= points[-1][0]:
        return points[-1][1]
    if x<points[0][0]:
        return points[0][1]
    return None
class not_enough_data:
    pass

# extract the effective mass from the system
def extractEffectiveMass(df,cut_low=0,cut_right=None,n_cuts=5,makePlot=False):
    if cut_right==None:
        cut_right=max(df["t"])
    window=(cut_right-cut_low)/n_cuts
    
    slope=np.zeros(n_cuts)
    intercept=np.zeros(n_cuts)
    slopeError=np.zeros(n_cuts)
    interceptError=np.zeros(n_cuts)
    
    for i in range(0,n_cuts):
        # select the right intervals for the linear fit
        cut_low_current=cut_low + i*window
        cut_high_current=cut_low + (i+1)*window
        df1=df[(df["t"]>cut_low_current) & (df["t"]<cut_high_current) ]
        if len(df1["t"]) <= 3:
            raise not_enough_data()
        
        params,covs=curve_fit(linear,df1["t"],df1["W"],sigma=df1["deltaW"],maxfev=100000)
        slope[i]=params[0]
        intercept[i]=params[1]
        slopeError[i]=sqrt(covs[0][0])
        interceptError[i]=sqrt(covs[1][1])
        
        if makePlot==True:
            up=(slope[i]+slopeError[i])+(intercept[i]+interceptError[i])/df1["t"]
            down=(slope[i]-slopeError[i])+(intercept[i]-interceptError[i])/df1["t"]
            
            plt.fill_between(df1["t"],up,down,alpha=0.4)
            plt.errorbar(np.array(df1["t"]),np.array(df1["W"])/np.array(df1["t"]),np.array(df1["deltaW"])/np.array(df1["t"]),fmt="or")
    return np.array([slope.mean(),sqrt( slopeError.mean()**2 + slope.var())])

# returns the minimum of a polynomial in a certain range
def getMinPol(p,xmin,xmax):

    # compute the roots of  first derivative
    r=p.deriv().r
    r=r[r.imag==0].real
    # select points with second derivative positive
    r=r[p.deriv(2)(r)> 0]

    # select only points within the boundary
    
    r=r[r<xmax]
    r=r[r>xmin]

    return r


def getOptimizationTable(filename):
    
    df=pd.read_csv(filename,sep=" ",index_col=False,header=None,names=["energy","parameter"])
    df["step"]=df.index
    return df
    
def getMinEnergyParameter(df,cutOff=0):
    parameters=df[df["step"]>cutOff]["parameter"]
    energies=df[df["step"]>cutOff]["energy"]
    
    weights=1./energies
    mean=np.average(parameters,weights=weights)
    error=np.sqrt(np.average(parameters**2,weights=weights) - mean**2)

    return [mean,error]

        

def readConfigurations(filename):
    
    with open(filename) as f:
        content = f.readlines()
    configurations=[]
    i=0
    positions=np.array([])
    spins=np.array([])
    
    nWalkers=int(string.split(content[i]," ")[1])
    print nWalkers
    for iW in range(0,nWalkers):
        i+=4
        
        nOrbitals=int(string.split(content[i]," ")[2])
        print nOrbitals
        positions=np.zeros(nOrbitals)
        spins=np.zeros(nOrbitals)
        
        
        orbitals=[]
        for iOrbital in range(0,nOrbitals) :
            i+=1
            orbital=np.array((string.split(content[i])),dtype=float)
            orbitals.append(orbital)
            
        configurations.append(np.array(orbitals))
    return configurations

    
def plot1DSpinorConfig(config):
    config2= config[config[:,2]==1]
    
    plt.plot(config2[:,0],config2[:,0]*0,"ob")
    config2= config[config[:,2]==-1]
    plt.plot(config2[:,0],config2[:,0]*0,"or")
    
    
    
