import inputFolders
import anal
import jastrow
import pandas as pd
import tools
import numpy as np
import scipy as sp
import models
import scipy.optimize
import scipy.interpolate as interp
import matplotlib.pylab as plt
from math import *
import datetime
import inspect

def gatherData(jumps=0,makePlot=True):
    tab=[]
    for f in inputFolders.folders:
        try:
            n=int(f.get_particles(index=0))*2
            g=abs(float(f.get_g()))
            tau=float(f.get_time_step())
            g_tilde=abs(float(f.get_g_tilde()))
            ratio=np.round(g_tilde/g,2)
            landa=n*4/g**2
            landa_tilde=n*4/g_tilde**2
            walkers=np.int(f.get_mean_walkers())
            
            l=2*abs(float(f.get_trap(0)["center"]))
            beta=2*abs(float(inputFolders.get_value("jastrow_parameter jastrowGauss.in position",f)))
            alpha=abs(float(inputFolders.get_value("jastrow_parameter jastrowGauss.in alpha",f)))
            datae=anal.anal_energy(filename=f.dir_name+"/energy.dat",jumps=jumps,makePlot=makePlot)
            #datae=[0,0]
            datacm=anal.anal_energy(filename=f.dir_name+"/center_of_mass_difference.dat",jumps=jumps,makePlot=False)
            datacm2=anal.anal_energy(filename=f.dir_name+"/center_of_mass_differenceSquared.dat",jumps=jumps,makePlot=False)
            #datacm=[0,0]
            
            tab.append([ratio,landa,n,datae[0]/n,datae[1]/n,datacm[0],datacm[1],l,beta,alpha,landa_tilde,tau,datacm2[0],datacm2[1]])
            
        except IOError:
            print "IO error in " + f.dir_path
        except ValueError:
            print "Empty " + f.dir_path + " "
            
    tab=np.array(tab)

    
    df=pd.DataFrame({"ratio":tab[:,0],"landa":tab[:,1],"n":tab[:,2],"e":tab[:,3],"deltae":tab[:,4],"cm":tab[:,5],"deltacm":tab[:,6],"l":tab[:,7],"beta":tab[:,8],"alpha":tab[:,9],"landa_tilde":tab[:,10],"tau":tab[:,11],"cm2":tab[:,12],"deltacm2":tab[:,13]})
    
    return df

def gatherDataOptBeta(jumps=0,last=False,makePlot=False):
    tab=[]
    for f in inputFolders.folders:
        try:
            n=int(f.get_particles(index=0))*2
            g=abs(float(f.get_g()))
            g_tilde=abs(float(f.get_g_tilde()))
            ratio=np.round(g_tilde/g,2)
            landa=n*4/g**2
            landa_tilde=n*4/g_tilde**2
            p=abs(float(f.getOptParameter(jumps=6,makePlot=makePlot,last=last)))*2
            l=2*abs(float(f.get_trap(0)["center"]))
            #beta=2*abs(float(inputFolders.get_value("jastrow_parameter jastrowGauss.in position",f)))
            alpha=abs(float(inputFolders.get_value("jastrow_parameter jastrowGauss.in alpha",f)))
            #datacm=[0,0]
            tab.append([ratio,landa,n,l,alpha,landa_tilde ,p])
            
        except IOError:
            print "IO error in " + f.dir_path
        except ValueError:
            print "Empty " + f.dir_path + " "
    if tab==[]:
        return pd.DataFrame({})

    
    tab=np.array(tab)

    
    df=pd.DataFrame({"ratio":tab[:,0],"landa":tab[:,1],"n":tab[:,2],"l":tab[:,3],"alpha":tab[:,4],"landa_tilde":tab[:,5],"beta":tab[:,6]})
    
    return df


def gatherDataOptAlpha(jumps=0,last=False,makePlot=False):
    tab=[]
    for f in inputFolders.folders:
        try:
            n=int(f.get_particles(index=0))*2
            g=abs(float(f.get_g()))
            g_tilde=abs(float(f.get_g_tilde()))
            ratio=np.round(g_tilde/g,2)
            landa=n*4/g**2
            landa_tilde=n*4/g_tilde**2
            p=float(f.getOptParameter(jumps=jumps,makePlot=makePlot,last=last))
            l=2*abs(float(f.get_trap(0)["center"]))
            beta=2*abs(float(inputFolders.get_value("jastrow_parameter jastrowGauss.in position",f)))
            #alpha=abs(float(inputFolders.get_value("jastrow_parameter jastrowGauss.in alpha",f)))
            #datacm=[0,0]
            tab.append([ratio,landa,n,l,beta ,landa_tilde ,p])
        except IOError:
            print "IO error in " + f.dir_path
        except ValueError:
            print "Empty " + f.dir_path + " "
    if tab==[]:
        return pd.DataFrame({})

    
    tab=np.array(tab)

    
    df=pd.DataFrame({"ratio":tab[:,0],"landa":tab[:,1],"n":tab[:,2],"l":tab[:,3],"beta":tab[:,4],"landa_tilde":tab[:,5],"alpha":tab[:,6]})
    
    return df

# performs a fit over the positions
def fitPositions(df,f,makePlot=True,lMax=1):
    
    values=df[["landa","landa_tilde"]].drop_duplicates()
    tab=[]
    deg=len(inspect.getargspec(f).args)-1
    x=np.linspace(0,1,num=10000)
    for row in values.iterrows():
        #subset of dataframe for unique landa and landa_tilde
        landa=row[1]["landa"]
        landa_tilde=row[1]["landa_tilde"]
        df1=df[ (abs(df["landa"]-landa)<1e-4) & (abs(df["landa_tilde"]-landa_tilde)<1e-4) ]
        # performs the fit
        df1=df1[ df1["l"]<=lMax]
        x1=np.array(df1["l"])
        y1=np.array(df1["cm"])
        sigma1=np.array(df1["deltacm"])
        fit=sp.optimize.curve_fit(f,x1,y1,sigma=sigma1)
        
        if makePlot:
            plt.errorbar(df1["l"],np.array(df1["cm"]),np.array(df1["deltacm"]),fmt="o"  )
            plt.plot(x,f(x,*fit[0]))
        row=[landa,landa_tilde]
        for i in range(0,deg):
            row.append(fit[0][i])
        for i in range(0,deg):
            row.append(np.sqrt(np.diag(fit[1]))[i])
        tab.append(row)

    #build dataframe of data from table
    tab=np.array(tab)
    dictData={}
    dictData["landa"]=tab[:,0]
    dictData["landa_tilde"]=tab[:,1]
    for i in range(0,deg):
        dictData["p"+str(i)]=tab[:,2+i]
        dictData["deltap"+str(i)]=tab[:,2+deg+i]
    return pd.DataFrame(dictData)
    
def plotSpinDipole(df,f,lMax=1):
    res=fitPositions(df,f,makePlot=False,lMax=lMax)
    x=np.logspace(2,5,num=1000)
    plt.xscale("log")
    plt.errorbar(np.array(res["landa_tilde"]),np.array(1./abs(res["p0"])),np.array(res["deltap0"]/abs(res["p0"])**2),fmt="o")
    plt.plot(x,  (1-np.sqrt(100./x))/(1+np.sqrt(100/x)))
    return res
    
# performs a fit over the positions
def fitBeta(df,f,makePlot=True,lMax=1):
    
    values=df[["landa","landa_tilde"]].drop_duplicates()
    tab=[]
    deg=len(inspect.getargspec(f).args)-1
    x=np.linspace(0,1,num=10000)
    for row in values.iterrows():
        #subset of dataframe for unique landa and landa_tilde
        landa=row[1]["landa"]
        landa_tilde=row[1]["landa_tilde"]
        df1=df[ (abs(df["landa"]-landa)<1e-4) & (abs(df["landa_tilde"]-landa_tilde)<1e-4) ]
        # performs the fit
        df1=df1[ df1["l"]<=lMax]
        x1=np.array(df1["l"])
        y1=np.array(df1["beta"])
        fit=sp.optimize.curve_fit(f,x1,y1)
        
        if makePlot:
            plt.plot(df1["l"],np.array(df1["beta"]),"o")
            plt.plot(x,f(x,*fit[0]))
        row=[landa,landa_tilde]
        for i in range(0,deg):
            row.append(fit[0][i])
        tab.append(row)

    #build dataframe of data from table
    tab=np.array(tab)
    dictData={}
    dictData["landa"]=tab[:,0]
    dictData["landa_tilde"]=tab[:,1]
    for i in range(0,deg):
        dictData["p"+str(i)]=tab[:,2+i]
    return pd.DataFrame(dictData)
    
def filterFutureWalkers(neg=False):
    tmps=[]
    for f in inputFolders.folders:
        isFW=f.get_measure_by_label("center_of_mass_difference").attrib["futureWalkers"]
        if neg:
            if isFW=="false":
                tmps.append(f)
        else:
            if isFW=="true":
                tmps.append(f)
                
    inputFolders.folders=tmps
