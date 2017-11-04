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
import PBS

reload(PBS)
inputFolders.set_anal_folder()

def plan_impurity_trapped():
    inputFolders.set_anal_folder()
    for method in ["svmc"]:
        for landa in [100]:
            for landa_tilde in [1000,200,600]:
              for n in [80]:
                for alpha in np.linspace(1.85,2,num=10):
                  for x in [0.5]:
                    for l in [0]:
                    #for l in [0]:
                      for nF in [0]:
                        #for beta in np.linspace(l/20.,2*l,num=10.):
                        for beta in [0]:
                         for t in [1e-3]:
                           for w in [1]:
                            if landa==0:
                                g="inf"
                            else:
                                g=sqrt(n*4./landa)
                                
                            g_tilde=sqrt(n*4./landa_tilde)
                            n1=int(n*x)
                            n2=n-n1
                            if (method=="dmc"):
                          	   f1=inputFolders.get_folder(["method","g","n index 0","n index 1","jastrow_parameter jastrowGauss.in position","trap center 0"],["svmc",g,n1,n2,-beta/2.,-l/2.])
                            f=inputFolders.new()
                            l=float(l)
                            f.set_particles(n1,index=0)
                            f.set_particles(n2,index=1)
                            f.set_lBox(30)
                            f.set_method(method)
                            f.set_g_tilde(g_tilde)
                            f.set_g(g)
                            if (nF ==0):
 				FW="false"
			    else:
 				FW="true"	
                            f.get_measure_by_label("center_of_mass_difference").attrib["bins"]=str(nF)
                            f.get_measure_by_label("center_of_mass_difference").attrib["futureWalkers"]=str(FW)
                            f.get_measure_by_label("center_of_mass_differenceSquared").attrib["bins"]=str(nF)
                            f.get_measure_by_label("center_of_mass_differenceSquared").attrib["futureWalkers"]=str(FW)
                            # if no interpolation function is provided
                            
                            f.set_mean_walkers(w)
                            f.set_time_step(t)
                            #f.set_winding_number(bins=500,time=2)
                            #f.make_jastrows()
                            f.set_trap(0,center=-l/2.)
                            f.set_trap(1,center=l/2.)
                            f.save()
                            
                            # inter particle jastrow
                            j=jastrow.jastrow_delta_in_trap(n,g_tilde)
                            
                            j.print_jastrow_parameters(f.dir_path + "/jastrowi.in")
                            # intra particle jastrow
                            j=jastrow.jastrow_delta_in_trap(n,g)

                            j.print_jastrow_parameters(f.dir_path + "/jastrow.in")

                            # gaussian 1
                            j=jastrow.jastrow_gaussian(n,position=-beta/2.,alpha=alpha)
                            j.print_jastrow_parameters(f.dir_path + "/jastrowGauss.in")
                            # gaussian 2
                            j=jastrow.jastrow_gaussian(n,position=beta/2.,alpha=alpha)
                            j.print_jastrow_parameters(f.dir_path + "/jastrowGaussi.in")
                            
                            if (method=="dmc"): 
                            	f.set_initial_condition(f1)

# minimize according to some parameter
def minimizeParameter(x,y,n,sigma=None,makePlot=True):
    x1=min(x)
    x2=max(x)
    
    if sigma is None:
        w=None
    else:
        w=1./np.array(sigma)
    p=np.poly1d(np.polyfit(x,y,n,w=w))
    
    if sigma is None:
        plt.plot(x,y,"or")
    else:
        plt.errorbar(x,y,sigma,fmt="or")
        
    xs=np.linspace(x1,x2,num=100)
    if makePlot==True:
            plt.plot(xs,p(xs))
        
            
    roots=anal.getMinPol(p,x1,x2)
    return min(roots)

def minimizeAlpha(df,makePlot=False,degree=8):
    
    headers=df[["ratio","landa","n","l","beta","landa_tilde"]].drop_duplicates()
    tab=[]
    for header in headers.iterrows():
            ratio=header[1]["ratio"]
            landa=header[1]["landa"]
            n=header[1]["n"]
            l=header[1]["l"]
            beta=header[1]["beta"]
            landa_tilde=header[1]["landa_tilde"]
            df1=df[ (df["l"]==l) & (df["n"]==n) & (df["landa"]==landa) & (df["ratio"]==ratio) & (df["beta"]==beta) & (df["landa_tilde"]==landa_tilde)  ]
            alpha=minimizeParameter(df1["alpha"],df1["e"],degree,sigma=df1["deltae"],makePlot=makePlot)
            tab.append([ratio,landa,n,l,beta,alpha,landa_tilde])
            
    tab=np.array(tab)
    dfOpt=pd.DataFrame({"ratio":tab[:,0],"landa":tab[:,1],"n":tab[:,2],"l":tab[:,3],"beta":tab[:,4],"alpha":tab[:,5],"landa_tilde":tab[:,6]})
    return dfOpt

def startAlpha():
    sch=PBS.PBSScheduler()
    sch.time=datetime.timedelta(minutes=30)
    sch.outFile="qmc.pbs"
    #sch.wait=datetime.datetime.now() + sch.time
    sch.setCommand("mpirun -np " + str(int(sch.ppn*sch.nodes))+" ./qmc");
    inputFolders.schedule_all(sch,maxN=60)
    
