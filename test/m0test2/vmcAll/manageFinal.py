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
reload(inputFolders)

inputFolders.set_anal_folder()
def plan_impurity_trapped(df):
    inputFolders.set_anal_folder()
    for method in ["dmc"]:
        for row in df.iterrows():    
            n=row[1]["n"]
            landa=row[1]["landa"]
            landa_tilde=row[1]["landa_tilde"]
            alpha=row[1]["alpha"]
            beta=row[1]["beta"]
            l=row[1]["l"]
            landa=row[1]["landa"]
            for nF in [0,100,500,1000]:
              for t in [5e-4]:
                        for w in [1000]:
                            x=0.5
                            g=sqrt(n*4./landa)
                            g_tilde=sqrt(n*4./landa_tilde)
                            n1=int(n*x)
                            n2=n-n1
                            if (method=="dmc"):
                                try:
                          	    f1=inputFolders.get_folder(["method","n index 0","n index 1"],["svmc",n1,n2])
                                except IOError:
                                    print "Could not find VMC calculation."
                                    continue

                            f=inputFolders.new()
                            print f.dir_path
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
                            
                            # if using a DMC monte carlo method
                            if (method=="dmc"):
                                try:
                            	    f.set_initial_condition(f1)
                                except IOError:
                                    print "Cannot set the initial condition for " + f.dir_path
                                    f.removeFolder()
                                    
                                    
def minimizeVMCDMCDifference(dfVMC,dfDMC,makePlot=False,deg=8):
    df=pd.merge(dfVMC,dfDMC,on=["ratio","landa","n","l","beta","landa_tilde","alpha"] )
    
    headers=df[["ratio","landa","n","l","landa_tilde","alpha"]].drop_duplicates()
    tab=[]
    for header in headers.iterrows():
        try:
            ratio=header[1]["ratio"]
            landa=header[1]["landa"]
            landa_tilde=header[1]["landa_tilde"]
            n=header[1]["n"]
            l=header[1]["l"]
            alpha=header[1]["alpha"]
            df1=df[ (df["l"]==l) & (df["n"]==n) & (df["landa"]==landa) & (df["ratio"]==ratio) & (df["landa_tilde"]==landa_tilde)]
            
            #m=interp.interp1d(df1["beta"],df1["cm_x"]-df1["cm_y"])
            m1=np.poly1d(np.polyfit(df1["beta"],df1["cm_x"],deg))
            m2=np.poly1d(np.polyfit(df1["beta"],df1["cm_y"],deg))
            beta1=np.min(df1["beta"])
            beta2=np.max(df1["beta"])
            
            if makePlot==True:
                x=np.linspace(beta1,beta2)
                plt.plot(x,m1(x),"r")
                plt.plot(x,m2(x),"b")
                plt.errorbar(np.array(df1["beta"]),np.array(df1["cm_x"]),np.array(df1["deltacm_x"]),fmt="or")
                plt.errorbar(np.array(df1["beta"]),np.array(df1["cm_y"]),np.array(df1["deltacm_y"]),fmt="ob")
            betaOpt=sp.optimize.brentq(lambda x: m1(x)-m2(x),beta1,beta2)
            tab.append([ratio,landa,n,l,betaOpt,landa_tilde,alpha])
        except ValueError:
            print "No minimum for " + str(header)   
    if tab==[]:
	return None
    
    tab=np.array(tab)
    dfOpt=pd.DataFrame({"ratio":tab[:,0],"landa":tab[:,1],"n":tab[:,2],"l":tab[:,3],"beta":tab[:,4],"landa_tilde":tab[:,5],"alpha":tab[:,6]})
    return dfOpt

def startBeta():
    sch=PBS.PBSScheduler()
    sch.time=datetime.timedelta(hours=8)
    sch.outFile="qmc.pbs"
    sch.ppn=5
    sch.nodes=2
    #sch.wait=datetime.datetime.now() + sch.time
    sch.setCommand("mpirun -np " + str(int(sch.ppn*sch.nodes))+" ./qmc");
    inputFolders.schedule_all(sch,maxN=10)
    
