#!/home/lucap/anaconda2/bin/python
import scipy as sp
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import math
import tools
#parameters of the integration
def lieb(gammas):
    if tools.check_if_string(gammas) or (not tools.check_if_iterable(gammas)):
        return lieb_value(float(gammas))
    else:
        rs=[]
        if tools.check_if_iterable(gammas):
            for gamma in gammas:
                rs.append(lieb_value(gamma))
            return np.array(rs)
   
def lieb_value(gamma):
    a=-1
    b=1
    n=1000
    N=10
    L=10
    c=1
    eps=1E-10
    dx=(b-a)*1./(n-1)
    gamma=float(gamma)
    #internal variables
    i=0
    err=1000
    #vector of the components of the system
    g_old=np.empty(n)
    g_new=np.empty(n)
    g_old.fill(1)
    x_vec=np.empty(n)
    for i in range(0,n):
        x_vec[i]=a+i*dx

        #self-consistent iteration to find the correct solution
    e_vec=[]
    g_old.fill(1)
    err=100
    while eps<err:
        landa=sp.integrate.simps(g_old,dx=dx)*gamma
        for i in range(0,n):
            buff=g_old/(landa**2+(x_vec - x_vec[i])**2)
            g_new[i]=(sp.integrate.simps(buff,dx=dx)*2*landa+1)/(2*math.pi)
        err=np.linalg.norm(g_new-g_old)
        
        g_old=np.copy(g_new)
    
    e=sp.integrate.simps(g_old*(x_vec**2),dx=dx)*1.0*(gamma/landa)**3
    return e/2

def lieb_derivative(gamma,h):
    return (lieb_value(gamma + h) - lieb_value(gamma - h))/(2*h)
def chemical_potential(gamma,h):
    return 3*lieb_value(gamma)-gamma*lieb_derivative(gamma,h)

    
