import scipy as sp
import numpy as np
from math import *
import matplotlib.pylab as plt
def errfcApproximation(x):
    #errf return np.sqrt(1-np.exp(-x**2))*(0.5*sqrt(pi) + 31./200.*np.exp(-x**2)-341./8000*np.exp(-2*x**2))*2/sqrt(pi)
    a1=0.254829592
    a2=-.284496736
    a3=1.421413741
    a4=-1.453152027
    a5=1.061405429
    p=0.3275911

    t=1./(1+p*x)
    return x*(a1*t + a2*t**2 + a3*t**3 + a4*t**4 + a5*t**5)
    
def errfcApproximation2(x):
    return 1/sqrt(pi)*(1-1./(2*x**2) + 3./(2*x**2)**2 - 15/(2*x**2)**3)

def errfcApproximation3(x):
    a=2.7889
    
    return a*x/( (a-1)*np.sqrt(pi*x**2) +np.sqrt(pi*x**2 + a**2) )



x1=np.linspace(1,20,num=10000)
x2=np.linspace(0,20,num=10000)
#plt.plot(x,np.sqrt(pi)*sp.special.erfc(x)*np.exp(x**2))
#plt.plot(x,sp.special.erfc(x))
plt.plot(x1,abs( errfcApproximation(x1) - x1*sp.special.erfc(x1)*np.exp(x1**2)))
plt.plot(x2,abs( errfcApproximation3(x2) - x2*sp.special.erfc(x2)*np.exp(x2**2)))

plt.yscale("log")
plt.show()
