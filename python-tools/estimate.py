import numpy as np
import matplotlib.pylab as plt
import scipy as sp
from scipy import optimize
from scipy import stats
from math import *
import tools

############ exceptions

class TooManyBlocks:
    pass

class InvalidInput:
    pass
class notEnoughData:
    pass
# maxes a blocked time series
def block(values,n,minBlockLength=1):

    # checking inputs
    if n==0:
        raise InvalidInput()
    blockLength=len(values)/int(n)

    if blockLength<minBlockLength:
        raise TooManyBlocks()

    # builds means of blocks
    values2=np.zeros(n)
    
    for i in range(0,n-1):
        values2[i]=np.mean(values[i*blockLength:(i+1)*blockLength])

    values2[n-1]=np.mean(values[(n-1)*blockLength:])
    return values2

# estimate error Jack Nife
def resampleJackNife(valuesIn,minBlocks=10,makePlot=True):
    
    if len(valuesIn)<minBlocks:
        raise NotEnoughData()
    
    n=len(valuesIn)
    Ms=[]
    errors=[]
    values=valuesIn
    while(n>=minBlocks):
        # bloc and computes the mean minus one block in values2
        values=block(values,n)
        sumV=np.sum(values)
        m=sumV/n
        values2=np.zeros(n)
        values2[0]=sumV/n
        for i in range(0,n):
            values2[i]=(sumV-values[i])/(n-1)
        # compute the deviation of the new means from original mean
        dev=sqrt(np.sum((values2-m)**2)*n/(n-1))
        errors.append(dev)
        n=n/2
        
    erors=np.array(errors)
    Ms=2**np.linspace(0,len(errors),num=len(errors))

    index2=1./np.sqrt(Ms)
    params,covs=optimize.curve_fit(lambda x,a,b: a*x + b,index2,errors)
    res=stats.linregress(index2,errors)
    errorsFit=[np.sqrt(covs[0][0]),np.sqrt(covs[1][1])]

    if makePlot:
        x=np.linspace(np.min(Ms),np.max(Ms),num=1000)
        plt.plot(Ms,errors,"o")
        plt.plot(x,params[1] + params[0]/np.sqrt(x))
        y1=params[1]-errorsFit[1] + (params[0]-errorsFit[0])/np.sqrt(x)
        y2=params[1]+errorsFit[1] + (params[0]+errorsFit[0])/np.sqrt(x)
        plt.fill_between(x, y1, y2, where=y2 >= y1, interpolate=True,alpha=0.2)
        plt.plot(x,params[1]+0*x ,"--")
    
    return [np.mean(valuesIn),params[1],errorsFit[1],abs(params[0]/(params[1]*sqrt(np.max(Ms))))]
    
# divide the samples in n blocks
def estimateErrorBlocking(valuesIn,minBlocks=10,makePlot=False):
    
    if len(valuesIn)<minBlocks:
        raise NotEnoughData()
        
    values=valuesIn
    errors=[]
    Ms=[]
    i=len(values)
    
    while i>=minBlocks:
        values=block(values,i)
        errors.append(np.sqrt(np.mean(values**2)-np.mean(values)**2)/sqrt(i))
        i=i/2

    errors=np.array(errors)
    Ms=2**np.linspace(0,len(errors),num=len(errors))
   
    index2=1./np.sqrt(Ms)
    params,covs=optimize.curve_fit(lambda x,a,b: a*x + b,index2,errors)
    res=stats.linregress(index2,errors)
    errorsFit=[np.sqrt(covs[0][0]),np.sqrt(covs[1][1])]
    
    if makePlot:
        x=np.linspace(np.min(Ms),np.max(Ms),num=1000)
        plt.plot(Ms,errors,"o")
        plt.plot(x,params[1] + params[0]/np.sqrt(x))
        y1=params[1]-errorsFit[1] + (params[0]-errorsFit[0])/np.sqrt(x)
        y2=params[1]+errorsFit[1] + (params[0]+errorsFit[0])/np.sqrt(x)
        plt.fill_between(x, y1, y2, where=y2 >= y1, interpolate=True,alpha=0.2)
        plt.plot(x,params[1]+0*x ,"--")
    
    return [params[1],errorsFit[1],abs(params[0]/(params[1]*sqrt(np.max(Ms))))]
        
# estimate the error using Jacknife
# divide the samples in n blocks

def distributionValues(y,bins=None):
    if bins==None:
        bins=int(sqrt(len(y)))
    
    hist,bins=np.histogram(y,bins=bins)

    plt.bar(bins[1:],hist,width=bins[1]-bins[0])
    

def getMean(filename,method="blocking",jumps=0,makePlot=False):
    y=np.array(tools.read_matrix_from_file(filename)[1])
    y=y[jumps:]
    err=estimateErrorBlocking(y,minBlocks=10,makePlot=makePlot)
    
    return [np.mean(y),err[0],err[1],err[2]]

#y=np.array(tools.read_matrix_from_file("center_of_mass_differenceSquared.dat")[1])
#x=np.linspace(0,len(y),num=len(y))
#plt.plot(x,y)

