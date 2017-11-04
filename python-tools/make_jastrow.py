#!/home/luca/software/anaconda/bin/python
import scipy as sp
from scipy import optimize
import numpy as np
from math import *
import re

class jastrow:
    def __init__(self):
        # jastrow parameters
        self.c=1
        self.R=5
        self.l_box=10
        
        # function whose root needs to be found
    def f(self,x) :
        return self.c*tan(x) - (pi/2 - x)*1./self.R
    def get_parameters(self):
        
        d={}
        #print str(f(2.7)) + " " +str(f(4))
        delta=optimize.brentq(f,2.0,4.4)
        d["delta"]=self.delta
        d["k"]=self.c*tan(self.delta)
        d["c"]=self.c
        d["a"]=1./sin(d["k"]*self.R+d["delta"])
        d["cut_off"]=self.R
        return

def write_dictionary_as_input(d,filename):
    f=open(filename,"w")
    for key in d:
        f.write(key+"="+ str(d[key])+"\n")
j=jastrow()

write_dictionary_as_input(j.get_parameters(),"jastrow.in")
