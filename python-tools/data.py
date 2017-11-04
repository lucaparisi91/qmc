#!/home/luca/software/anaconda/bin/python
import scipy as sp
import numpy as np
import pandas as pd
from math import *
import os
import os.path
import tools

def collect_data(dirname):
    d={}
    if os.path.isfile(dirname+"/input.in"):
        d.update(tools.read_settings(dirname+"/input.in"))
        d["n"]=float(d["particles"].strip().split()[0])
    else:
        return {}
    if os.path.isfile(dirname + "/energy/status.log"):
        d.update( tools.read_settings(dirname+"/energy/status.log"))
    else:
        return {}
    if os.path.isfile(dirname+"/jastrow.in"):
            b=tools.read_settings(dirname+"/jastrow.in")
            d["g"]=b["g"]
            if b["delta"]==0:
                d["g"]="inf"
    else:
        return {}
    if os.path.isfile(dirname+"/jastrowi.in"):
            b=tools.read_settings(dirname+"/jastrowi.in")
            d["g_tilde"]=b["g"]
            if b["delta"]==0:
                d["g_tilde"]="inf"
    else:
        return {}
    return d

def make_db(main_dir="."):
    dataframes=[]
    j=0
    for dirname in os.listdir(main_dir):
        dirname=main_dir+"/"+dirname
        if os.path.isdir(dirname):
            
            
            d=collect_data(dirname)
            if d != {}:
                dataframes.append(pd.DataFrame(d,index=[j]))
                j=j+1
    if dataframes != []:
        return pd.concat(dataframes)
    else:
        return None

