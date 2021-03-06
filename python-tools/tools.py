#!/home/luca/software/anaconda/bin/python
import re
import os
import numpy as np
import matplotlib.pylab as plt
import  xml.etree.ElementTree as ET
from scipy.optimize import curve_fit
import anal
from math import *

class absent_file:
        def warn(self):
                print "Data not found."
class bad_table:
        def warn(self):
                print "Wrong format of the table."


def write_dictionary_in_file(d,filename):
    f=open(filename,"w+")
    for key in d:
        f.write(key+"="+ str(d[key])+"\n")
    f.close()
def dictionary_to_xml(d,root):
    for key in d:
        ET.SubElement(root, str(key)).text =str(d[key])
def write_array_in_file(a,filename):
    if (len(a) == 0 ):
        return
    f=open(filename,"w+")
    for el in a:
        f.write(str(el)+"\n")
    f.close()
def write_x_y(d,filename):
    f=open(filename,"w+")
    for i in range(0,len(d[0])):
        f.write(str(d[0][i])+ " " + str(d[1][i])+"\n")
def write_matrix_in_file(m,filename):
    m=np.array(m)
    if (m.shape == (0,0)):
        return
    
    f=open(filename,"w+")
    
    for i in range (0,m.shape[0]):
        for j in range(0,m.shape[1]):
            f.write(str(m[i][j]) + " " )
        f.write("\n")
    f.close()
        
def get_name(filename):
    if ( re.match(".*/.+",filename)):
             return re.match(".*/(.*)$",filename).group(1)
    else:
        return filename
def read_settings(filename):
    d={}
    line_name=""
    line_value=""
    inline=0
    f=open(filename)
    for line in f :
        line=line.strip()
        if line=="":
            continue
        if (line[0] == "#") :            
            continue
        if inline == 1:
            if bool(re.search("=",line)):
                    
                    d[line_name]=line_value
                    line_name=""
                    line_value=""
                    inline=0
            else:
                    
                    line_value=line_value+line
                    continue
        
        if bool(re.search("=>",line)):
            #it is a vector
            line=line.split('=>')
            inline=1
            line_name=line[0]
            
        else:
            if bool(re.search("=",line)):
                # is a scalar
                line=line.split('=')
                d[line[0]]=line[1]
            
            
    return d
        
    f.close()

def check_pid(pid):
    pid=int(pid)
    """ Check For the existence of a unix pid. """
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    else:
        return True

def read_matrix_from_file(filename,nMax=None):
    rows=[]
    first_line=True
    line_len=0
    n=0
    # read e.dat
    f=open(filename)
    for line in f :
        if (nMax!=None ):
                if(n>nMax):
                        break
        try:
            
            line=re.sub(" +",' ',line)
            line=re.sub("^ *",'',line)
            line=re.sub(" *$",'',line)
            line_list=line.split()

        except:

            print "could not format the line.Skipping"

        if first_line:
            # allocate a vector of observables
            first_line=False
            line_len=len(line_list)
            for i in range(0,line_len):   
                rows.append([])
        else:
            try:
                if len(line_list) != line_len:
                    raise bad_line_format
                
            except:
                print "Wrong number of elements in the line. Skipping"
                continue
        
        for i  in range(0,line_len):
            try:
                e=float(line_list[i])
            except:
                print "Cannot convert" + line_list[i] + "element to float. Skipping value."
                continue
            
            rows[i].append(e)
            n=n+1
            
    f.close()
    return np.array(rows)
    
def read_array_from_file(filename,begin=0,end=0,e_min=None,e_max=None):
    
    values=[]
    i=0

    if e_min=="None":
        e_min=None
    if e_max=="None":
        e_max=None
    # read e.dat
    f=open(filename)
    for line in f :
        
        i=i+1
        #print i
        if i<begin:
            continue
        else:
            
            if (i>end and end > 0):
                break
        try:
            e=float(line.split()[0])
        except:
            print "cannot convert '"+line.split()[0]+"' to float"
            continue
        if e_min !=None :
            if e < e_min:
                continue
        if e_max !=None :
            if e > e_max:
                continue
            
        values.append(e)
        
    f.close()
     
    return np.array(values)
def plot_function(function_to_plot,x1=0,x2=10,bins=1000):
    vec=[]
    xs=[]
    for x in np.arange(x1,x2,(x2-x1)*1./bins):
        vec.append(function_to_plot(x))
        xs.append(x)
    plt.plot(xs,vec)

def dictionary_to_string(d):
    s=""
    for key in d:
        s=s+str(key)+"="+str(d[key])+";"
    return s
def get_scalar_value(filename):
    if (os.path.exists(filename)):
        
        with open(filename,"r") as f:
            params=f.readline().split();
    else:
        raise absent_file
    
    if len(params) >=3:
            return [float(params[0]),float(params[1]),float(params[2])]
    else:
            return [float(params[0]),float(params[1])]

def plot_table(table):
        table1=np.array(table)
        table1=table[table1[:,0].argsort()]
        if len(table1.shape) != 2:
            raise bad_table
        if table1.shape[1]>=3:
            plt.errorbar(table1[:,0],table1[:,1],yerr=table1[:,2],fmt='o')
        else:
                if table.shape[1]>=2:
                        plt.plot(table1[:,0],table1[:,1],marker='o')
                else:
                        raise bad_table

def print_table(table,filename):
    table=np.array(table)
    
    with open(filename,"w") as f:
        
        f.write(str(table.shape[0]) + " " + str(table.shape[1]) + "\n")
        for i in range(0,table.shape[0]):
            for j in range(0,table.shape[1]):
                f.write(str(table[i,j]) + " ")
            f.write("\n")
def linear(x,a,b):
    return a*x+b

def fit_table(table,func=linear,cut_low=0,cut_heigh=None):
    table=np.array(table)
    
    if cut_heigh==None:
        cut_heigh=table.shape[0]-1
    if len(table.shape)!=2:
            raise bad_table
    if table.shape[1]==3:
            
            params,covs=curve_fit(func,table[cut_low:cut_heigh,0],table[cut_low:cut_heigh,1],sigma=table[cut_low:cut_heigh,2],maxfev=100000)
    else:
            if table.shape[1]==2:
                params,covs=curve_fit(func,table[cut_low:cut_heigh,0],table[cut_low:cut_heigh,1],maxfev=100000)
            else:
                    raise bad_table

    return [params,covs]
# performs the fit over a certain table

def plot_fit(x,fit_res,bins=100):
    x=np.append(x,0)
    x=np.sort(x)
    plt.plot(x,fit_res[0][0]*np.array(x)+fit_res[0][1])
    plt.fill_between(x, (fit_res[0][0]-sqrt(fit_res[1][0][0]))*np.array(x)+fit_res[0][1]-sqrt(fit_res[1][1][1]),(fit_res[0][0]+sqrt(fit_res[1][0][0]))*np.array(x)+fit_res[0][1]+sqrt(fit_res[1][1][1]),interpolate=True,facecolor="green",edgecolor=None,alpha=0.4)
    
def build_label(a1,a2):
    
    
    label=""
    if check_if_string(a1):
        label=str(a1) + "=" + str(a2)
    else:
        if check_if_iterable(a1):
            
            for i in range(0,len(a1)):
                if i!= 0:
                    label=label+","
                label=label + str(a1[i]) + "=" + str(a2[i])
        else:
                    return None
    
    return label 
    
                
def to_float(vec):
    a=[]
    if check_if_string(vec):
        a=float(vec)
    else:
        if check_if_iterable(vec):
            for el in vec:
                a.append(float(el))
            return a
        else:
            a=None
        
    return a
def check_if_string(a):
    if isinstance(a,basestring):
        return True
    else:
        return False
def check_if_iterable(a):
    if  hasattr(a, "__len__"):
        return True
    else:
        return False

class string_not_a_bool:
    pass

def string_to_bool(s):
    if s=="True":
        return True
    else:
        if s=="False":
            return False
        else:
            raise string_not_a_bool
