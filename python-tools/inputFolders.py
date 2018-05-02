import re
import os
import shutil
import hashlib
import tools
import jastrow
import os.path
import random
import filecmp
import anal
import xml
import numpy as np
import  xml.etree.ElementTree as ET
import subprocess
import exceptions
import matplotlib.pylab as plt
import time
import pandas as pd
import datetime
from math import *

class nothing_to_do:
    pass
class bad_jastrow_parameter:
    pass
class missing_jastrow_parameter:
    pass
class no_process:
    def warn(self):
            print "No process\n"
            
class process_running:
    def warn(self):
            print "Process already running\n"
            
class absent_file:
        def __init__(self,filename=""):
                print "file: " + filename + "\n"
        def warn(self):
                print "Data not found.\n"

class bad_directory:
        pass

code_dir="/home/luca/source/cqmc"
prevRun_dir=["/home/luca/source/cqmc/test"]
run_dir="/home/luca/source/cqmc/test"
default_files=["input.xml","qmc"]
required_files=["input.xml"]
folders=[]


def get_jastrow_paramer(f,filename,parameter):
    settings_file=f.dir_path + "/"+filename
    root=ET.parse(settings_file).getroot()
    # finds the optimal value for some kind of parameter
    return root.find(parameter).text

# returns the variables i was looking for during the simulation
def get_values(qs,f):
    rs=[]
    if tools.check_if_string(qs):
        return get_value(qs,f)
    else:
        if tools.check_if_iterable(qs):
            for q in qs:
                rs.append(get_value(q,f))
            return rs
        else:
            return None
        
def get_label(qs,f):
    rs=get_values(qs,f)
    #print rs
    return tools.build_label(qs,rs)

def get_value(q,f):
    
    if q=="g_tilde":
        return get_g_tilde(f)
    else:
        if q=="g":
            return f.get_g()
        if q=="ratio":
            return float(f.get_g_tilde())/float(f.get_g())
        else:
            if q=="n":
                return get_n_particles(f)
            else:
                if q=="time_step":
                    return f.get_time_step()
                else:
                    if q=="walkers":
                        return f.get_mean_walkers()
                    else:
                        if q=="method":
                            return f.get_method()
                        
                        else:     
                                p=re.compile("jastrow_parameter.*")
                                if p.match(q)!=None:
                                
                                    filename=q.split(" ")[1]
                                    param_name=q.split(" ")[2]
                                    return get_jastrow_paramer(f,filename,param_name)
                                
                                # return the nunumber of particles in the ith set
                                p=re.compile("n index .*")
                                if p.match(q)!=None:
                                    return f.get_particles(index=q.split(" ")[2])
                                # returns the center of the ith trap
                                p=re.compile("trap center .*")
                                if p.match(q)!=None:
                                    return f.get_trap(int(q.split(" ")[2]))["center"]
                                        
                                
    raise not_found

def get_n_particles(f):
    return f.root.find("system").find("particles").attrib["n"]

def get_g_tilde(f):
    return f.root.find("system").find("g_tilde").text


def match_target_single(f,target,src_target,eps=1e-4):
    try :
        src_target_tmp=float(src_target)
        dst_target=float(get_values(target,f))
        
        # handles the infinity case
        if src_target_tmp==float("inf"):
            
            if dst_target==float("inf"):
                return True
            else:
                return False
            
        if abs(src_target_tmp-dst_target) < eps:
            return True
        else:
            return False
    except:
        src_target_tmp=str(src_target)
        dst_target=str(get_values(target,f))
    #print str(src_target_tmp) + " " + str(dst_target) + "\n"
        if src_target_tmp==dst_target:
            return True
        else:
            return False
    
def match_target(f,target,src_target,eps=1e-4):
    res=True
    
    if tools.check_if_string(src_target) or (not tools.check_if_iterable(src_target) ):
        
        return match_target_single(f,target,src_target,eps=eps)
    else:
        
        for i in range(0,len(src_target)):
            
            if not match_target_single(f,target[i],src_target[i],eps=eps):
                res=False
    return res
        
    
# set the folder to be analized

def set_anal_folder(dirname="."):
    global run_dir
    
    
    dirname=os.path.abspath(dirname)
    #must not end with a slash
    
    if dirname[-1]=="/":
        dirname=dirname[:-1]
    run_dir=dirname
    
    scan()
    
def get_anal_folder():
    return run_dir

def check_valid(folder_path):
    
        for file_name in required_files:
            if (os.path.exists(folder_path+"/"+file_name)):
                return True
            else:
                return False
            
def scan():
    
    global folders
    folders=[]
    #print run_dir
    for dirname in os.listdir(run_dir):
        
        
        if ( dirname=="time_step" or dirname=="walkers"):
                continue
        dirname=run_dir+"/"+dirname        
        
        if os.path.isdir(dirname) and check_valid(dirname):
            folders.append(folder())
            folders[len(folders)-1].load_folder(dirname)
def new(dirname=None):
    global folders
    folders.append(folder())
    folders[len(folders)-1].create(dirname=dirname)
    return folders[len(folders)-1]

class folder:
    
    def __init__(self):
        self.tree=None
        self.root=None
        self.vmc_optimize_folder=""
        self.dir_path=None
        self.dir_name=None
        self.opened=False
        self.running=False
        self.pid=0
        self.qued=False

        # optimize a variational parameter from a set of possible values
    def update_vmc_optimize(self,time_step=1e-3):
        prevFolder=get_anal_folder()
        set_anal_folder(self.vmc_optimize_folder)
        
        # for each value create a folder for the estimation and start the evaluation
        
        try:
            for f in folders:
                
                f.set_time_step(time_step)
                f.save()
        except:
            print "Something went wrong"
        
        set_anal_folder(prevFolder)
        self.save()
        return self
    def rm_vmc_optimize(self):

        
        prevFolder=get_anal_folder()
        set_anal_folder(self.vmc_optimize_folder)
        
        # for each value create a folder for the estimation and start the evaluation
        
        try:
            for f in folders:
                f.removeFolder()
        
        except no_process:
            pass

        self.vmc_optimize_folder=""
        set_anal_folder(prevFolder)
        self.save()
        return self
    def create_vmc_optimize(self,values):
        # create a sub directory where to perform the variational procedure
        prevFolder=get_anal_folder()
        
        
        self.vmc_optimize_folder=self.dir_path + "/var_par_est"
            
        create_anal_folder(self.vmc_optimize_folder)
        try:
            shutil.copy2(self.dir_path+"/input.xml",self.vmc_optimize_folder)
        except:
            raise absent_file("input.xml")
        
        set_anal_folder(self.vmc_optimize_folder)
        
        # for each value create a folder for the estimation and start the evaluation
        
        try:
            for value in values:
                f=new()
                f.set_method("svmc")
                f.set_mean_walkers("1")
                f.set_time_step("0.001")
                f.make_jastrows(rs=value)
                f.save()
                # return to the previous folder i was working with
        except:
            print "Something went wrong"
        
        set_anal_folder(prevFolder)
        self.save()
        return self
    
    def start_vmc_optimize(self):
        
        prevFolder=get_anal_folder()
        try:
            if self.vmc_optimize_folder=="" or None:
                raise nothing_to_do
            
            set_anal_folder(self.vmc_optimize_folder)
            
            #print len(folders)
            if len(folders)==0:
                raise nothing_to_do
            for f in folders:
                f.start()
        except nothing_to_do:
            print "Nothing done"
        except no_process as e:
            e.warn()
        
        set_anal_folder(prevFolder)
        
    # stop to optimize in a vmc optimization
    def stop_vmc_optimize(self):
         
         prevFolder=get_anal_folder()
         try:
             if self.vmc_optimize_folder==" ":
                 raise nothing_to_do
         
             set_anal_folder(self.vmc_optimize_folder)
             if len(folders)==0:
                 raise nothing_to_do
             
             for f in folders:
                 f.stop()
         except nothing_to_do:
             pass
         except no_process as e:
             e.warn()
         except:
             print "Something went wrong."
  
         set_anal_folder(prevFolder)
    def get_vmc_energy_table(jumps=0):
        tab=[]
        prevFolder=get_anal_folder()
        try:
             if self.vmc_optimize_folder==" ":
                 raise nothing_to_do
         
             set_anal_folder(self.vmc_optimize_folder)
             if len(folders)==0:
                 raise nothing_to_do
             
             for f in folders:
                 f.anal_scalar(filename="energy.dat",jumps=jumps)

             tab=anal.get_energy_table(target="jastrow_parameter jastrowi.in cut_off")
             
        except nothing_to_do:
             pass
         
        except:
             print "Something went wrong."
  
        set_anal_folder(prevFolder)
        return np.array(tab)
        
    # saves everything to an external file
    def save(self):
        # append to the ouput xml file the vmc optimization folder
        
        self.o_tree.find("vmc_optimize_folder").text=self.vmc_optimize_folder
        
        if (self.tree != None):
            self.tree.write(self.dir_path+"/input.xml");
            
        self.o_tree.write(self.dir_path+ "/out.xml")    
    def create(self,dirname=None):
        #creates a unique dirname(to use as a label)
        if dirname==None:
            dirname=str(random.randint(0,1000000))
            md5=hashlib.md5()
            md5.update(dirname)
            dirname=md5.hexdigest()[-8:]
            dirname=re.sub(" ","_",dirname)
            
        self.dir_name=dirname
        dirpath=run_dir+"/"+dirname
        self.dir_path=dirpath
        #creates the fold if necessary
        if (not os.path.exists(dirname)) and  (dirname!="."):
            os.makedirs(dirname)
        # copy default files in global variables
        
        
        if os.path.abspath(run_dir) != os.path.abspath(dirname):
            for file_2c in default_files:
                shutil.copy2(run_dir + "/" + file_2c,dirname)
        self.opened=True
        self.load_input()
        #self.load_ouput()
        
    def load_input(self,filename=""):
        if (filename == ""):
            filename=self.dir_path + "/input.xml"
        self.tree=ET.parse(filename)
        
        self.root=self.tree.getroot()
        # load process status
        process_status_file=self.dir_path + "/process_status.log"
        if (os.path.exists(process_status_file)):
            d=tools.read_settings(process_status_file)
            self.pid=d["pid"]
            self.running=tools.string_to_bool(d["running"])
            if "qued" in d:
                self.qued=tools.string_to_bool(d["qued"])
        self.load_ouput()
        
    def load_ouput(self,filename=""):
        if (filename == ""):
            filename=self.dir_path + "/out.xml"
            # open/creates xml tree
        
        if (os.path.exists(filename)):
            self.o_tree=ET.parse(filename)
            self.o_root=self.o_tree.getroot()
            if (self.o_root.find("vmc_optimize_folder")!= None):
                
                self.vmc_optimize_folder=self.o_root.find("vmc_optimize_folder").text
            
        else:
            self.o_root = ET.Element("out")
           
            self.o_tree= ET.ElementTree(self.o_root)
            
        
        cur=self.o_root.find("mean_energy")
        if (cur == None):
            ET.SubElement(self.o_root,"mean_energy")
            
        cur=self.o_root.find("error_energy")
        if (cur == None):
            ET.SubElement(self.o_root,"error_energy")
            
        cur=self.o_root.find("conv_energy")
        if cur==None:
            ET.SubElement(self.o_root,"conv_energy")

        cur=self.o_root.find("sigma_energy")
        if cur==None:
            ET.SubElement(self.o_root,"sigma_energy")
        
        cur=self.o_root.find("vmc_optimize_folder")
        if (cur == None):
            ET.SubElement(self.o_root,"vmc_optimize_folder")
        
    def start(self):
        
            if not self.running:
                    #p=subprocess.Popen(run_dir + "/launch.sh " + self.dir_path,shell=True,stdout=subprocess.PIPE)
                    
                    
                    self.pid=subprocess.check_output([run_dir + "/launch.sh",self.dir_path])
                    self.running=True
                    # saves info to a file
                    self.save_process_status()
            else:
                    # exception: the file is still running
                    raise process_running
            
    # place on a que
    def que(self):
        self.qued=True
    
    def check_pid(self):
        try:
            os.kill(self.pid, 0)
        
        except OSError:
            return False
        else:
            return True
    def update_status(self):
        self.running=self.check_pid()
        
    def save_process_status(self):
            d={"pid":self.pid,"running":self.running,"qued":self.qued}
            tools.write_dictionary_in_file(d,self.dir_path + "/process_status.log")
            
    def stop(self):
            is_zero=False
            try:
                if int(self.pid)==0:
                    is_zero=True
            except:
                pass
            
            if self.running and (not is_zero) :
                    
                    p=subprocess.Popen(run_dir + "/stop.sh " + self.pid,shell=True)
                    
                    self.pid=0
                    self.running=False
                    self.save_process_status()
            else:
                    self.pid=0
                    self.running=False
                    self.save_process_status()
                    
                    raise no_process
                    
    # copy the initial configuration from a different run
    def set_initial_condition(self,f):
        list_init=["save.xml","walkers.in"]
        for file_init in list_init:
            shutil.copy2(f.dir_path + "/" + file_init,self.dir_path)
    
    
    
    def is_stationary(self,filename="energy.dat",bins=10,eps=0.01):
        return anal.is_stationary(self.dir_path + "/" + filename,bins=bins,eps=eps)
            
    def set_mean_walkers(self,w):
            self.root.find("method").find("mean_walkers").text=str(w)
    def get_mean_walkers(self):
            return self.root.find("method").find("mean_walkers").text    
    def set_time_step(self,t):
            self.root.find("method").find("delta_tau").text=str(t)
    def get_time_step(self):
            return self.root.find("method").find("delta_tau").text
    def set_lBox(self,l):
        
        self.root.find("system").find("lBox").text=str(l)
    def get_lBox(self):
        if( self.root.find("system").find("lBox") is not None):
            return self.root.find("system").find("lBox").text
        else:
            return self.get_particles()
    def set_particles(self,n,index=None):
        if (index==None):
            self.root.find("system").find("particles").attrib["n"]=str(n)
        else:
            index=int(index)
            for el in self.root.find("system").iter("particles"):
                if int(el.attrib["index"])==index:
                    el.attrib["n"]=str(n)
    
    def get_particles(self,index=None):
        if index==None:
            return self.root.find("system").find("particles").attrib["n"]
        else:
            index=int(index)
            for el in self.root.find("system").iter("particles"):
                if int(el.attrib["index"])==index:
                    return el.attrib["n"]
    def get_g_tilde(self):
            return self.root.find("system").find("g_tilde").text
    def get_g(self):
            return self.root.find("system").find("g").text
    def set_g_tilde(self,g):
            self.root.find("system").find("g_tilde").text=str(g)
    def set_g(self,g):
            self.root.find("system").find("g").text=str(g)
    def get_method(self):
            return self.root.find("method").attrib["kind"]
    def set_trap(self,i,center=0,omega=1):
        traps=self.root.find("system").find("oneBodyPotential").findall("oneBodyPotential")
        for trap in traps:
            if int(trap.attrib["setA"])==i:
                trap.attrib["omega"]=str(omega)
                trap.attrib["center"]=str(center)
                
    def get_trap(self,i):
        traps=self.root.find("system").find("oneBodyPotential").findall("oneBodyPotential")
        trap=traps[i]
        return {"omega":float(trap.attrib["omega"]),"center":float(trap.attrib["center"])}
        
    def set_method(self,meth):
            self.root.find("method").attrib["kind"]=meth
    # return (bins,skip,max_time) of the winding number
    
    def get_measure_by_label(self,label):
        for measure in self.root.find("measures"):
            if "label" in measure.attrib:
                label2=measure.attrib["label"]
            else:
                label2=measure.tag
            if label==label2:
                return measure
        return None
    
    def get_winding_number(self):
        skip=self.root.find("measures").find("winding_number").attrib["skip"]
        bins=self.root.find("measures").find("winding_number").attrib["bins"]
        return (bins,skip, float( self.get_time_step() ) *(float(skip) + 1)*int(bins) )
    
    def set_winding_number(self,label,bins,time):
        
            steps=int( float(time)/float( self.get_time_step() ) ) 
            skip=steps/int(bins) - 1
            if (skip<0):
                    skip=0
                    bins=steps
            
            self.get_measure_by_label(label).attrib["bins"]=str(bins)
            self.get_measure_by_label(label).attrib["skip"]=str(skip)
            
    def load_folder(self,dir_path):
        self.dir_path=dir_path
        self.dir_name=os.path.basename(dir_path)
        self.opened=True
        if check_valid(dir_path):
            self.load_input()
            #self.load_ouput()
            return True
        else:
            return False
    def get_array(self,filename):
        values=tools.read_matrix_from_file(self.dir_path+"/"+filename)
        x=[k for k in range(0,len(values[0]))]
        tab=np.zeros((len(values[0]),3))
        tab[:,0]=x
        tab[:,1]=values[1]
        tab[:,3]=values[2]
        return tab

    def removeFolder(self):
        
        if not os.path.exists(run_dir+"/backup"):
            os.makedirs(run_dir+"/backup")
            
        shutil.move(self.dir_path,run_dir+"/backup/"+self.dir_name)
    def plot_file(self,filename,error=True,label="",reset_index=False,jumps=0):
            anal.plot_file(self.dir_path+"/"+filename,error,label,reset_index,jumps=jumps)
            
            return self
        
    # returns a dataframe with the evolution in imaginary time
    def winding_imaginary_time(self,label,jumps=0):
        attrib=self.get_measure_by_label(label).attrib
        winding=self.get_measure_by_label(label)
        setA=int(attrib["setA"])
        skip=int(attrib["skip"])
        if "op" in attrib:
            op=attrib["op"]
            setB=int(attrib["setA"])
            nA=int(get_value("n index "+str(setA),self))
            nB=int(get_value("n index "+str(setB),self))
            n=nA+nB
        else:
            n=int(get_value("n index "+str(setA),self))
        
        time_step=float(self.get_time_step());
        
        res=np.array(anal.analVectorHistory(self.dir_path+"/"+label+".dat",makePlot=False,jumps=jumps))
        res[:,0]=res[:,0]*time_step*(skip+1)
        res[:,1]=n*res[:,1]
        res[:,2]=n*res[:,2]
        df=pd.DataFrame({"t":res[:,0],"W":res[:,1],"deltaW":res[:,2]})
        
        return df
        
    def plot_effective_mass(self,filename="winding_number_t.dat",label=""):
            
            anal.plot_effective_mass(self.dir_path+"/"+filename,settings_file=self.dir_path + "/input.xml",param_file=self.dir_path+"/effective_mass_t.dat",label=label)
            
    def anal_effective_mass(self,data="winding_number_t.dat",settings="input.xml",out="effective_mass_t.dat",cut_low=0,cut_heigh=None,n_cuts=10,factor=1):            
            anal.anal_effective_mass(self.dir_path+"/" + data,settings_file=(self.dir_path+"/"+settings),out_file=(self.dir_path+"/"+out),cut_low=cut_low,cut_heigh=cut_heigh,n_cuts=n_cuts,factor=factor)

    def anal_scalar(self,filename="energy.dat",jumps=0,bins=None,eps=0.01):
        
        # analiza a scalar value
        dirName=self.dir_path+"/"+filename[:-4]
        res=anal.anal_energy(self.dir_path,jumps=jumps,bins=bins,filename=filename,sub_directory=dirName)
        
        
        error="?"
        error=res[1]
        mean=res[0]
        conv=res[2]
        sigma=res[3]
        
        self.o_root.find("mean_energy").text=str(mean)
        self.o_root.find("error_energy").text=str(error)
        self.o_root.find("sigma_energy").text=str(sigma)
        self.o_root.find("conv_energy").text=str(conv)
        
        
        with open(self.dir_path+"/"+filename[:-4]+"_res_anal.dat","w") as f:
            #print self.dir_path+"/"+filename[:-4]+"_t.dat"
            f.write(str(mean) + " " + str(error) + " " + str(conv)  + " " + str(sigma) )
            #print str(r.mean) + " +/ " + error
        
            #saves the reblocking info
            #m=np.zeros((len(r.indices),2))
            #m[:,0]=r.indices
            #m[:,1]=r.errors
            #tools.write_matrix_in_file(m,self.dir_path+ "/"+ filename + ".err")
            self.save()

    # make the jastrow of a system with an impurity
    def make_jastrows(self,g=None,g_tilde=None,rs=0.2):
        #print d
        if (self.root.find("system").find("lBox") != None):
            n=int(self.root.find("system").find("lBox").text)      
        else:
            n=int(self.root.find("system").find("particles").attrib["n"])
        if g==None:
            
            if (self.root.find("system").find("g") != None):
                g=self.root.find("system").find("g").text
        if g_tilde==None:
            if (self.root.find("system").find("g_tilde") != None):
                g_tilde=self.root.find("system").find("g_tilde").text

        if g==None or g_tilde==None:
            raise missing_jastrow_parameter
        
        for jastrowElement in self.root.iter("jastrow"):
            filename=jastrowElement.get("file")
            kind=jastrowElement.get("kind")
            if filename=="jastrow.in":
                inter=g
            else:
                inter=g_tilde
            
            if kind=="delta":
                
                if inter=="inf":
                    j=jastrow.TG_jastrow(int(n))
                else:
                    j=jastrow.jastrow_general(n,float(inter))
            else:
                
                if kind=="delta_bound_state":
                    if inter=="inf":
                        raise bad_jastrow_parameter
                    
                    j=jastrow.delta_bound_state(float(inter),int(n),float(rs))
        
            j.print_jastrow_parameters(self.dir_path+"/"+filename)
                
        
def set(filename,prop,value):
    f=open("tmp_input.txt","w")
    g=open(filename)
    prop=str(prop)
    value=str(value)
    inline=0
    for line in g:
        
        line=line.strip()
        if inline==1 :
            if re.match(".*=.*",line.strip()) :
                inline=0
            else:
                continue
        
        if re.match("^"+prop+"=>.*$",line.strip()) :
            f.write(prop+"=>"+"\n")
            f.write(value+"\n")
            inline=1
                
        else:
                 
            line=re.sub("^"+prop+"=.*",prop+"="+value,line.strip())
            
            f.write(line+"\n")
    shutil.move("tmp_input.txt",filename)
    

# creates the need folder by a dictionary d
# must contain jastrow parameters (g,g_tilde)


def check_dir_valid(dirname):
    
    if (os.path.exists(dirname+"/input.xml") and os.path.exists(dirname+"/jastrow.in") and os.path.exists(dirname+"/jastrowi.in") ):
            return True
    return False
            
'''
def show_input(main_dir=run_dir_default,model={}):
    
    for dirname in os.listdir(main_dir):
        shortname=dirname
        dirname=main_dir+"/"+dirname
        if os.path.isdir(dirname) and check_dir_valid(dirname):
            
            d=tools.read_settings(dirname+"/input.in")
            d["g"]=tools.read_settings(dirname+"/jastrow.in")["g"]
            d["g_tilde"]=tools.read_settings(dirname+"/jastrowi.in")["g"]

            if model=={}:
                model=d
            else:
                for key in model:
                    model[key]=d[key]
            print shortname,"||",tools.dictionary_to_string(model)

def anal_all(main_dir=run_dir_default,jumps=0,bins=10,filename="energy.dat",sub_directory="energy",e_min=None,e_max=None):
    for dirname in os.listdir(main_dir):
        shortname=dirname
        dirname=main_dir+"/"+dirname
        if os.path.isdir(dirname) and check_dir_valid(dirname):
            anal.anal_energy(dirname=dirname,jumps=jumps,bins=bins,filename=filename,sub_directory=sub_directory,e_max=e_max,e_min=e_min);    
'''

def plot_table(table):
        plt.figure(figsize=(20,10))
        plt.errorbar(table[0],table[1],yerr=table[2],fmt='o')
        plt.show()

def time_step_study_set(time_steps,mean_walkers):
        for t in time_steps:
                f=input.new()
                f.set_method("dmc")
                f.set_mean_walkers(mean_walkers)
def walkers_study_set(walkers,time_step):
        for w in walkers:
                f=input.new()
                f.set_method("dmc")
                f.set_mean_walkers(w)
                f.set_time_step(time_step)

def start_all(max=None):
    i=0
    
    for f in folders:
        
        try:
            if (max!=None):
                if (i<max):
                    f.start()
                    i=i+1
                else:
                    pass
            else:
                f.start()
                
        except process_running as exc:
            exc.warn()
        except:
            print "Error in " + f.dir_path
            
def startArray():
    with open("jobArray.dat","w") as of:
        for f in folders:
            of.write(f.dir_path + "\n")
    
    
def schedule_all(scheduler,maxN=None):
    i=0
    j=0
    scheduler.wait=None
    for f in folders:
        try:
            i+=1
            
            if (maxN!=None):
                if (i>maxN):
                    j+=1
                    scheduler.wait=(datetime.datetime.now() + j*scheduler.time  )
                    i=1
            scheduler.out()
            f.start()
        except process_running as exc:
            exc.warn()
      
            
def update_status_all():
    for f in folders:
        f.update_status()
        
def stop_all():
        for f in folders:
                try:
                        f.stop()
                except no_process as exc:
                        exc.warn()
                
                        
# create a folder where to accomplish all the data analisis
def create_anal_folder(dirname="."):
        print run_dir
        try:
                if not os.path.exists(dirname):
                        os.makedirs(dirname)
        except:
                raise bad_directory

        file_list=["input.xml","launch.sh","stop.sh","qmc"]
        for file_to_copy in file_list:
                if not os.path.exists(run_dir+"/"+file_to_copy):
                        raise absent_file(file_to_copy)

        
        for file_to_copy in file_list:
                
                shutil.copy2(run_dir + "/" + file_to_copy,dirname)
                
        try_file=run_dir + "/qmc.pbs"
        if os.path.exists(try_file):
            shutil.copy2(try_file,dirname)

def filter_file(filename,reverse=False):
    global folders
    tmp=[]
    for i in range(0,len(folders)):
        f=folders[i]
        if reverse:
            if (not os.path.exists(f.dir_path + "/"+str(filename))):
                tmp.append(f)
        else:
            if (os.path.exists(f.dir_path + "/"+str(filename))):
                tmp.append(f)
    folders=tmp
    
def plot_all_energy(target="g_tilde",reset_index=False,jumps=0,linear_increment=False):
    handles=[]
    global folders
    plt.figure(figsize=(20,10))
    # plot all energies
    print "Energies"
    for f in folders:
        label=get_label(target,f)
        
        try:
            p=anal.plot_file(f.dir_path+"/energy.dat",label=label ,reset_index=reset_index,jumps=jumps,linear_increment=linear_increment)
            handles.append(p)
        except:
            print ("Warning. Cannot plot " + f.dir_path )
            
        
    plt.legend(handles=handles)
    
    handles=[]
    
def plot_all_density(dirname=".",target="g_tilde",files=["density_t.dat"],reset_index=False,jumps=0,linear_increment=False):
    
    handles=[]
    global folders
    plt.figure(figsize=(20,10))
    # plot all energies
    print "Densities"
    global folders
    for f in folders:
        try:
            label=get_label(target,f)
            for file in files:
                handle=anal.plot_file(f.dir_path+"/"+file,label=label,reset_index=reset_index,jumps=0,linear_increment=linear_increment)
                handles.append(handle)
        except IOError:
            print "Error in input file in directory" + f.dir_path
        except :
            print "Error in " + f.dir_path
    plt.legend(handles=handles)
    handles=[]
            

def plot_all_pair_correlation(dirname=".",target="g_tilde"):
    
    handles=[]
    global folders
    plt.figure(figsize=(20,10))
    global folders
    # plot all energies
    print "Pair correlations"
    for f in folders:
        label=get_label(target,f)
        handles.append(anal.plot_file(f.dir_path+"/pair_correlation_t.dat",label=label,error=True ))
        
    plt.legend(handles=handles)
    plt.axhline(y=1, xmin=0,hold=True)
    
    handles=[]

def anal_all(dirname=".",jumps=0,bins=10,cut_low=1,cut_heigh=None):
        anal_all_energy(jumps=jumps,bins=bins)
        anal_all_effective_mass(cut_low=cut_low,cut_heigh=cut_heigh)

def plot_all_effective_mass(dirname=".",target="g_tilde",label="winding_number",factor=1,fit_lines=True):
    
    handles=[]
    winding_file=label+"_t.dat"
    plt.figure(figsize=(20,10))
    #plot the winding number and the effective mass
    print "Effective mass"
    for f in folders:
        label=get_label(target,f)
        handles.append(anal.plot_effective_mass(f.dir_path+"/" + winding_file,settings_file=f.dir_path + "/input.xml",param_file=f.dir_path+"/effective_mass_t.dat",label=label,factor=factor,fit_lines=fit_lines))
        
        # plot the fit function if present
        
    plt.legend(handles=handles)
    
    handles=[]

    
def anal_all_energy(jumps=0,bins=10):
    for f in folders:
        
        f.anal_scalar("energy.dat",jumps,bins)
        
def anal_all_effective_mass(cut_low=1,cut_heigh=None,n_cuts=10,factor=1):
    for f in folders:
        f.anal_effective_mass(cut_low=cut_low,cut_heigh=cut_heigh,n_cuts=n_cuts,factor=factor)
        print f.get_g()

def get_energy_table(target="g_tilde"):
    
    values=[]
    # cycle over all possible folders
    for f in folders:
        
        label=get_values(target,f)
        
        if (os.path.isfile(f.dir_path+"/energy_res_anal.dat")):
            value=tools.get_scalar_value(f.dir_path+"/energy_res_anal.dat")
            print value
            #if label !="inf":
            
            values.append(np.hstack( (tools.to_float(label),value) ))
    return values

def get_variance_tab(tab):
    tabv=tab
    tabv[:,1]=tab[:,3]
    return tabv[:,0:2]
    
    
def get_effective_mass_table(dirname=".",target="g_tilde"):
    #dirname=os.path.abspath(dirname);
    #if (dirname[-1]=="/"):
    #     dirname=dirname[:-1]
    #  analize the energy
    #input.run_dir=dirname
    #input.scan();
    global folders
    values=[]
    # cycle over all possible folders
    for f in folders:
        label=get_values(target,f)
        if (os.path.isfile(f.dir_path+"/effective_mass_t.dat")):
            value=anal.get_effective_mass(f.dir_path+"/effective_mass_t.dat")
            
            values.append([tools.to_float(label),float(value[0]),float(value[1])])
    return values

def get_folder(target,src_target):
    
    for f in folders:
            
            if match_target(f,target,src_target):
                
                return f
    return None

def filter_folders(target,src_target,eps=1e-4):
    
    global folders
    fs=[]
    for f in folders:
        
        
        if match_target(f,target,src_target,eps=eps):
            fs.append(f)
    folders=fs

def deleteFolder(dirname):
    global folders

    for i in range(0,len(folders)):
        
        if folders[i].dir_name==str(dirname):
            del folders[i]

            
# builds a table from all of them
def build_table(target,filename="energy.dat"):
    global folders
    tab=[]
    for f in folders:
        tmp=anal.anal_energy(filename=f.dir_name+"/"+filename,jumps=200)
        tab.append([float(get_value(target,f)),tmp[0],tmp[1],tmp[2]])
    return np.array(tab)

def structure_factor_q_step(f,label):

        mes=f.get_measure_by_label(label)
        if mes==None:
            print "No measurement Found ?"
            return None

        l=float(f.get_lBox())
        # decide the number of bins value
        
        bins=int(mes.attrib["bins"])

        # decide the max value
        
        if "max" in mes.attrib:
            maxX=float(mes.attrib["max"])
        else:
            maxX=bins*2*pi/l
        
        deltaQ=floor(maxX/((2*pi/l)*bins))*(2*pi/l)
        if (deltaQ==0):
            deltaQ=2*pi/l
        return deltaQ
    
# returns the structure factor using the usual procedure
def analStructureFactor(f,label,bins=None,makePlot=False,jumps=0):
    res=anal.analVectorHistory(f.dir_path+"/"+label+".dat",bins=bins,jumps=jumps,makePlot=makePlot)
    df=pd.DataFrame({"x":res[:,0],"value":res[:,1],"error":res[:,2]})
    qStep=structure_factor_q_step(f,label)
    df["x"]= df["x"]*qStep + 2*pi/float(f.get_lBox()) 
    return df

def getDensity(f,label,normalize=1,makePlot=False):
    densityMeasure=f.get_measure_by_label(label)
    df=pd.read_csv(f.dir_path + "/" +label + "_t.dat",delim_whitespace=True,names=["x","value","delta","conv"],header=None,index_col=False)
    lBox=float(f.get_lBox())
    df["x"]=df["x"]/len(df["x"])*lBox - lBox/2.
    if normalize is not None:
        step=abs(df["x"][1]-df["x"][0])
        A=np.sum(df["value"])*step
        if normalize=="N":
            normalize=float(f.get_n_particles(index=0)) + float(f,get_n_particles(index=1))
        df["value"]/=A
        df["value"]*=normalize

    if makePlot:
        plt.errorbar(df["x"],df["value"],df["delta"],fmt="or")
        
    return df
    

def getPairCorrelation(f,label,makePlot=False):
    densityMeasure=f.get_measure_by_label(label)
    setA=float(densityMeasure.attrib["setA"])
    setB=float(densityMeasure.attrib["setB"])
    df=pd.read_csv(f.dir_path + "/" +label + "_t.dat",delim_whitespace=True,names=["x","value","delta","conv"],header=None,index_col=False)
    lBox=float(f.get_lBox())
    n=float(f.get_particles(index=0)) + float(f.get_particles(index=1))
    
    df["x"]=df["x"]/len(df["x"])*lBox/2.
    
    df["x"]*=(n/lBox)

    
    nA=float(f.get_particles(index=setA))
    nB=float(f.get_particles(index=setB))
    if (setA == setB):
        
        normalize=nA*(nA-1)*1./n
    else:
        normalize=nA*nB*1./n
        
    
    if normalize is not None:
        step=abs(df["x"][1]-df["x"][0])
        A=np.sum(df["value"])*step
        if normalize=="N":
            normalize=float(f.get_n_particles(index=0)) + float(f,get_n_particles(index=1))
        df["value"]/=A
        df["value"]*=normalize

        df["delta"]=df["delta"]/A*normalize

    if makePlot:
        plt.errorbar(df["x"],df["value"],df["delta"],fmt="or")
        
    return df
    
    
    
