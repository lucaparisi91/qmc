import re
import os
import shutil
import hashlib
import tools
import jastrow
import os.path
import random
import anal
import xml
import numpy as np
import  xml.etree.ElementTree as ET
import subprocess
import exceptions

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
def new():
    global folders
    folders.append(folder())
    folders[len(folders)-1].create()
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
    def create(self):
        #creates a unique dirname(to use as a label)
        dirname=str(random.randint(0,1000000))
        md5=hashlib.md5()
        md5.update(dirname)
        dirname=md5.hexdigest()[-8:]
        dirname=re.sub(" ","_",dirname)
        self.dir_name=dirname
        dirname=run_dir+"/"+dirname
        self.dir_path=dirname
        #creates the fold if necessary
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        # copy default files in global variables
        #print run_dir
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
                    p=subprocess.Popen(run_dir + "/launch.sh " + self.dir_path,shell=True,stdout=subprocess.PIPE)
            #p.wait()
                    self.pid=p.stdout.read()
                    self.running=True
                    # saves info to a file
                    self.save_process_status()
            else:
                    # exception: the file is still running
                    raise process_running
            
    def save_process_status(self):
            d={"pid":self.pid,"running":self.running}
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
        list_init=["save.xml","walkers.xml"]
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
        
    def set_winding_number(self,bins=1000,time=1):
            steps=int( float(time)/float( self.get_time_step() ) ) 
            skip=steps/int(bins) -1
            if (skip<0):
                    skip=0
                    bins=steps
            
            self.root.find("measures").find("winding_number").attrib["bins"]=str(bins)
            self.root.find("measures").find("winding_number").attrib["skip"]=str(skip)
            
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
    
    def plot_file(self,filename,error=True,label="",reset_index=False,jumps=0):
            anal.plot_file(self.dir_path+"/"+filename,error,label,reset_index,jumps=jumps)
            
            return self
    def plot_effective_mass(self,filename="winding_number_t.dat",label=""):
            
            anal.plot_effective_mass(self.dir_path+"/"+filename,settings_file=self.dir_path + "/input.xml",param_file=self.dir_path+"/effective_mass_t.dat",label=label)

    def anal_effective_mass(self,data="winding_number_t.dat",settings="input.xml",out="effective_mass_t.dat",cut_low=0,cut_heigh=None,n_cuts=10):            
            anal.anal_effective_mass(self.dir_path+"/" + data,settings_file=(self.dir_path+"/"+settings),out_file=(self.dir_path+"/"+out),cut_low=cut_low,cut_heigh=cut_heigh,n_cuts=n_cuts)

    def removeFolder(self):
        
        if not os.path.exists(run_dir+"/backup"):
            os.makedirs(run_dir+"/backup")
            
        shutil.move(self.dir_path,run_dir+"/backup/"+self.dir_name)
            
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

def start_all():
        for f in folders:
                try:
                        f.start()
                        
                except process_running as exc:
                        exc.warn()
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

