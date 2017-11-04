#!/home/luca/software/anaconda/bin/python
import scipy as sp
from scipy import optimize
import numpy as np
import string
import os
import glob
import shutil
import subprocess
import threading
import jastrow
import tools
import signal
import threading
from multiprocessing import Process
import pandas as pd
from math import *
import re
import matplotlib.pylab as plt
# classes to analize a job
import anal

def write_dictionary_as_input(d,filename):
    f=open(filename,"w+")
    for key in d:
        f.write(key+"="+ str(d[key])+"\n")
class func_thread(threading.Thread):
    def __init__(self,target,*args):
        self._target=target
        self._args=args
        threading.Thread.__init__(self)
        
    def run(self):
        self._target(*self._args)
        
        

class settings:
    #name of the file containing the settings
    
    def __init__(self,filename,dirname):
        
        self.filename=tools.get_name(filename)
        self.variables={}
        self.dirname=dirname
        self.comments={}
        self.order_variables=[]
    def import_settings(self,filename):
        
        comment=""
        
        f=open(filename)
        for line in f :
            line=line.strip()
            if (line[0] == "#") :
                comment=comment+line
                continue
            
            line=line.split('=')
            
            self.variables[line[0]]=line[1]
            self.order_variables.append(line[0])
            self.comments[line[0]]=comment
            comment=""
        
        
        filename=tools.get_name(filename)
            
        self.filename=filename
        f.close()
    
    def export_settings(self,filename):
        
        f=open(filename,"w+")
        for key in self.order_variables:
            
            if (string.strip(self.comments[key]) != ""):
                f.write(self.comments[key]+"\n")
            f.write(key+"="+self.variables[key]+"\n")
        f.close()
    def export_data(self):
        self.export_settings(self.dirname+"/"+self.filename)
    def import_data(self):
        self.import_settings(self.dirname+"/"+self.filename)
    def update_setting(self,key,value):
        if key in self.variables.keys():
            self.variables[key]=str(value)
        else:
            print "The key" + key+ " was not found"
        return self
    def get_setting(self,key):
        if key in self.variables.keys():
            return self.variables[key]
        else:
            print "The key was not found."
            return None
    def check_setting(key):
        if key in self.variables.keys():
            return True
        else:
            return False
    def show(self):
        msg="----"+self.filename+"\n"
        for key in self.variables:
            msg=msg+key+"="+self.variables[key]+","
        return msg
        
class run:
    
    def __init__(self,key):
        self.last_job_id=0
        self.files=[]
        self.job_list=[]
        self.setting_list=[]
        self.key=int(key)
        self.job_del_list=[]
        # name of the directory 
        self.dirname="run_" + str(key)
        # create the directory if it does not exists
        if not os.path.exists(self.dirname):
                os.makedirs(self.dirname)
     #import settings in a file
    def add_job(self,key=None):
        if (key == None):
            self.last_job_id=self.last_job_id+1
            key=self.last_job_id
        else:
            if (key > self.last_job_id):
                self.last_job_id=int(key)
        j=job(key,self.dirname)
        self.job_list.append(j)
        return j
    
    def add_qmc(self):
        j=self.add_job().command="qmc"
        return j
    def add_del_job(self,key=None):
        if (key == None):
            self.last_job_id=self.last_job_id+1
            key=self.last_job_id
        else:
            if (key > self.last_job_id):
                self.last_job_id=key
        j=job(key,self.dirname)
        self.job_del_list.append(j)
        return j
    
    def show_jobs(self):
        msg=""
        for job in self.job_list:
            msg=msg + job.show() + "\n"
        return msg
    def get_job(self,key):
        for job in self.job_list:
            
            if job.key == key:
                
                return job
            
        print "No job found"
        return None
    
    def get_db(self):
        d={}
        for setting in self.setting_list:
            d.update(setting.variables)
        d["run_key"]=self.key
        j=self.get_job(1)
        if j==None:
            d["run_command"]=""
        else:
            d["run_command"]=self.get_job(1).command
        if os.path.exists(self.dirname+"/energy/status.log"):
            d.update(tools.read_settings(self.dirname+"/energy/status.log"))
        else:
            d["mean"]=0
            d["error_mean"]=0
        #d["id"]=self.get_job(1).job_id
        if len(self.job_list) !=0:
            d["job_pid"]=self.get_job(1).job_pid
            d["job_status"]=self.get_job(1).job_status    
        return pd.DataFrame(d,index=[0])
    
    def add_file(self,filename,overwrite=True):

        file_exist=os.path.exists(self.dirname+"/"+filename)
        if overwrite or (not file_exist) :
            shutil.copy(filename,self.dirname+"/"+filename)
            if not (filename in self.files):
                self.files.append(filename)
        return self
    def add_settings(self,filename):
        s=settings(filename,self.dirname)
        if os.path.exists(filename):
            s.import_settings(filename)
        self.setting_list.append(s)
    def default(self):
        #self.import_settings("qmc.in")
        #self.import_settings("system.in")
        #self.import_settings("measures.in")
        #self.import_settings("jastrow.in")
        #self.files=["VMC","DMC"]
        pass
    
    def get_settings(self,filename):
        for setting in self.setting_list:
            if setting.filename==filename:
                return setting
        print "setting not found"
        return None
    
    def show(self):
        msg="----------------key="+str(self.key)+"\n"
        msg=msg + "dir: " + self.dirname + "\n"
        msg=msg+"files: "
        for file_in_run in self.files:
            msg=msg+file_in_run+ ","
        msg=msg+"\n"
        for job in self.job_list:
            msg=msg+job.show()
        
        msg=msg+"Setting files:"
        for setting in self.setting_list:
            msg=msg+ setting.filename + ","
        
        msg=msg+ "\n------Settings:"
        for setting in self.setting_list:
            msg=msg+setting.show()+"\n"
        msg=msg+ "\n----------------"
        return msg
        
    def export_data(self,dirname="."):
        
        filename=self.dirname + "/run.log"
        f=open(filename,"w+")
        f.write("files=>")
        for file_to_copy in self.files:
            f.write(file_to_copy + ",")
        f.write("\n")
        f.write("settings=>")
        for setting in self.setting_list:
            f.write(setting.filename + ",")
        f.write("\n")

        f.write("jobs=>")
        for job in self.job_list:
            f.write(str(job.key) + ",")
            job.export_data()
        f.write("\n")

        f.write("jobs_del=>")
        for job in self.job_del_list:
            f.write(str(job.key) + ",")
            job.export_data()
        f.write("\n")
        
        f.close()

        for setting in self.setting_list:
            setting.export_data()
    
    def import_data(self,dirname="."):
        d=tools.read_settings(self.dirname+"/run.log")
        self.files=d["files"]
        for filename in d["settings"]:
            s=settings(filename,self.dirname)
            s.import_data()
            self.setting_list.append(s)
        for job_key in d["jobs"]:
            j=self.add_job(int(job_key))
            j.import_data()
        for job_key in d["jobs_del"]:
            j=self.add_del_job(int(job_key))
            j.import_data()
    # make some premilinary operations on inputfiles
    def make_jastrows(self):
        # evaluates the jastrow paramaters based on inputfile for bosons
        j1=jastrow.jastrow()
        d={}
        d["g"]=self.get("g")
        d["l_box"]=self.get("n_particles")
        d["cut_off"]=self.get("cut_off")
        j1.import_data_dict(d)
        j1.process()
        
        j1.print_jastrows(self.dirname)
        # impurity boson jastrow
        j2=jastrow.jastrow()
        d={}
        d["g"]=self.get("g_tilde")
        d["l_box"]=self.get("n_particles")
        d["cut_off"]=self.get("cut_off")
        j2.import_data_dict(d)
        j2.process()
        j2.print_jastrows(self.dirname,"_i")
        
        return self
    
    def anal_energy(self,jumps=0,bins=10,filename="e.dat",e_min=None,e_max=None):
        if self.get_job(1).command=="VMC":
            filename="e.dat"
        
        p=subprocess.Popen(["./anal_directory.py",self.dirname,filename,str(jumps),str(bins),str(e_min),str(e_max)])
        
        #anal.anal_energy(self.dirname,filename=filename,jumps=jumps,bins=bins)
        #t1.start()
        #t1.join()
    def plot_energy(self,begin=0,alpha=1,end=0,filename="e.dat"):
        
        anal.plot_single(self.dirname+"/"+filename,begin=begin,alpha=alpha,end=end)
        
    def launch(self):
        if len(self.job_list)==0:
            print "No jobs on the list"
        self.get_job(1).launch()
        self.export_data()
    def kill(self):
        if len(self.job_list)==0:
            print "No jobs on the list"
        self.get_job(1).kill()
        self.export_data()
    def check(self):
        if len(self.job_list)==0:
            print "No jobs on the list"
        self.get_job(1).check()
        self.export_data()
    def clean(self):
        #remove all job files not listed
        pass
    def delete_job(self,key):
        for i in range(0,len(self.job_list)):
            if self.job_list[i].key == key:
                self.job_del_list.append(self.job_list[i])
                self.job_list.remove(self.job_list[i])

    def set(self,setting_name,setting_value):
        setting_name=str(setting_name)
        setting_value=str(setting_value)
        found=0
        for setting in self.setting_list:
            if setting_name in setting.variables:
                setting.update_setting(setting_name,setting_value)
                found=found+1
        if (found == 0):
            print "Setting "+ setting_name + " was not found." 
        
        return self

    def get(self,setting_name):
        setting_name=str(setting_name)
        for setting in self.setting_list:
            if setting_name in setting.variables.keys():
                return setting.variables[setting_name]
        return None
        
class job:
    
    def __init__(self,key,dirname):
        self.key=int(key)
        self.job_id=0
        self.job_pid=0
        self.job_status=0
        self.job_priority=0
        self.command=""
        self.dirname=dirname
        self.job_output="out.log"
    # start analizing the file containting energy    
    
    
        #write_dictionary_as_input(j.get_output_parameters(),self.dirname+"/jastrow.in")
    def check(self):
        # check if the programming is still running
        if int(self.job_status)==1:
            if not tools.check_pid(self.job_pid):
                self.job_status=2
            
    def launch(self):
        if (self.job_status == 1):
            print "Job already active"
            return self
        f=open(self.job_output,"w+")
        dirname=os.path.abspath(self.dirname)
        command=dirname+"/"+self.command
        print command
        try:
            process=subprocess.Popen([command],shell=False,cwd=self.dirname,stdout=f)
        except:
            print "The process did not start"
            
        else:
            try:
                self.job_pid=process.pid
            except:
                print 'could not find the PID'
            else:
                self.job_status=1
        f.close()
        return self
    def set_command(self,command):
        self.command=command
        
    def kill(self,signal_to_send=signal.SIGTERM):
        try:
            os.kill(self.job_pid, signal_to_send)
        except:
            print "Could not kill the process."
            return self
        
        if (signal_to_send == signal.SIGTERM):
            self.job_status=2 # killed status
        else:
            self.job_status=3 # unkown status
        return self
            
    def export_data(self):
        # make log
        
        dirname=self.dirname + "/jobs"
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        
        f=open(dirname+"/job_"+str(self.key)+".log","w+")
        f.write("job_priority="+str(self.job_priority)+"\n")
        f.write("job_id="+str(self.job_id)+"\n")
        f.write("job_status="+str(self.job_status)+"\n")
        f.write("job_pid="+str(self.job_pid)+"\n")
        f.write("command="+str(self.command)+"\n")
        f.close()

        return self
        
    def import_data(self):
        
        filename=self.dirname+"/jobs/job_"+str(self.key)+".log"
        if not os.path.exists(filename):
            print "No job log found at "+ filename
            return self
        
        d=tools.read_settings(filename)
        self.job_priority=int(d["job_priority"])
        self.job_id=int(d["job_id"])
        self.job_pid=int(d["job_pid"])
        self.job_status=int(d["job_status"])
        self.command=d["command"]
    def show(self):
        msg="Job "+str(self.key)+"-> " + "Status: "+str(self.job_status) + ",Pid: "+str(self.job_pid) + ",Command: " + str(self.command) + "\n"
        return msg
            
class runs():
    
    
    def __init__(self):
        self.last_job_id=0
        self.run_list=[]
        self.run_del_list=[]
    def get(self,key):
        for run_ob in self.run_list:
            if run_ob.key==key:
                return run_ob
    def get_db(self):
        dbs=[]
        for r in self.run_list:
            dbs.append(r.get_db())
        return pd.concat(dbs)
            
    def add_run(self,key=None):
        if (key == None):
            self.last_job_id=self.last_job_id+1
            key=self.last_job_id
        else:
            if (key > self.last_job_id):
                self.last_job_id=key
        r=run(key)
        self.run_list.append(r)
        return r
    
    def add_del_run(self,key=None):
        if (key == None):
            self.last_job_id=self.last_job_id+1
            key=self.last_job_id
        else:
            if (key > self.last_job_id):
                self.last_job_id=key
        r=run(key)
        self.run_del_list.append(r)
        return r
    
    def export_data(self,dirname='.'):
        f=open("runs.log","w+")
        f.write("keys=>")
        for run_ob in self.run_list:
            f.write(str(run_ob.key)+",")
            run_ob.export_data()
        

        f.write("\nkeys_del=>")
        for run_ob in self.run_del_list:
            f.write(str(run_ob.key)+",")
            run_ob.export_data()
        f.close()
        
    def import_data(self,dirname='.'):
        if not os.path.exists("runs.log"):
            print "No runs found."
            return
        d=tools.read_settings("runs.log")
        for key in d["keys"]:
            self.add_run(int(key)).import_data()
        
        for key in d["keys_del"]:
            self.add_del_run(int(key)).import_data()
    def delete(self,key):
        for i in range(0,len(self.run_list)):
            if self.run_list[i].key == key:
                # mark a run for deletion
                self.run_del_list.append(self.run_list[i])
                self.run_list.remove(self.run_list[i])
                return self
        print "Not found"
        
        self.export_data()
        
    # delete folders marked for deletion
    def clean(self):
        for element in self.run_del_list:
            
            shutil.rmtree(element.dirname)
            self.run_del_list.remove(element)
            print "removed directory "+ element.dirname    
        self.export_data    
        return self
    
    
    def show(self):
        for run_ob in self.run_list:
            print run_ob.show()
    
    
    def add_qmc(self):
        r=self.add_run()
        r.add_file("qmc")
        r.add_settings("qmc.in")
        r.add_settings("measures.in")
        l_box=int(r.get("n_particles"))
        #l_box=float(r.get_settings("qmc.in").get_setting("l_box"))
        # set the default maximum value of the correlation function
        r.get_settings("measures.in").update_setting("g_r_max",l_box/2)
        r.add_settings("jastrow.in")
        #r.add_settings("jastrow_i.in")
        r.get_settings("jastrow.in").update_setting("cut_off",l_box/2)
        #r.get_settings("jastrow_i.in").update_setting("cut_off",l_box/2)
        r.make_jastrows()
        self.export_data()
        return r
    
    
ps=runs()
ps.import_data()
#ps.delete(3)
#ps.delete(1)
#ps.add_qmc()

#ps.clean()
#ps.get(1).preprocess()
#ps.get(1).get_job(1).check()
#ps.export_data()
ps.show()
