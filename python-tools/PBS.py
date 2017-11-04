import string
import datetime
from math import *

class noWallTime:
    pass
class noWaitTime:
    pass

class PBSScheduler:
    def __init__(self):
        self.que="allq"
        self.output="qmc.out"
        self.name="qmc"
        self.wait=None
        self.err="qmc.e"
        self.nodes=1.
        self.outFile="qmc.pbs"
        self.ppn=1.
        self.time=None
    def stringWallTime(self):
        if self.time != None:
            seconds=self.time.seconds
            h=int(floor(seconds/3600.))
            m=int(floor((seconds - h*3600)/60.))
            s=seconds - h*3600 - m*60.
            
            return "%0.2d:%0.2d:%0.2d" % (h,m,s)
        
        else:
            raise noWallTime()
    
    def stringWaitTime(self):
        if self.wait != None:
            return self.wait.strftime("%Y%m%d%H%M.%S")
        else:
            raise noWaitTime()
        
    def setWallTime(self,time):
        self.time=time
        
    def setCommand(self,command):
        self.command=command
        
    def setNodes(self,n):
        self.nodes=n

    def setPPN(self,n):
        self.ppn=n
        
    def process(self):
        script="#!/bin/bash\n"
        script+="#PBS -l nodes="+str(int(self.nodes))+":ppn="+str(int(self.ppn))
        if self.time!=None:
            script+=",walltime="+self.stringWallTime()+"\n"
        else:
            script+="\n"

        script+="#PBS -q " + str(self.que) + "\n"

        script+="#PBS -o " + str(self.output) + "\n"

        script+="#PBS -N " + str(self.name) + "\n"
        
        script+="#PBS -e " + str(self.err) + "\n"
        
        if self.wait != None:
            script+="#PBS -a " + self.stringWaitTime() + "\n"
            
        script+="cd $PBS_O_WORKDIR\n"

        script+=self.command + "\n"
        
        return script
    
    # output script to file
    def out(self):
        f=open(self.outFile,"w+")
        f.write(self.process())
        f.flush()
        f.close()
    
    
        
