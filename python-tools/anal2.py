#!/home/luca/software/anaconda/bin/python
import os
import subprocess
import glob
import tools

class anal:
    def __init__(self):
        self.inputfile='inputfile'
        self.outputfile='out.log'
        self.update=False
    def processOutputFile(dataframe):
        return dataframe

    def processInputFile(dataFrame):
        return dataframe
'''    
    def readoutputfile(self,filename,header=[]):
        f=open(filename)
        data={}
        if (not header):
            #assume the first row as the header
            header=f.readline().split()
        #set the dictionary
        for col in header:
            data[col]=[]
        #add values to the columns of the dictionary
        for line in f:
            row=line.strip().split()
            i=0
            for col in header:
                #check if it is a number
                try:
                    float(row[i])
                except:
                    pass
                else:
                    row[i]=float(row[i])
                data[col].append(row[i])
                i=i+1
        return data
'''
'''
    def readinputfile(self,filename):
        data={}
        for line in open(filename):
            if line[0] != '#':
                row=line.strip().split('=')
                try:
                    float(row[1])
                except:
                    pass
                else:
                    row[1]=float(row[1])
                data[row[0]]=[row[1]]
        return data
'''
    def merge(self,data1,data2):
        data={}
        keys1=data1.keys()
        keys2=data2.keys()
        n2=len(data2[keys2[0]])
        n1=len(data1[keys1[0]])
        #initialize
        for key in keys1+keys2:
            data[key]=[]
        for i in range(0,n1):
            for key in keys2:
                data[key]=data[key]+data2[key]
            for key in keys1:
                data[key]=data[key]+[ data[key][i] ]*n2
        return data
    def lookup(self):
        #buid list of all directories starting with run
        dirs_files=glob.glob('run*')
        dirs=[]
        datalist=[]
        for file in dirs_files:
            if os.path.isdir(file):
                dirs.append(file)
        for dir in dirs:
            print "analize directory "+dir+" ..."
            #make the input dataframe
            inputdata=pd.DataFrame(self.readinputfile(dir+'/'+self.inputfile))
            #select input fields
            
            inputdata.loc[:,'dir']=dir
            #make the output dataframe(select last column)
            outputdata=pd.DataFrame(self.readoutputfile(dir+'/'+self.outputfile))
            #purge output data
            outputdata=self.processOutputFile(outputdata)
            outputdata.loc[:,'dir']=dir
            #############################make statistics on OPpopulation.dat
           
            
            if (not os.path.isfile(dir+"/OPpopulation.sta")) or ( self.update):
                subprocess.call("sed -n '0~1000p' "+dir+"/OPpopulation.dat > "+dir+"/OPpopulation.p;",shell=True)
                subprocess.call("./checkEquilibrium.py "+dir+"/"+"OPpopulation.p "+">"+dir+"/OPpopulation.sta",shell=True)
            OPpopulationdata=pd.DataFrame(self.readoutputfile(dir+'/'+'OPpopulation.sta'))
            OPpopulationdata.loc[:,'dir']=dir
            #merge OPpopulation data mith input data
            inputdata=pd.merge(inputdata,OPpopulationdata,on='dir',how='inner')
            ############################################################
            ##########################################make the distribution
            distributionsdata=[]
            for file in outputdata['DISTRIBUTION']:
                disfile=dir+"/"+file[:-3]+'dis'
                if (not os.path.isfile(disfile)) or ( self.update):
                    subprocess.call("./anal "+dir+"/"+file+" > "+disfile,shell=True)
                distributiondata=pd.DataFrame(self.readoutputfile(disfile))
                distributiondata.loc[:,'FILE']=file
                #print distributiondata
                distributionsdata.append(distributiondata)
            #merge distribution with ouput
            outputdata=pd.merge(outputdata,pd.concat(distributionsdata),left_on='DISTRIBUTION',right_on='FILE')
            ##################################################################
            #merge input and output
            data=pd.merge(inputdata,outputdata,on='dir',how='inner')
            data=data[['Temp','HFlux','TIME','NUMATOMS','a_fin','cov_fin','cov_fin_std','ClusterSize','NumClusters']]
            #data=data[['Temp','HFlux','TIME','NUMATOMS','a_fin','cov_fin','cov_fin_std']]
            datalist.append(data)
        
        return pd.concat(datalist,ignore_index=True)
    

