import pandas as pd
import manageBeta
import manage
import inputFolders
import matplotlib.pylab as plt

#gather VMC data
#inputFolders.set_anal_folder();
#inputFolders.set_anal_folder();inputFolders.filter_folders("method","svmc")
#inputFolders.set_anal_folder();inputFolders.filter_folders("method","vmc");inputFolders.filter_folders("landa_tilde",0.1);inputFolders.filter_folders("n index 0",10);

#dfV=manage.gatherData(jumps=100)
#dfV.to_csv("dataFinalV2.dat",sep=" ")

# gather DMC data
#inputFolders.set_anal_folder();inputFolders.filter_folders("method","dmc");inputFolders.filter_folders("landa_tilde",100);inputFolders.filter_folders("n index 0",10);inputFolders.filter_folders("trap center 0",-0.05)

def filterFutureWalkers(neg=False):
    tmps=[]
    for f in inputFolders.folders:
        isFW=f.get_measure_by_label("center_of_mass_difference").attrib["futureWalkers"]
        if neg:
            if isFW=="false":
                tmps.append(f)
        else:
            if isFW=="true":
                tmps.append(f)
                
    inputFolders.folders=tmps

inputFolders.set_anal_folder();inputFolders.filter_folders("method","dmc");
filterFutureWalkers(neg=True)

#inputFolders.set_anal_folder();inputFolders.filter_folders("method","dmc");inputFolders.filter_folders("landa_tilde",0.1);inputFolders.filter_folders("n index 0",10);
dfD=manage.gatherData(jumps=500)
dfD.to_csv("dataFinalLandaTilde2.dat",sep=" ")
