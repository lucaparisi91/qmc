import pandas as pd
import manageBeta
import manage
import inputFolders
import matplotlib.pylab as plt

#gather VMC data
#inputFolders.set_anal_folder();
inputFolders.set_anal_folder();inputFolders.filter_folders("method","svmc")
#inputFolders.set_anal_folder();inputFolders.filter_folders("method","vmc");inputFolders.filter_folders("landa_tilde",0.1);inputFolders.filter_folders("n index 0",10);

dfV=manage.gatherData(jumps=10)
dfV.to_csv("dataFinalV2.dat",sep=" ")

# gather DMC data
#inputFolders.set_anal_folder();inputFolders.filter_folders("method","dmc");inputFolders.filter_folders("landa_tilde",100);inputFolders.filter_folders("n index 0",10);inputFolders.filter_folders("trap center 0",-0.05)

#inputFolders.set_anal_folder();inputFolders.filter_folders("method","dmc");

#inputFolders.set_anal_folder();inputFolders.filter_folders("method","dmc");inputFolders.filter_folders("landa_tilde",0.1);inputFolders.filter_folders("n index 0",10);
#dfD=manage.gatherData(jumps=0)
#dfD.to_csv("dataFinalLandaTildeD2.dat",sep=" ")
