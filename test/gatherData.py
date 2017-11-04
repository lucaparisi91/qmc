import manage
import pandas as pd
import inputFolders

inputFolders.set_anal_folder()
inputFolders.filter_folders("method","svmc")
df=manage.gatherData(jumps=20.)
df.to_csv("dataBetaV.dat",sep=" ")

