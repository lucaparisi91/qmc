import pandas as pd
import manageBeta
import manage
import inputFolders
import matplotlib.pylab as plt

inputFolders.set_anal_folder();
inputFolders.filter_folders("n index 0",20)
df=manage.gatherDataOptAlpha(jumps=0,last=False,makePlot=False)
df.to_csv("alphaOptN40.dat")
plt.show()
