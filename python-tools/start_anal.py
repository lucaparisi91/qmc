#!/home/luca/software/anaconda/bin/python

import anal
from run import *

for run_ob in ps.run_list:
    #if (run_ob.key ==1):
    #    continue
    anal.anal_energy(run_ob.dirname)
