#!/home/luca/software/anaconda/bin/python
import numpy as np
import anal
import sys
# import an analizing system for the particles in the system
#print sys.argv
if (sys.argv[5]=="None"):
	e_min=None
else:
	e_min=float(sys.argv[5])
if (sys.argv[6]=="None"):
        e_max=None
else:
        e_max=float(sys.argv[6])

anal.anal_energy(sys.argv[1],filename=sys.argv[2],jumps=int(sys.argv[3]),bins=int(sys.argv[4]),e_min=e_min,e_max=e_max)
