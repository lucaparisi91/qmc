#!/home/luca/software/anaconda/bin/python
# [E]=h_bar**2/m
# [L]=1/n
# binding energy for g_B-> infinity for all g_I(coupling of the impurity)
import tools
import matplotlib.pylab as plt
from math import *

def TG_binding_energy(g_tilde):
    #g->g_tilde
    kf=pi/2
    return 1/pi*g_tilde**2*(kf/g_tilde-atan(kf/g_tilde))
def perturbative_binding_energy(g_tilde,gamma):
    # g-> 0,g_tilde -> 0, infinite mass
    return g_tilde**2/2*gamma*(-1/2)
def TG_inf_binding_energy(g_tilde):
    ef=pi**2/4
    return 2*ef*g_tilde**2/pi**3*(pi/g_tilde + (pi/g_tilde)**2*atan(g_tilde/pi)-atan(pi/g_tilde))

def TG_mob_binding_energy(g,N):
    kf=(N+0.5)*pi/N
    x=g/(2*kf)
    return kf**2/pi*(  x + atan(x) - x**2*(pi/2-atan(x)))

def TG_mass(g):
    x=2*kf/g
    return 2/pi*(atan(x))**2/(atan(x)-x/(1+x**2))cd 
