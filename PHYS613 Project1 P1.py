"""
Created on Wed Oct 02 18:53:32 2013
PHYS 613, Project 1
Nick Crump
"""

# Project 1: Minimum Energy Molecules
"""
Finding the optimal atomic positions that minimize the interatomic potential 
energy of Potassium molecules (K). 
"""

import numpy as np
from math import sqrt,exp
import matplotlib.pyplot as plt


# Part 1
# define function for interatomic potential Ucoh = Urep - Uel
# returns value of function for a constant separation distance d
#*******************************************************************
def P(n,d):
    # constants for Potassium molecule (K)
    eps = (1.51/1000)*(27.21138386/2.0)  # (eV)
    gam = (19.3/1000)*(27.21138386/2.0)  # (eV)
    r0 = (8.253*0.52917721)              # (angstroms)
    p = 10.58
    q = 1.34
    
    # terms for iteration
    Urep = 0
    Uel = 0
    
    # iterate through double sum of energy function
    for i in range(n):
        Uel0 = 0
        
        for j in range(n):
            if j != i:
                # calculate repulsive and attractive energy terms
                Urep = Urep + eps*(exp(-p*((d/r0)-1.0)))
                Uel0 = Uel0 + (gam**2.0)*(exp(-2.0*q*((d/r0)-1.0)))
        
        # take care of sqrt in sum over j
        Uel = Uel + sqrt(Uel0)
    
    # calculate total cohesive potential energy
    Ucoh = Urep - Uel

    return Ucoh
#*******************************************************************    

# get potential energy for 3 molecules as a function of distance
PE = []
dist = np.arange(3.2,11,0.01)

# loop to get energies
for i in dist:
    PEi = P(3,i)
    PE.append(PEi)
        
# plot potential energy vs separation distance
plt.figure()
plt.plot([3.2,11],[0,0],'k--')
plt.plot(dist,PE,label='N=3 atoms')
plt.xlabel('Atomic Separation (angstroms)')
plt.ylabel('Interatomic Potential (eV)')
plt.legend()