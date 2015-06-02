"""
Created on Sat Sep 21 14:03:02 2013
PHYS 613, Assignment 4
Nick Crump
"""

# Problem 1 (EBB4)
"""
Finding best fit parameters using a non-linear least squares method to 
molecule data with a fit model following the Morse potential. 
"""

from LeastSquares import gaussNewton3
from math import exp
import numpy as np
import matplotlib.pyplot as plt


# import sample data set "P1Molecule.txt"
# NOTE: Positions in angstroms, Energy in eV
pos1, pos2, E = np.loadtxt('P1Molecule.txt', skiprows=3, unpack=True)

# get distance R between atom positions
Rdist = [abs(pos1[i] - pos2[i]) for i in range(len(E))]

    
# function for Morse potential to model molecule interatomic potential
#*******************************************************************
f = lambda a,b,c,x: (a - a*exp(b*(c-x)))**2 - a
dfa = lambda a,b,c,x: 2*a*(exp(b*(c-x)))**2 - 1
dfb = lambda a,b,c,x: 2*(a**2)*(c-x)*exp(b*(c-x))*(exp(b*(c-x)) - 1)
dfc = lambda a,b,c,x: 2*(a**2)*b*exp(b*(c-x))*(exp(b*(c-x)) - 1)
#*******************************************************************

# call gauss newton method from my least squares module
x,Morse = gaussNewton3(Rdist,E,f,dfa,dfb,dfc,2.6,1.1,1.7,2)


# plot data set with fractal dimension fit model
plt.xlabel('Energy (eV)')
plt.ylabel('Atom Separation (angstrom)')
plt.plot(x,E,'bo', linewidth=1.8, label='Data Set')

plt.plot(x,Morse, 'r', linewidth=1.8, label='Morse Potential')
plt.annotate('$variance=0.0015$',fontsize=18,xy=(0.63,0.50),xycoords='figure fraction')
plt.annotate('$D_e=2.6138$',fontsize=18,xy=(0.67,0.70),xycoords='figure fraction')
plt.annotate(r'$\beta=1.2578$',fontsize=18,xy=(0.68,0.65),xycoords='figure fraction')
plt.annotate('$R_0=1.6453$',fontsize=18,xy=(0.67,0.60),xycoords='figure fraction')
plt.legend()
