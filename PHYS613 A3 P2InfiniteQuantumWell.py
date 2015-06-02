"""
Created on Mon Sep 16 18:31:35 2013
PHYS 613, Assignment 3
Nick Crump
"""

# Problem 2 (EBB2)
"""
Generating the density of states (DOS) by finding energy eigenvalues of 
a particle in a 2D infinite square well. 
"""

import numpy as np
from math import exp
from LeastSquares import gaussNewton3
import matplotlib.pyplot as plt

# Part A
# generate the lowest 500,000 energy eigenvalues from E(nx,ny) = nx**2 + ny**2
# ----------------------------------------------------------------------------
nstates = 500000  # keep this many states
nxy = 1000        # iterate over this many to make sure to get 500,000 lowest
states = []       # stores all E(nx,ny)

# iterate to get all E(nx,ny)
for i in range(1,nxy):
    for j in range(1,nxy):
        s = i**2 + j**2
        states.append(s)

# sort E(nx,ny) list and keep only lowest 500,000
states = np.sort(states)
states = states[0:nstates]    

# plot histogram of density of states (DOS) vs energy
plt.figure(1)
plt.hist(states,bins=20,label='2D Infinite Well')
plt.xlabel('Energy')
plt.ylabel('Density of States')
plt.legend()


# Part B & C
# find energy differences E(i+1)-E(i) and E(i+2)-E(i) up to i = 100000
# ----------------------------------------------------------------------------
# keep duplicates (degeneracies) from array of energy states
num = 100000

# get E(i+1) - E(i)
delta1 = [states[i+1] - states[i] for i in range(num)]
# get E(i+2) - E(i)
delta2 = [states[i+2] - states[i] for i in range(num)]

# plot histogram of energy differences vs energy
plt.figure(2)
plt.subplot(211)
# gets the histogram data and plots histogram
n1, bins1, patches1 = plt.hist(delta1,bins=20,label='First Neighboring Levels')
plt.xlim(0,20)
plt.ylabel(r'Frequency of $\Delta E_{1}$')
plt.legend()

plt.subplot(212)
# gets the histogram data and plots histogram
n2, bins2, patches2 = plt.hist(delta2,bins=20,label='Second Neighboring Levels')
plt.xlim(0,20)
plt.xlabel(r'$\Delta\ $ Energy')
plt.ylabel(r'Frequency of $\Delta E_{2}$')
plt.legend()


# Part D
# fitting parameters of a known function to distributions
# ----------------------------------------------------------------------------
# guess at parameters for exponential decay f = A*exp(-Bx) to fit dataset 1
A1 = 75600
B1 = 1.1
x1 = np.linspace(min(bins1),max(bins1),1000)
exp1 = [A1*exp(-B1*i) for i in x1]
    
    
# guess at parameters for exponential decay f = A*exp(-Bx) to fit dataset 1
A2 = 52100
B2 = 0.7
x2 = np.linspace(min(bins2),max(bins2),1000)
exp2 = [A2*exp(-B2*i) for i in x2]
    
    
# plot histogram of energy differences vs energy with exponential fits
plt.figure(3)
plt.subplot(211)
plt.hist(delta1,bins=20,label='First Neighboring Levels')
plt.plot(x1+0.7,exp1,'r-',linewidth=3,label='Exponential Fit')
plt.xlim(0,20)
plt.ylabel(r'Frequency of $\Delta E_{1}$')
plt.legend()

plt.subplot(212)
plt.hist(delta2,bins=20,label='Second Neighboring Levels')
plt.plot(x2+0.7,exp2,'r-',linewidth=3,label='Exponential Fit')
plt.xlim(0,20)
plt.xlabel(r'$\Delta\ $ Energy')
plt.ylabel(r'Frequency of $\Delta E_{2}$')
plt.legend()
