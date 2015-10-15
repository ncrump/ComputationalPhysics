"""
Created on Fri Sep 06 21:03:27 2013
PHYS 613, Assignment 2
Nick Crump
"""

# Exercise 2.15
"""
From Computational Physics by Devries
"""

from math import sin,cos
import numpy as np
import matplotlib.pyplot as plt


# define function - even state solutions of 1D finite square well potential
#***********************************************************************  
def evenFunc(a,m,V0,E):
    hbarSq = 0.076199682    # eV(nm**2)(Melectron)
    
    alpha = (2.0*m*E/hbarSq)**0.5
    beta = ((2.0*m*(V0-E))/hbarSq)**0.5
    
    fEven = beta*cos(alpha*a) - alpha*sin(alpha*a)
    return fEven
#***********************************************************************  


# define function - odd state solutions of 1D finite square well potential
#***********************************************************************  
def oddFunc(a,m,V0,E):
    hbarSq = 0.076199682    # eV(nm**2)(Melectron)
    
    alpha = (2.0*m*E/hbarSq)**0.5
    beta = ((2.0*m*(V0-E))/hbarSq)**0.5
    
    fOdd = alpha*cos(alpha*a) + beta*sin(alpha*a)
    return fOdd
#***********************************************************************  


# enter root finding algorithm by Bisection method
#*********************************************************************** 
def rootBisection(f, xI, xF, Tol, nMax):
        # initialize variables    
        error = 1
        n = 1
        xiMid = 0  # initial midpoint value to store the n-1 value
        
        # loop until error is less than input tolerance
        while error > Tol:
            xMid = 0.5*(xI+xF)
            
            # set up main Bisection method:
            # make bracket interval smaller each iteration until root is found
            # check conditions and update bracket points     
            if f(xI)*f(xMid) > 0:
                xI = xMid
                error = abs(xMid - xiMid)  # calculate approx error
                n = n + 1
                xiMid = xMid               # store the n-1 midpoint
                
            elif f(xI)*f(xMid) < 0:
                xF = xMid
                error = abs(xMid - xiMid)  # calculate approx error
                n = n + 1
                xiMid = xMid               # store the n-1 midpoint
        
        # output results to user         
        print 'Root =', xMid
        print 'Iterations = ', n-1
            
# end rootBisection function
#***********************************************************************


# main program that calls functions, finds roots and does plotting
#***********************************************************************

# setup plotting of allowed energy equation as function of energy
#--------------------------
m = 1.0                 # Melectron
a = 0.50                # nm
V0 = 10.0               # eV
hbarSq = 0.076199682    # eV(nm**2)(Melectron)

E = np.arange(0,V0,0.1)    
evenF = []
oddF = []

for i in E:
    fEven = evenFunc(a,m,V0,i)
    fOdd = oddFunc(a,m,V0,i)
    evenF.append(fEven)
    oddF.append(fOdd) 
    
plt.figure(1)
plt.plot([0,V0],[0,0],'k',linewidth='2')
plt.plot(E,evenF,'b',label='Even States')
plt.plot(E,oddF,'r',label='Odd States')
plt.xlabel('Energy (eV)')
plt.ylabel('$f\ (E)$')
plt.legend(loc=9)


# setup root finder routine
#--------------------------
sfEven = lambda E: ((2.0*m*(V0-E)/hbarSq)**0.5)*cos(((2.0*m*E/hbarSq)**0.5)*a) - ((2.0*m*E/hbarSq)**0.5)*sin(((2.0*m*E/hbarSq)**0.5)*a)
sfOdd = lambda E: ((2.0*m*E/hbarSq)**0.5)*cos(((2.0*m*E/hbarSq)**0.5)*a) + ((2.0*m*(V0-E)/hbarSq)**0.5)*sin(((2.0*m*E/hbarSq)**0.5)*a)

rootBisection(sfEven, 0, 0.5, 10e-8, 30)
rootBisection(sfOdd, 0.8, 1.2, 10e-8, 30)


# setup plotting of ground state energies as function of V0 and a
# these found manually by adjusting V0 and a and finding roots
#--------------------------
even1 = [0.47646,0.54218,0.58707,0.62058,0.64695,0.66849,0.68655,0.70201,0.71545]
odd1 = [1.67230,2.01417,2.23445,2.39306,2.51496,2.61278,2.69373,2.76228,2.82139]
V = range(2,11)

even2 = [3.26434,2.05163,1.35689,0.96104,0.71545,0.55299,0.44008,0.35847,0.29760]
odd2 = [9.89048,7.46026,5.19897,3.75528,2.82139,2.19143,1.74897,1.42720,1.18627]
a = [0.104,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50]

plt.figure(2)
fig = plt.subplot(211)
plt.plot(V,even1,'bo-',label='Even States')
plt.plot(V,odd1,'ro-',label='Odd States')
plt.annotate('$a=0.3\ nm$',fontsize=14,xy=(0.67,0.80),xycoords='figure fraction')
plt.xlabel('Potential $V_0 (eV)$')
plt.ylabel('Ground State Energies')
plt.legend(loc=7)

fig = plt.subplot(212)
plt.plot(a,even2,'bo-',label='Even States')
plt.plot(a,odd2,'ro-',label='Odd States')
plt.annotate('$V_0=10\ eV$',fontsize=14,xy=(0.67,0.34),xycoords='figure fraction')
plt.xlabel('Distance $a (nm)$')
plt.ylabel('Ground State Energies')
plt.legend(loc=7)
#***********************************************************************
