"""
Created on Fri Sep 06 21:03:27 2013
PHYS 613, Assignment 2
Nick Crump

"""
# Exercise 2.13
# Exercise 2.14
"""
From Computational Physics by Devries
"""

from math import sin,cos,exp,sqrt
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
        return round(xMid,5)
            
# end rootBisection function
#***********************************************************************


# main program that calls functions, finds roots and does plotting
#***********************************************************************

# setup root finder routine
#--------------------------
a = 0.3                 # nm
m = 1.0                 # Melectron
V0 = 10.0               # eV
hbarSq = 0.076199682    # eV(nm**2)(Melectron)

sfEven = lambda E: ((2.0*m*(V0-E)/hbarSq)**0.5)*cos(((2.0*m*E/hbarSq)**0.5)*a) - ((2.0*m*E/hbarSq)**0.5)*sin(((2.0*m*E/hbarSq)**0.5)*a)
sfOdd = lambda E: ((2.0*m*E/hbarSq)**0.5)*cos(((2.0*m*E/hbarSq)**0.5)*a) + ((2.0*m*(V0-E)/hbarSq)**0.5)*sin(((2.0*m*E/hbarSq)**0.5)*a)

Eeven = rootBisection(sfEven, 0, 2.0, 10e-5, 30)
Eodd = rootBisection(sfOdd, 2.0, 4.0, 10e-5, 30)
print 'Eigenvalues = ', Eeven, Eodd


# setup plotting of allowed energy equation as function of energy
#--------------------------
E = np.arange(0,10.1,0.1)    
evenF = []
oddF = []

for i in E:
    fEven = evenFunc(0.3,1.0,10.0,i)
    fOdd = oddFunc(0.3,1.0,10.0,i)
    evenF.append(fEven)
    oddF.append(fOdd) 
    
plt.figure(1)
plt.plot(E,evenF,'b',label='Even States')
plt.plot(E,oddF,'r',label='Odd States')
plt.plot(Eeven,0,'bo',Eodd,0,'ro')
plt.xlabel('Energy (eV)')
plt.ylabel('$f\ (E)$')
plt.legend(loc=9)


# setup wavefunction plotting as function of distance & plot potential well
#--------------------------
# x arrays for regions around well
R1 = np.arange(-0.6,-0.29,0.01)  # region 1 left of well 
R2 = np.arange(-0.3,0.301,0.01)  # region 2 inside well
R3 = np.arange(0.3,0.601,0.01)   # region 3 right of well

# alpha & beta values for even states
alphEven = sqrt(2*m*Eeven/hbarSq)
betaEven = sqrt(2*m*(V0-Eeven)/hbarSq)

# even state wavefunctions for 3 regions (arbitrary normalization coefficients)
# wavefunctions shifted to make energy eigenvalues the zero baseline
psiR1even = [30*exp(betaEven*i)+Eeven for i in R1]
psiR2even = [cos(alphEven*i)+Eeven for i in R2]
psiR3even = [30*exp(-betaEven*i)+Eeven for i in R3]

# alpha & beta values for odd states
alphOdd = sqrt(2*m*Eodd/hbarSq)
betaOdd = sqrt(2*m*(V0-Eodd)/hbarSq) 

# odd state wavefunctions for 3 regions (arbitrary normalization coefficients)
# wavefunctions shifted to make energy eigenvalues the zero baseline
psiR1odd = [-30*exp(betaOdd*i)+Eodd for i in R1]
psiR2odd = [sin(alphOdd*i)+Eodd for i in R2]
psiR3odd = [30*exp(-betaOdd*i)+Eodd for i in R3]

    
plt.figure(2)
# plot lines for potential V(x)
plt.plot([-0.6,-0.3],[10,10],'k',linewidth='4')
plt.plot([-0.3,-0.3],[10,0],'k',linewidth='4')
plt.plot([-0.3,0.3],[0,0], 'k',linewidth='4')
plt.plot([0.3,0.3], [0,10], 'k',linewidth='4')
plt.plot([0.3,0.6],[10,10], 'k',linewidth='4')
plt.xticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6])
plt.annotate('$V_0$',fontsize=16,xy=(0.23,0.82),xycoords='figure fraction')
# plot lines for energy eigenvalues
plt.plot([-0.6,0.6],[Eeven,Eeven],'g',linewidth='2',linestyle='--')
plt.plot([-0.6,0.6],[Eodd,Eodd],'g',linewidth='2',linestyle='--')
plt.annotate('Ground State Energies',fontsize=12,xy=(0.39,0.27),xycoords='figure fraction')
plt.annotate('$E_{even}=0.71545$',fontsize=12,xy=(0.75,0.20),xycoords='figure fraction')
plt.annotate('$E_{odd}=2.82139$',fontsize=12,xy=(0.755,0.40),xycoords='figure fraction')
# plot wavefunctions for each ground state energy
plt.plot(R1,psiR1even,'b',label='$\psi_{even}\ ({x})$')
plt.plot(R2,psiR2even,'b')
plt.plot(R3,psiR3even,'b')
plt.plot(R1,psiR1odd,'r',label='$\psi_{odd}\ ({x})$')
plt.plot(R2,psiR2odd,'r')
plt.plot(R3,psiR3odd,'r')
plt.annotate(r'$\psi_{1}=C\ \exp({\beta x})$',fontsize=12,xy=(0.15,0.625),xycoords='figure fraction')
plt.annotate(r'$\psi_{2odd}=A\ \sin({\alpha x})$',fontsize=12,xy=(0.42,0.65),xycoords='figure fraction')
plt.annotate(r'$\psi_{2even}=B\ \cos({\alpha x})$',fontsize=12,xy=(0.42,0.60),xycoords='figure fraction')
plt.annotate(r'$\psi_{3}=F\ \exp({-\beta x})$',fontsize=12,xy=(0.73,0.625),xycoords='figure fraction')
plt.yticks(range(-2,14,2))
# set titles
plt.xlabel('Distance (nm)')
plt.ylabel('Ground State Wavefunctions')
plt.legend(loc=9)
#***********************************************************************