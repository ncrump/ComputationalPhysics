"""
Created on Mon Sep 16 18:31:35 2013
PHYS 613, Assignment 3
Nick Crump
"""

# Problem 1 (EBB1)
"""
Finding bound state eigenvalues of a wavefunction in 1D square well. 
Well has one infinite wall and one finite wall. 
"""

from math import sin,cos,exp,sqrt
import sympy
import numpy as np
import matplotlib.pyplot as plt


# define function - equation of allowed energies of single wall finite well
#***********************************************************************  
def wavFunc(L,m,V0,E):
    hbarSq = 0.076199682    # eV(nm**2)(Melectron)
    
    psi = sqrt(V0-E)*sin(L*sqrt(2.0*m*E/hbarSq)) + sqrt(E)*cos(L*sqrt(2.0*m*E/hbarSq))
    return psi
#***********************************************************************  


# enter root finding algorithm by hybrid Bisection/Newton-Raphson method
#*********************************************************************** 
def rootNewtonRaphson(g, f, xI, xF, Tol, nMax):
        # initialize variables    
        error = 1
        n = 1
        
        # get symbolic derivative of input function      
        x = sympy.Symbol('x')        # define x as symbolic variable
        sdf = sympy.diff(g, x)       # get symbolic derivative of input function g
        df = sympy.lambdify(x, sdf)  # turn symbolic derivative into numeric function
        
        # check condition for starting x value
        if f(xI) < f(xF): xn0 = xI
        else: xn0 = xF
        
        # loop until error is less than input tolerance
        while error > Tol:
            
            # set up main Newton-Raphson method:
            # use derivative of function as tangent line to get new x point 
            # calcuate new x value from f(x) and df(x)
            xn1 = xn0 - (f(xn0)/df(xn0))
            
            # if new x point is outside bracket interval then use midpoint instead
            if xn1 < xI or xn1 > xF:
                xn1 = 0.5*(xI+xF)
                
                # check conditions and update bracket points     
                if f(xI)*f(xn1) > 0:
                    xI = xn1
                    
                elif f(xI)*f(xn1) < 0:
                    xF = xn1
                
            error = abs(xn1 - xn0)  # calculate approx error
            n = n + 1
            xn0 = xn1               # store the n-1 value for x
            
        # output results to user
        return round(xn1,5)
            
# end rootNewtonRaphson function
#***********************************************************************


# main program that calls functions, finds roots and does plotting
#***********************************************************************

# setup plotting of allowed energy equation as function of energy
#--------------------------
L = 2.00000             # nm  (20 angsroms)
m = 1.09770             # Melectron  (10**-27 grams)
V0 = 0.93623            # eV  (15x10**-13 ergs)
hbarSq = 0.076199682    # eV(nm**2)(Melectron)

E = np.arange(0,V0,0.01) 
num = len(E)   
psiArrE = np.zeros(num)

for i in range(num):
    psi = wavFunc(L,m,V0,E[i])
    psiArrE[i] = psi 
    
plt.figure(1)
plt.plot([0,1],[0,0],'k',linestyle='--')
plt.plot([0.07110,0.28213,0.62268],[0,0,0],'bo')
plt.plot(E,psiArrE,'b',label='Allowed Energy Function')
plt.xlabel('Energy (eV)')
plt.ylabel('$f\ (E)$')
plt.legend()


# setup root finder routine
#--------------------------
g = 'sqrt(0.93623-x)*sin(2.0*sqrt(2.0*1.09770*x/0.076199682)) + sqrt(x)*cos(2.0*sqrt(2.0*1.09770*x/0.076199682))'
f = lambda x: eval(g)
E1 = rootNewtonRaphson(g, f, 0, 0.2, 10e-5, 30)
E2 = rootNewtonRaphson(g, f, 0.2, 0.4, 10e-5, 30)
E3 = rootNewtonRaphson(g, f, 0.6, 0.8, 10e-5, 30)
print 'Eigenvalues = ',E1,E2,E3


# setup wavefunction plotting as function of distance & plot potential well
#--------------------------
# x arrays for regions around well
R1 = np.arange(0,2.1,0.1)     # region 1 inside well 
R2 = np.arange(2.0,6.0,0.1)  # region 2 right of well

# value of alpha as function of En
alphE1 = sqrt(2.0*m*E1/hbarSq)
alphE2 = sqrt(2.0*m*E2/hbarSq)
alphE3 = sqrt(2.0*m*E3/hbarSq)

# value of beta as function of En
betaE1 = sqrt(2*m*(V0-E1)/hbarSq)
betaE2 = sqrt(2*m*(V0-E2)/hbarSq)
betaE3 = sqrt(2*m*(V0-E3)/hbarSq)

# wavefunctions for the 3 energy levels in both regions (arbitrary normalization coefficients)
# wavefunctions shifted to make energy eigenvalues the zero baseline
#--------------------------
# wavefunctions for E1
psiR1E1 = [0.1*sin(alphE1*i)+E1 for i in R1]
psiR2E1 = [500*exp(-betaE1*i)+E1 for i in R2]
# wavefunctions for E2
psiR1E2 = [0.1*sin(alphE2*i)+E2 for i in R1]
psiR2E2 = [-350*exp(-betaE2*i)+E2 for i in R2]
# wavefunctions for E3
psiR1E3 = [0.1*sin(alphE3*i)+E3 for i in R1]
psiR2E3 = [30*exp(-betaE3*i)+E3 for i in R2]

# probability density functions |Psi|**2
#--------------------------
# probability function for E1
psiR1E1Sq = [0.1*(sin(alphE1*i)**2)+E1 for i in R1]
psiR2E1Sq = [5000*(exp(-betaE1*i)**2)+E1 for i in R2]
# probability function for E2
psiR1E2Sq = [0.1*(sin(alphE2*i)**2)+E2 for i in R1]
psiR2E2Sq = [1000000*(exp(-betaE2*i)**2)+E2 for i in R2]
# probability function for E3
psiR1E3Sq = [0.1*(sin(alphE3*i)**2)+E3 for i in R1]
psiR2E3Sq = [10000*(exp(-betaE3*i)**2)+E3 for i in R2]


# plot wavefunctions
#--------------------------
plt.figure(2)
# plot potential V0
plt.plot([0,0],[0,1.4],'k',linewidth='4')
plt.plot([0,L],[0,0],'k',linewidth='4')
plt.plot([L,L],[0,V0], 'k',linewidth='4')
plt.plot([L,10], [V0,V0], 'k',linewidth='4')
# plot lines for energies E1, E2, E3
plt.plot([-2,6],[E1,E1],'g',linewidth='2',linestyle='--')
plt.plot([-2,6],[E2,E2],'g',linewidth='2',linestyle='--')
plt.plot([-2,6],[E3,E3],'g',linewidth='2',linestyle='--')
# plot wavefunctions
# wavefunction psiE1
plt.plot(R1,psiR1E1,'b',label='$\psi_{E1}\ (x)$')
plt.plot(R2,psiR2E1,'b')
# wavefunction psiE2
plt.plot(R1,psiR1E2,'r',label='$\psi_{E2}\ (x)$')
plt.plot(R2,psiR2E2,'r')
# wavefunction psiE3
plt.plot(R1,psiR1E3,'g',label='$\psi_{E3}\ (x)$')
plt.plot(R2,psiR2E3,'g')
# set x,y plot limits
plt.xlim(-2.0,6.0)
plt.ylim(-0.1,1.4)
# potential annotations
plt.annotate(r'$V=\infty$',fontsize=16,xy=(0.18,0.85),xycoords='figure fraction')
plt.annotate(r'$V=0$',fontsize=16,xy=(0.38,0.85),xycoords='figure fraction')
plt.annotate(r'$V=V_0$',fontsize=16,xy=(0.55,0.85),xycoords='figure fraction')
# energy annotations
plt.annotate('$E_{1}=0.07110$',fontsize=12,xy=(0.14,0.205),xycoords='figure fraction')
plt.annotate('$E_{2}=0.28213$',fontsize=12,xy=(0.14,0.32),xycoords='figure fraction')
plt.annotate('$E_{3}=0.62268$',fontsize=12,xy=(0.14,0.50),xycoords='figure fraction')
# label axes
plt.xlabel('Distance (nm)')
plt.ylabel('Wavefunctions')
plt.legend()


# plot probability density functions
#--------------------------
plt.figure(3)
# plot potential V0
plt.plot([0,0],[0,1.4],'k',linewidth='4')
plt.plot([0,L],[0,0],'k',linewidth='4')
plt.plot([L,L],[0,V0], 'k',linewidth='4')
plt.plot([L,10], [V0,V0], 'k',linewidth='4')
# plot lines for energies E1, E2, E3
plt.plot([-2,6],[E1,E1],'g',linewidth='2',linestyle='--')
plt.plot([-2,6],[E2,E2],'g',linewidth='2',linestyle='--')
plt.plot([-2,6],[E3,E3],'g',linewidth='2',linestyle='--')
# plot wavefunctions
# wavefunction psiE1
plt.plot(R1,psiR1E1Sq,'b',label='$|\psi_{E1}|^2\ (x)$')
plt.plot(R2,psiR2E1Sq,'b')
# wavefunction psiE2
plt.plot(R1,psiR1E2Sq,'r',label='$|\psi_{E2}|^2\ (x)$')
plt.plot(R2,psiR2E2Sq,'r')
## wavefunction psiE3
plt.plot(R1,psiR1E3Sq,'g',label='$|\psi_{E3}|^2\ (x)$')
plt.plot(R2,psiR2E3Sq,'g')
# set x,y plot limits
plt.xlim(-2.0,6.0)
plt.ylim(-0.1,1.4)
# potential annotations
plt.annotate(r'$V=\infty$',fontsize=16,xy=(0.18,0.85),xycoords='figure fraction')
plt.annotate(r'$V=0$',fontsize=16,xy=(0.38,0.85),xycoords='figure fraction')
plt.annotate(r'$V=V_0$',fontsize=16,xy=(0.55,0.85),xycoords='figure fraction')
# energy annotations
plt.annotate('$E_{1}=0.07110$',fontsize=12,xy=(0.14,0.205),xycoords='figure fraction')
plt.annotate('$E_{2}=0.28213$',fontsize=12,xy=(0.14,0.32),xycoords='figure fraction')
plt.annotate('$E_{3}=0.62268$',fontsize=12,xy=(0.14,0.50),xycoords='figure fraction')
# label axes
plt.xlabel('Distance (nm)')
plt.ylabel('Probability Density Functions')
plt.legend()
#***********************************************************************