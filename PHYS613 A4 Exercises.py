"""
Created on Sat Sep 21 20:58:16 2013
PHYS 613, Assignment 4
Nick Crump
"""

# Exercise 3.6
# Exercise 3.10
# Exercise 3.20
# Exercise 3.21 & 3.22
"""
From Computational Physics by Devries
"""

import Interpolation as Interp
import FirstDerivative as FD
import LeastSquares as LS
import numpy as np
from math import exp
from scipy.special import jn
import matplotlib.pyplot as plt


# Exercise 3.6
#*******************************************************************
# generate points of Bessel function J(1,x)
x = np.arange(0,12,0.5)
npts = len(x)
bess = [jn(1,i) for i in x]

# call cubic spline method from my interpolation module
xx,yy = Interp.cubicSpline(x,bess,npts,0.1)

# plot data and interpolation
plt.figure(1)
plt.plot([0,12],[0,0],'k--')
plt.plot(3.83,0,'gx',markersize=20)
plt.plot(x,bess,'bo',label='Bessel Points')
plt.plot(xx,yy,'r',label='Cubic Spline Interpolation')
plt.xlabel('x')
plt.ylabel('y')
plt.annotate('$x=3.83$',fontsize=18,xy=(0.38,0.46),xycoords='figure fraction')
plt.legend()
#*******************************************************************


# Exercise 3.10
#*******************************************************************
# define function and actual derivative
func = lambda x: x*exp(x)
dfunc = lambda x: (x+1)*exp(x)

# define range of h and x value to evaluate derivative
h = np.arange(0.05,0.5,0.05)
x = 2

# arrays to store derivatives
bd = []
fd = []
cd = []
actual = []

# get first derivatives by backward, forward, and central difference
for i in h:
    bd.append(FD.backwardDiff(func,x,i))
    fd.append(FD.forwardDiff(func,x,i))
    cd.append(FD.centralDiff(func,x,i))
    actual.append(dfunc(x))
    
# make plots of derivative approximation vs step size h
plt.figure(2)
plt.subplot(311)
plt.plot(h,bd,'bo-',label='Backward')
plt.plot(h,actual,'ro-')
plt.ylim(16,24)
plt.legend(loc=3)

plt.subplot(312)
plt.plot(h,fd,'bo-',label='Forward')
plt.plot(h,actual,'ro-',label='Actual')
plt.ylabel('First Derivative Approximations')
plt.ylim(20,32)
plt.legend(loc=2)

plt.subplot(313)
plt.plot(h,cd,'bo-',label='Central')
plt.plot(h,actual,'ro-')
plt.xlabel('Step Size h')
plt.yticks([20,21,22,23,24])
plt.legend(loc=3)
#*******************************************************************


# Exercise 3.20
#*******************************************************************
# position vs time data from pg136 of Devries
t = np.arange(0,1.1,0.1)
y = [1.67203,1.79792,2.37791,2.66408,2.11245,2.43969,1.88843,1.59447,1.79634,1.07810,0.21066]

# call polynomial fit method from my least squares module
xx,yy = LS.polyFit(t,y,2)

# plot data and least squares quadratic fit
plt.figure(3)
plt.plot(t,y,'bo',label='Data Set')
plt.plot(xx,yy,'r',label='Quadratic Fit')
plt.xlabel('Time')
plt.ylabel('Height')
plt.annotate(r'$y(t)=1.654+3.903t-\frac{1}{2}10.404t^{2}$',fontsize=18,xy=(0.20,0.40),xycoords='figure fraction')
plt.legend()
#*******************************************************************


# Exercise 3.21 & 3.22
#*******************************************************************
# intensity vs wavelength data from pg144 of Devries
x = np.arange(435.784,435.884,0.005)
y = [40.0,44.0,41.0,46.0,47.0,54.0,66.0,97.0,129.0,153.0,165.0,168.0,143.0,111.0,79.0,64.0,52.0,51.0,44.0,46.0,41.0]


# function for Lorentzian lineshape to model light intensities vs wavelength
#*******************************************************************
f = lambda a,b,c,d,x: d + a/(1.0 + 4.0*(x - b)**2/c**2)
dfa = lambda a,b,c,d,x: 1.0/(1.0 + 4.0*(x - b)**2/c**2)
dfb = lambda a,b,c,d,x: -4.0*a*(-2*x + 2*b)/(c**2*(1.0 + 4.0*(x - b)**2/c**2)**2)
dfc = lambda a,b,c,d,x: 8.0*a*(x - b)**2/(c**3*(1.0 + 4.0*(x - b)**2/c**2)**2)
dfd = lambda a,b,c,d,x: 1.0
#*******************************************************************

# call gauss newton method from my least squares module
xx,yy = LS.gaussNewton4(x,y,f,dfa,dfb,dfc,dfd,165.0,435.84,0.03,40.0,2)

# plot data and least squares quadratic fit
plt.figure(4)
plt.plot(x,y,'bo',label='Data Set')
plt.plot(xx,yy,'r',label='Lorentzian Fit')
plt.xlabel('$\lambda\ (nm)$')
plt.ylabel('Intensity')
plt.xticks(np.arange(435.78,435.90,0.02),['435.78','435.80','435.82','435.84','435.86','435.88','435.90'])
plt.yticks(range(0,250,50))
plt.legend()
#*******************************************************************