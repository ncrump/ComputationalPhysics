"""
Created on Wed Sep 25 23:14:56 2013
PHYS 613, Assignment 5
Nick Crump
"""

# Exercise 4.4
# Exercise 4.6
"""
From Computational Physics by Devries
"""

import Integration as Int
import numpy as np
from math import sqrt,sin,pi
import matplotlib.pyplot as plt


# Exercise 4.4
#*******************************************************************
v = np.linspace(-4,4,256)  # angles for Fresnel equations

# arrays to store Fresnel Eq values and intensities
Cv = []
Sv = []
Iv = []

# generate Fresnel Eqns and intensities
for i in v:
    # call Romberg method from my integration module
    Cvi = Int.romberg('cos((pi/2)*x**2)',0,i,10)
    Svi = Int.romberg('sin((pi/2)*x**2)',0,i,10)
    Ivi = 0.5*((Cvi+0.5)**2 + (Svi+0.5)**2)
    
    # append values
    Cv.append(Cvi)
    Sv.append(Svi)
    Iv.append(Ivi)


# plot Intensity
plt.figure(1)
plt.subplot(211)
plt.plot(v,Iv,'g',label='Intensity')
plt.ylabel('Intensity')
plt.legend(loc=4)

# plot Fresnel Eqns C(v) and S(v)
plt.subplot(212)
plt.plot(v,Cv,'b',label='Fresnel $C(v)$')
plt.plot(v,Sv,'r',label='Fresnel $S(v)$')
plt.xlabel('Distance $v$')
plt.ylabel('Fresnel Equations')
plt.legend(loc=4)

# plot C(v) vs S(v) to get Cornu Spiral
plt.figure(2)
plt.plot(Cv,Sv,label='Cornu Spiral')
plt.plot([-0.5,0.5],[-0.5,0.5],'rx',markersize=8)
plt.xlabel('Fresnel $C\ (v)$')
plt.ylabel('Fresnel $S\ (v)$')
plt.legend(loc=2)
#*******************************************************************


# Exercise 4.6
#*******************************************************************
# set initial parameters
L = 1.0  # length of string (m)
g = 9.8  # gravity (m/s**2)

# period for simple pendulum (using small angle approx)
Tsimp = (2*pi)*sqrt(L/g)  # (sec)

# arrays 
period = []                  # calculated periods (sec)
theta0 = np.arange(0,180,10)  # initial angles (deg)

# iterate over angle theta0 (initial angle at t=0)
for i in theta0:
    # argument of elliptical integral (converted to radians)
    k = sin((i/2)*(pi/180))
    
    # integrand
    f = '1.0/(1-('+str(k)+'**2)*sin(x)**2)**0.5'
    
    # call Romberg method from my integration module
    romb = Int.romberg(f,0,pi/2,10)
    
    # compute and store periods for each initial angle
    p = 4.0*sqrt(L/g)*romb
    period.append(p)
    
    # check when period differs from simple pendulum period by more than 1%
    if abs(p-Tsimp) > 0.01*Tsimp: print i,p,Tsimp
        
  
# plot results   
plt.figure(3)
plt.plot(theta0,period,'bo-')
plt.xlabel('Amplitude (deg)')
plt.ylabel('Period (sec)')
#*******************************************************************