"""
Created on Fri Nov 29 21:27:07 2013
PHYS 613, Assignment 11
Nick Crump
"""

# Exercise 7.1
"""
Calculate the traveling wave motion of a string by
solving the 1-D wave equation as a difference equation. 
"""

import numpy as np
import matplotlib.pyplot as plt


# set string parameters
# --------------------------------------------
L = 1.0          # length (meters)
T = 10.0         # tension (newtons)
m = 0.001        # mass (kg)
mu = m/L         # mass density (kg/m)
c = (T/mu)**0.5  # wave velocity (m/s)
# --------------------------------------------

# set iteration parameters
# --------------------------------------------
dt = 0.0001      # seconds
dx = 0.01        # meters
Nt = 20          # time steps
Nx = 100         # space steps
# --------------------------------------------

# check for stability condition
# --------------------------------------------
# sqrt(epsilon) needs to be <= 1 for stability
eps = (c*dt/dx)**2

if eps**0.5 > 1.0:
    print '\n'
    print '******************************'
    print '  WARNING: UNSTABLE SOLUTION'
    print '******************************'
# --------------------------------------------

# step through wave equation solutions
# --------------------------------------------
# initial wave shape U(x,0) and derivative dU(x,0) at t=0
x = np.arange(0,1+dx,dx)
U0 = np.exp(-100*(x-0.5)**2)  # initial wave shape
#dU0 = np.zeros(Nx+1)         # no initial motion
dU0 = -200*c*(x-0.5)*U0       # initially a traveling wave

# calculate wave U(x,1) at first time step t=1
U1 = np.zeros(Nx+1)
for i in range(1,Nx):
    U1[i] = (0.5*eps*(U0[i+1] + U0[i-1]) +(1 - eps)*U0[i] + dt*dU0[i])

# now calculate wave U(x,t) over all time steps from U(x,0) and U(x,1)
U = [U0,U1]
# iterate over time-steps
for i in range(Nt-2):
    Unew = np.zeros(Nx+1)
    
    # iterate over x-steps
    for j in range(1,Nx):
        Unew[j] = eps*(U1[j+1] + U1[j-1]) + 2*(1 - eps)*U1[j] - U0[j]

   # store and rename arrays for next iteration
    U.append(Unew)
    U0 = U1
    U1 = Unew
# --------------------------------------------
    

# plot results
# --------------------------------------------
for i in range(0,20,4):
    plt.plot(x,U[i], label='t = '+str(i))
    plt.xlabel('x')
    plt.ylabel('Amplitude')
    plt.ylim(0,1)
    plt.legend()
# --------------------------------------------