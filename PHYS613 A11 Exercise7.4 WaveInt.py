"""
Created on Fri Nov 29 21:27:07 2013
PHYS 613, Assignment 11
Nick Crump
"""

# Exercise 7.4
"""
Calculate the interference wave motion of a string by
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
Nt = 300          # time steps
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
x = np.arange(0,1+dx,dx)

# sign to show constructive/destructive interference
a = 1.0
b = 1.0

# traveling waves f(x,0) and g(x,0)
f0 = a*np.exp(-100*(x-0.75)**2)
g0 = b*np.exp(-100*(x-0.25)**2)

# derivatives of traveling waves df(x,0) and dg(x,0)
df0 = -200*c*(x-0.75)*f0
dg0 = 200*c*(x-0.25)*g0

# initial wave shape U(x,0) and derivative dU(x,0) at t=0
U0 = f0 + g0
dU0 = df0 + dg0

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
    
    
# plot results (3x3 subplot the brute force way...)
# --------------------------------------------
plt.subplot(3,3,1)
plt.plot(x,U[0])
plt.xticks([])
plt.yticks([])
plt.annotate('0',fontsize=15,xy=(0.14,0.86),xycoords='figure fraction',color='r')
plt.subplot(3,3,2)
plt.plot(x,U[20])
plt.xticks([])
plt.yticks([])
plt.annotate('20',fontsize=15,xy=(0.415,0.86),xycoords='figure fraction',color='r')
plt.subplot(3,3,3)
plt.plot(x,U[40])
plt.xticks([])
plt.yticks([])
plt.annotate('40',fontsize=15,xy=(0.685,0.86),xycoords='figure fraction',color='r')
plt.subplot(3,3,4)
plt.plot(x,U[60])
plt.xticks([])
plt.yticks([])
plt.annotate('60',fontsize=15,xy=(0.14,0.57),xycoords='figure fraction',color='r')
plt.subplot(3,3,5)
plt.plot(x,U[80])
plt.xticks([])
plt.yticks([])
plt.annotate('80',fontsize=15,xy=(0.415,0.57),xycoords='figure fraction',color='r')
plt.subplot(3,3,6)
plt.plot(x,U[100])
plt.xticks([])
plt.yticks([])
plt.annotate('100',fontsize=15,xy=(0.685,0.57),xycoords='figure fraction',color='r')
plt.subplot(3,3,7)
plt.plot(x,U[120])
plt.xticks([])
plt.yticks([])
plt.annotate('120',fontsize=15,xy=(0.14,0.29),xycoords='figure fraction',color='r')
plt.subplot(3,3,8)
plt.plot(x,U[140])
plt.xticks([])
plt.yticks([])
plt.annotate('140',fontsize=15,xy=(0.415,0.29),xycoords='figure fraction',color='r')
plt.subplot(3,3,9)
plt.plot(x,U[160])
plt.xticks([])
plt.yticks([])
plt.annotate('160',fontsize=15,xy=(0.685,0.29),xycoords='figure fraction',color='r')