"""
Created on Sat Nov 09 13:53:23 2013
PHYS 613, Assignment 8 
Nick Crump
"""

# Exercise 5.2
"""
Implement the Euler methods to solve the first order
ODE below and compare results for different step sizes. 
"""

import ODEsolve as ode
import matplotlib.pyplot as plt


# function for ODE y'(x) = y**2 + 1, y(0) = 0
# ***************************************************
def f(x,y): 
    func = y**2 + 1
    return func
# ***************************************************

# define step sizes
h = [0.20, 0.10, 0.05]
n = len(h)

# loop through step sizes to plot solutions
for i in range(n):
    # call Euler methods from my ODE solve module
    x1,y1 = ode.simpleEuler(f,0,0,1,h[i])
    x2,y2 = ode.modifiedEuler(f,0,0,1,h[i])
    x3,y3 = ode.improvedEuler(f,0,0,1,h[i])
    
    # plot solutions for each step size
    plt.subplot(3,1,i+1)
    plt.plot(x1,y1,'b.-',linewidth=1.5,label='h='+str(h[i]))
    plt.plot(x2,y2,'r.-',linewidth=1.5)
    plt.plot(x3,y3,'g.-',linewidth=1.5)
    plt.legend(loc=2)