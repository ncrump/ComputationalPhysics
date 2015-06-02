"""
Created on Sat Nov 09 13:53:23 2013
PHYS 613, Assignment 8 
Nick Crump
"""

# Exercise 5.16
"""
Implement the Runge-Kutta-Fehlberg method to solve the second 
order ODE below and compare to the analytic result. 
"""

import numpy as np
import ODEsolve as ode
import matplotlib.pyplot as plt

# function for ODE y''(t) = -4y(t), y(0)=1, y'(0)=0
# ***************************************************
def f(ti,IC): 
    ui = IC[0]  # y'(t)
    yi = IC[1]  # y (t)
    
    func = np.array([-4*yi, ui])
    return func
# ***************************************************

# set initial conditions as y'(0), y(0)
IC = [0, 1]
# call RKF45 method from my ODE solve module
t,s = ode.RKF45HO(f,IC,0,2*np.pi,0.01,1e-5,stp=1)

# get analytic solution for comparison
# this ODE is the harmonic oscillator for w = 2
yAct = IC[1]*np.cos(2*t) + 0.5*np.sin(2*t)  

# plot numerical and analytic solutions
plt.plot(t,s[1],'b', label='RKF45 Solution')
plt.plot(t,yAct,'r', label='Analytic Solution')
plt.xlabel('Time')
plt.ylabel('Position')
plt.legend()