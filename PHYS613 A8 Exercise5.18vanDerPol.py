"""
Created on Sat Nov 09 13:53:23 2013
PHYS 613, Assignment 8 
Nick Crump
"""

# Exercise 5.18
# Exercise 5.19
"""
Implement the Runge-Kutta-Fehlberg method to solve the 2nd-order
non-linear differential equation of the van der Pol oscillator.
"""

import numpy as np
import ODEsolve as ode
import matplotlib.pyplot as plt


# function for van der Pol oscillator
# ***************************************************
def f(ti,IC): 
    eps = 1.0
    ui = IC[0]  # x'(t)
    xi = IC[1]  # x (t)
    
    func = np.array([-xi-eps*ui*(xi**2 - 1), ui])
    return func
# ***************************************************

# set initial conditions as x'(0), x(0)
IC = [0, 0.5]
# call RKF45 method from my ODE solve module
t,s = ode.RKF45HO(f,IC,0,8*np.pi,0.01,1e-5,stp=1)

# plot position and velocity vs time
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,s[1],'b',label='Position')
plt.ylabel('Position')
plt.legend(loc=4)
plt.subplot(2,1,2)
plt.plot(t,s[0],'r',label='Velocity')
plt.xlabel('Time')
plt.ylabel('Velocity')
plt.legend(loc=4)

# plot phase space trajectory
plt.figure(2)
plt.plot(s[1],s[0],'g',label='Phase Space')
plt.xlabel('Position')
plt.ylabel('Velocity')
plt.legend(loc=2)