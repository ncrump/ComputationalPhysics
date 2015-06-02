"""
Created on Sat Nov 09 13:53:23 2013
PHYS 613, Assignment 8 
Nick Crump
"""

# Exercise 5.9: Projectile Motion
"""
Implement the Runge-Kutta-Fehlberg method to solve the 1st-order
vector differential equation of projectile motion. 
"""

import numpy as np
import ODEsolve as ode
import matplotlib.pyplot as plt


# function for projectile motion w/air resistance
# ***************************************************
def f(ti,IC): 
    g = 9.8
    c = 0.1
    vx = IC[0]
    vy = IC[1]  
    
    vxy = (vx**2 + vy**2)**0.5
    func = np.array([-c*vx*vxy, -g-c*vy*vxy])
    return func
# ***************************************************

# set initial conditions as Vx0, Vy0
IC = [20, 0]
# call RKF45 method from my ODE solve module
t,s = ode.RKF45HO(f,IC,0,10,0.01,1e-5,stp=1)

# plot components of velocity
plt.plot(t,s[0],'b.-',label='$V_x$')
plt.plot(t,s[1],'r.-',label='$V_y$')
plt.xlabel('Time (sec)')
plt.ylabel('Velocity (m/s)')
plt.legend()