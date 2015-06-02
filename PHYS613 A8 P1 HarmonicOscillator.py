"""
Created on Sat Nov 09 13:53:23 2013
PHYS 613, Assignment 8 
Nick Crump
"""

# Problem 1 (EBB11): 1D Harmonic Oscillator
"""
Solves the 1D harmonic oscillator using simple and modified Euler methods. 
"""

import numpy as np
import matplotlib.pyplot as plt


# solve 1D harmonic oscillator by simple Euler method
# ***************************************************
#----------------------------------------------------
# t0,ft = time to oscillate over
# x0 = initial position
# v0 = initial velocity
# h = time step
# w = oscillation freq (omega)
# m = mass
#----------------------------------------------------
def simpleEulerSHM(t0,tf,x0,v0,h,w,m):
    
    # define initial parameters
    k = m*(w**2)        # spring constant
    pe = 0.5*k*(x0**2)  # potential energy
    ke = 0.5*m*(v0**2)  # kinetic energy
    e = pe + ke         # total energy
    
    # setup array for storing solutions
    t = np.arange(t0,tf,h)  # time
    x = [x0]                # position
    v =[v0]                 # velocity
    PE = [pe]               # kinetic energy
    KE = [ke]               # potential energy
    E = [e]                 # total energy

    pts = len(t)-1
    # loop to solve harmonic oscillator
    for i in range(pts):
        # simple Euler method
        xi = x[i] + h*(w**2)*v[i]
        vi = v[i] - h*x[i]
        
        # get energies
        pe = 0.5*k*(xi**2)
        ke = 0.5*m*(vi**2)
        e = pe + ke

        # store to arrays
        x.append(xi)
        v.append(vi)
        PE.append(pe)
        KE.append(ke)
        E.append(e)
        
    return t,x,v,PE,KE,E
# ***************************************************
    
    
# solve 1D harmonic oscillator by modified Euler method
# ***************************************************
#----------------------------------------------------
# t0,ft = time to oscillate over
# x0 = initial position
# v0 = initial velocity
# h = time step
# w = oscillation freq (omega)
# m = mass
#----------------------------------------------------
def modifiedEulerSHM(t0,tf,x0,v0,h,w,m):
    
    # define initial parameters
    k = m*(w**2)        # spring constant
    pe = 0.5*k*(x0**2)  # potential energy
    ke = 0.5*m*(v0**2)  # kinetic energy
    e = pe + ke         # total energy
    
    # setup array for storing solutions
    t = np.arange(t0,tf,h)  # time
    x = [x0]                # position
    v =[v0]                 # velocity
    PE = [pe]               # kinetic energy
    KE = [ke]               # potential energy
    E = [e]                 # total energy

    pts = len(t)-1
    # loop to solve harmonic oscillator
    for i in range(pts):
        # modified Euler method
        xii = x[i] + 0.5*h*(w**2)*v[i]
        vii = v[i] - 0.5*h*x[i]
        xi = x[i] + h*vii
        vi = v[i] - h*xii
        
        # get energies
        pe = 0.5*k*(xi**2)
        ke = 0.5*m*(vi**2)
        e = pe + ke

        # store to arrays
        x.append(xi)
        v.append(vi)
        PE.append(pe)
        KE.append(ke)
        E.append(e)
        
    return t,x,v,PE,KE,E
# ***************************************************

# set initial parameters
t0 = 0
tf = 10
x0 = 0
v0 = 1
h = 0.1
w = 1
m = 1
# call Euler SHM methods above
t1,x1,v1,PE1,KE1,E1 = simpleEulerSHM(t0,tf,x0,v0,h,w,m) 
t2,x2,v2,PE2,KE2,E2 = modifiedEulerSHM(t0,tf,x0,v0,h,w,m) 

# get analytic solution for comparison
xAct = x0*np.cos(w*t1) + (1.0/w)*np.sin(w*t1)    
vAct = -x0*w*np.sin(w*t1) + np.cos(w*t1)
# get errors between numerical and analytic solutions
Err1 = np.abs(xAct-np.array(x1))
Err2 = np.abs(xAct-np.array(x2))

# plot position vs time for both methods
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t1,x1,'b.-', label='Position')
plt.plot(t1,v1,'r.-', label='Velocity')
plt.plot(t1,E1,'g.-', label='Energy')
plt.ylim(-1.5,1.5)
plt.legend(loc=3)
plt.subplot(2,1,2)
plt.plot(t2,x2,'b.-', label='Position')
plt.plot(t2,v2,'r.-', label='Velocity')
plt.plot(t2,E2,'g.-', label='Energy')
plt.ylim(-1.5,1.5)
plt.xlabel('Time')

# plot errors between numerical and analytic solutions
plt.figure(2)
plt.subplot(2,1,1)
plt.plot(t1,xAct,'r.-', label='Analytic Solution')
plt.plot(t1,x1,'b.-')
plt.plot(t1,x2,'g.-')
plt.ylabel('Position')
plt.legend(loc=3)
plt.subplot(2,1,2)
plt.plot(t1,Err1,'b.-', label='Simple Euler')
plt.plot(t1,Err2,'g.-', label='Modified Euler')
plt.xlabel('Time')
plt.ylabel('Absolute Error')
plt.ylim(-0.1,0.5)
plt.legend(loc=2)

# plot phase space velocity vs position
plt.figure(3)
plt.plot(x1,v1,'b', label='SHM Phase Space')
plt.xlabel('Position')
plt.ylabel('Velocity')
plt.legend(loc=9)

# plot potential energy vs position
plt.figure(4)
plt.plot(x1,PE1,'b', label='SHM Potential')
plt.xlabel('Position')
plt.ylabel('Potential Energy')
plt.legend(loc=9)