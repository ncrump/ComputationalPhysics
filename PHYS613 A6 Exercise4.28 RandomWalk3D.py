"""
Created on Thu Oct 10 23:37:58 2013
PHYS 613, Assignment 6
Nick Crump
"""

# Exercise 4.28
# Exercise 4.29
"""
From Computational Physics by Devries
"""

import numpy as np
from math import pi,sin,cos,acos
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Exercise 4.28
#*******************************************************************
# define initial conditions
velocity = 500      # m/sec
meanfreepath = 1.0  # meters
Nwalks = 100        # number of molecules to average over

# stores displacements after each collision
disp = []          
# stores final positions of molecules for plotting position distribution
mx = []
my = []
mz = []

# define number of collisions molecule with encounter after 1 sec
collisions = int(velocity/meanfreepath)
       
# run the random walk this many times
for i in range(Nwalks):
    # define start position [r,theta,phi] at origin of sphere
    pos = [0,0,0]  
    
    # stores [x,y,z] positions of molecule path for plotting in sphere
    xpos = []    
    ypos = []
    zpos = []
    
    # loop to generate random positions to step molecule through path
    for j in range(collisions):
        # generate a random uniform number between [0,1) for r (meters)
        r = np.random.uniform(0,1,1)
        # generate a random uniform number between [-1,1) for theta (radians)
        theta0 = np.random.uniform(-1,1,1)
        theta = acos(theta0)
        # generate a random uniform number between [0,2pi) for phi (radians)
        phi = np.random.uniform(-2*pi,2*pi,1)

        # check position r to move a constant +/- 1 meter in radius
        if 0.0 <= r < 0.5: pos[0] = pos[0]+1
        if 0.5 <= r < 1.0: pos[0] = pos[0]-1
      
        # move molecule in (theta,phi) direction from random number generated
        pos[1] = pos[1]+theta
        pos[2] = pos[2]+phi
        
        # convert positions from spherical [r,theta,phi] to cartesian [x,y,z]
        # for plotting path of molecule inside sphere
        xpos.append(pos[0]*sin(pos[1])*cos(pos[2]))
        ypos.append(pos[0]*sin(pos[1])*sin(pos[2]))
        zpos.append(pos[0]*cos(pos[1]))
        
        
    # convert final displacement from spherical [r,theta,phi] to cartesian [x,y,z]
    x = pos[0]*sin(pos[1])*cos(pos[2])
    y = pos[0]*sin(pos[1])*sin(pos[2])
    z = pos[0]*cos(pos[1])
    
    # store final molecule positions for plotting position distribution
    mx.append(x)
    my.append(y)
    mz.append(z)
    
    # get total displacement traveled and store to array
    dist = (x**2 + y**2 + z**2)**0.5
    disp.append(dist)

# get mean displacement over Nwalks
mean = np.average(disp)
print 'Mean Displacement = ',mean, ' meters'


# plot 3-D diffusive path of molecule inside sphere
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(xpos,ypos,zpos,'bo',linestyle='--',linewidth=0.5,label='Diffusive Path')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.legend()
#*******************************************************************


# Plots for Exercise 4.29
#*******************************************************************
# plot molecule distribution after each random walk
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(mx,my,mz,'bo',label='Molecule Distribution')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.legend()

# plot histogram of net displacements representing many different molecule paths
plt.figure()
plt.hist(disp,bins=10,label=str(Nwalks)+' Molecules')
plt.xlabel('Net Displacement (m)')
plt.ylabel('Frequency')
plt.legend()
#*******************************************************************