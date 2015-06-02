"""
Created on Thu Oct 10 23:37:58 2013
PHYS 613, Assignment 6
Nick Crump
"""

# Exercise 4.30
# Exercise 4.31
"""
From Computational Physics by Devries
"""

import numpy as np
from math import pi,sin,cos,acos
import matplotlib.pyplot as plt


# Exercise 4.30
#*******************************************************************
# define number of collisions to run and paths to average over
collisions = range(100,2600,200) 
paths = 100    

# stores mean displacements over N paths
Mdisp1 = []   

# loop through collisions
for i in collisions:
    # stores mean displacements over N collisions
    Mdisp0 = [] 
    
    # loop through number of paths to average over
    for k in range(paths):
        # define start position [r,theta,phi] at origin of sphere
        pos = [0,0,0]  
        
        # stores displacements after each collision
        disp = [] 
    
        # loop to generate random positions to step molecule through path
        for j in range(i):
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
            
            # convert final displacement from spherical [r,theta,phi] to cartesian [x,y,z]
            x = pos[0]*sin(pos[1])*cos(pos[2])
            y = pos[0]*sin(pos[1])*sin(pos[2])
            z = pos[0]*cos(pos[1])
        
            # get total displacement traveled and store to array
            dist = (x**2 + y**2 + z**2)**0.5
            disp.append(dist)

        # get mean displacement over N collisions
        mean0 = np.average(disp)
        Mdisp0.append(mean0)
    
    # get mean displacement over N paths
    mean1 = np.average(Mdisp0)
    Mdisp1.append(mean1)
    
# relation should be a power function so take logs of both sides
logx = np.log(np.array(collisions))
logy = np.log(np.array(Mdisp1))
# calculate slope of log-log plot to get power coefficient
rise = logy[len(logy)-1] - logy[0]
run = logx[len(logx)-1] - logx[0]
slope = round(rise/run,3)


# plot mean displacement vs N collisions 
plt.figure()
plt.subplot(211)
plt.plot(collisions,Mdisp1,'bo-',label='Mean Displacement')
plt.xlabel('N Collisions')
plt.ylabel('Mean Displacement')
plt.legend(loc=2)
    
# plot log-log of mean displacement vs N collisions
plt.subplot(212)
plt.plot(logx,logy,'bo-',label='Log-Log Plot')
plt.xlabel('log N Collisions')
plt.ylabel('log Mean Displacement')
plt.annotate('Slope = '+str(slope),fontsize=15,xy=(0.15,0.34),xycoords='figure fraction')
plt.legend(loc=2)