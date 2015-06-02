"""
Created on Thu Oct 10 23:37:58 2013
PHYS 613, Assignment 6
Nick Crump
"""

# Problem 1 (EBB6)
"""
Simulation of a random walk in 2-Dimensions on a lattice of evenly spaced points.  
"""

import numpy as np
import matplotlib.pyplot as plt


# Part 1
# calculate mean displacement of random walker on a square 2-D grid
#*******************************************************************
# define range of N steps to compute random walk for and number of sample walks per N
Nsteps = range(1000,4200,200)
Nwalks = 100

Mdisp = []  # stores mean displacement for each N over number of samples

# loop to iterate over random walk steps and samples to run per step
for i in Nsteps:    
    disp = []  # stores displacements for each N per walk
   
   # run the random walk for number of samples per N
    for j in range(Nwalks):
        # define start position at origin of 2-D x,y grid
        pos = [0,0]
        
        # store x,y positions for plotting path
        xp = [0]
        yp = [0]
                
        # generate array of random uniform numbers between [0,1)
        rand = np.random.uniform(0,1,i)
        
        # loop through array of random numbers for random walk
        for k in range(i):
            # depending on random number go right, left, up, or down by 1 unit
            if 0.00 <= rand[k] < 0.25: pos[0] = pos[0]+1
            if 0.25 <= rand[k] < 0.50: pos[0] = pos[0]-1
            if 0.50 <= rand[k] < 0.75: pos[1] = pos[1]+1
            if 0.75 <= rand[k] < 1.00: pos[1] = pos[1]-1
                
            # append x,y position for plotting path
            xp.append(pos[0])
            yp.append(pos[1])
        
        # get total displacement traveled and store to array
        dist = (pos[0]**2 + pos[1]**2)**0.5
        disp.append(dist)
    
    # get mean displacement for each N steps over 100 samples
    mean = np.average(disp)
    Mdisp.append(mean)
    
# relation should be a power function so take logs of both sides
x = np.log(np.array(Nsteps))
y = np.log(np.array(Mdisp))
# calculate slope of log-log plot to get power coefficient
rise = y[len(y)-1] - y[0]
run = x[len(x)-1] - x[0]
slope = round(rise/run,3)


# plot random walk path
plt.figure()
plt.plot(xp,yp,'b',label='2-D Random Walk')
plt.xlabel('x')
plt.ylabel('y')
plt.annotate('1000 Steps',fontsize=15,xy=(0.15,0.78),xycoords='figure fraction')
plt.legend(loc=2)

# plot mean displacement vs N steps 
plt.figure()
plt.subplot(211)
plt.plot(Nsteps,Mdisp,'bo-',label='2-D Random Walk')
plt.xlabel('N Steps')
plt.ylabel('Mean Displacement')
plt.annotate('Samples = '+str(Nwalks)+' Walks',fontsize=15,xy=(0.15,0.78),xycoords='figure fraction')
plt.legend(loc=2)
    
# plot log-log of mean displacement vs N steps 
plt.subplot(212)
plt.plot(x,y,'bo-',label='Log-Log Plot')
plt.xlabel('log N Steps')
plt.ylabel('log Mean Displacement')
plt.annotate('Slope = '+str(slope),fontsize=15,xy=(0.15,0.34),xycoords='figure fraction')
plt.legend(loc=2)