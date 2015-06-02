"""
Created on Sun Nov 17 21:46:06 2013
PHYS 613, Project 3 
Nick Crump
"""

# Anomalous Diffusion
"""
This plots the relationship of mean displacement vs time
of a random walk over a spanning cluster. 
"""

import numpy as np
import random as rand
import matplotlib.pyplot as plt
import ClusterGrowth as cluster
import ClusterFind as find
from datetime import datetime

# problem 3
# ------------------------
t1 = datetime.now()

# start conditions
N = 200
p = 0.60
steps = 40000
iters = 15

mid = int(0.5*N)
pos = [mid,mid]

step = range(steps)
disp = []
x = [mid]
y = [mid]

# call my module to generate percolation cluster
Ox,Oy,Px,Py = cluster.Percolation(N,p)
# call my module to find spanning cluster
clusters,spanning,labels,spanLab = find.ClusterFindFast(Ox,Oy,N,iters,msg='on')

# start at middle of spanning cluster for walk
if spanning > 0:
    while labels[pos[0],pos[1]] != spanLab:
        pos[1] = pos[1]+1
        
    for i in step:
        # get sites left, right, up, down
        lt = [pos[0], pos[1]-1]
        rt = [pos[0], pos[1]+1]
        up = [pos[0]+1, pos[1]]
        dn = [pos[0]-1, pos[1]]
        
        # get labels of neighboring sites
        dirs = [lt,rt,up,dn]
        labs = [labels[lt[0]][lt[1]],labels[rt[0]][rt[1]],labels[up[0]][up[1]],labels[dn[0]][dn[1]]]
        
        # check which neighboring sites are part of cluster
        spnPos = np.where(labs == spanLab)[0]
        lenPos = len(spnPos)
         
        # take random step on spanning cluster
        if lenPos > 0:
            stpLab = rand.sample(spnPos,1)[0]
            stpPos = dirs[stpLab]
            
            pos = stpPos
            x.append(pos[0])
            y.append(pos[1])
            dist = ((mid - pos[0])**2 + (mid - pos[1])**2)**0.5
            disp.append(dist)
        
        else:
            print 'breakout'
            break 
    
    # get logs of both sides
    logS = np.log(step)
    logD = np.log(disp)
    
    # calculate fractal dimension
    k = round(logD[len(logD)-1]/logS[len(logS)-1],3)  
    print 'k = ',k, ' z = ',1.0/k
    
    # plot random walk path on percolating cluster
    plt.figure()
    plt.plot(x,y,'b-', label='Anomalous Diffusion')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(loc=2)    
    
    # plot mean displacement vs step
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(step,disp,'b-', label='Mean Displacement')
    plt.xlabel('Step')
    plt.ylabel('Mean Displacement')
    plt.legend(loc=2)
    
    # plot log-log of mean displacement vs step
    plt.subplot(2,1,2)
    plt.plot(logS,logD,'r-', label='Log-Log Plot')
    plt.xlabel('log Step')
    plt.ylabel('log Mean Displacement')
    plt.legend(loc=2)

t2 = datetime.now()
print 'runtime = ',t2-t1