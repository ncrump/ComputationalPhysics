"""
Created on Sun Nov 17 21:46:06 2013
PHYS 613, Project 3 
Nick Crump
"""

# Box Counting Method
"""
This plots the relationship of the number of occupied sites
in a cluster when p=pc as a function of box size to obtain 
the fractal dimension of spanning clusters.  
"""

import numpy as np
import matplotlib.pyplot as plt
import ClusterGrowth as cluster
import ClusterFind as find
from datetime import datetime

# problem 2
# ------------------------
t1 = datetime.now()

# start conditions
N = 600
p = 0.60
step = 5
iters = 40

Lsize = []
Nsite = []

# call my module to generate percolation cluster
Ox,Oy,Px,Py = cluster.Percolation(N,p)
# call my module to find spanning cluster
clusters,spanning,labels,spanLab = find.ClusterFindFast(Ox,Oy,N,iters,msg='on')

if spanning > 0:
    for i in range(0, int(0.5*N), int(step)):
        # get number of occupied sites in smaller and smaller boxes
        count = len(np.where(labels[i:N-i+1, i:N-i+1] > 0)[0])
        side = N-i-i
        Lsize.append(side)
        Nsite.append(count)
        
    # reverse order of arrays for plotting
    Lsize = Lsize[::-1]
    Nsite = Nsite[::-1]
    
    # get logs of both sides
    logL = np.log(Lsize)
    logN = np.log(Nsite)
    
    # calculate fractal dimension
    df = round(logN[len(logN)-1]/logL[len(logL)-1],3)  
    print 'df = ',df      
    
    # plot M vs l trend - should follow power law M = l**b
    plt.subplot(2,1,1)
    plt.plot(Lsize,Nsite,'b.-', label='$M\ (l)$')
    plt.xlabel('Box Length $l$')
    plt.ylabel('Number of Sites $M$')
    plt.annotate('$600x600$',fontsize=18,xy=(0.15,0.75),xycoords='figure fraction')
    plt.legend(loc=2)
    
    # plot log-log trend - should be line with slope b = 1.9
    plt.subplot(2,1,2)
    plt.plot(logL,logN,'r.-', label='Log-Log')
    plt.xlabel('$log(l)$',fontsize=16)
    plt.ylabel('$log(M)$',fontsize=16)
    plt.annotate('$d_{f} = 1.915$',fontsize=18,xy=(0.15,0.32),xycoords='figure fraction')
    plt.legend(loc=2)

t2 = datetime.now()
print 'runtime = ',t2-t1