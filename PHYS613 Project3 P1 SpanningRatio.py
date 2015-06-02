"""
Created on Sun Nov 17 21:46:06 2013
PHYS 613, Project 3 
Nick Crump
"""

# Spanning Ratio
"""
This plots relationship of ratio of sites in spanning cluster
to total occupied sites as a function of occupation probability
for a given cluster size.  
"""

import numpy as np
import matplotlib.pyplot as plt
import ClusterGrowth as cluster
import ClusterFind as find
from datetime import datetime

# problem 1
# ------------------------
t1 = datetime.now()

# lattice sizes and occupation probabilities to loop over
N = [50,100,200]
iters = [10,10,10]
OccProb = np.arange(0.54,0.802,0.02)
Pinfty = []

# loop over lattice sizes
for i in range(3):
    Pinfty0 = []
    
    # loop through probabilities
    for j in OccProb:
        temp = []
        
        # average over this many trials
        for k in range(10):
            # call my module to generate percolation cluster
            Ox,Oy,Px,Py = cluster.Percolation(N[i],j)
            # call my module to find spanning cluster
            clusters,spanning,labels,spanLab = find.ClusterFindFast(Ox,Oy,N[i],iters[i])
            # get P_infinity ratio
            temp.append(spanning/(N[i]*N[i]*j))
        
        # average over trials
        Pinfty0.append(np.average(temp))
    
    # get array of P_infinity values for N
    Pinfty.append(Pinfty0)
    
# plot trends Pinfty vs p
plt.figure()
plt.plot([0.593,0.593],[0,1],'k--')
plt.plot(OccProb,Pinfty[0],'b.-', label='50x50')
plt.plot(OccProb,Pinfty[1],'r.-', label='100x100')
plt.plot(OccProb,Pinfty[2],'g.-', label='200x200')
plt.xlabel('Site Occupation Probability')
plt.ylabel('$P_{\infty}$',fontsize=16)
plt.xlim(0.4,1.0)
plt.ylim(0,1)
plt.legend(loc=2)

t2 = datetime.now()
print 'runtime = ',t2-t1