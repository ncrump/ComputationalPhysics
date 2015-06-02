"""
Created on Sun Nov 17 21:46:06 2013
PHYS 613, Project 3 
Nick Crump
"""

# Clusters vs Occupation Probability
"""
This plots relationship of number of clusters 
as a function of occupation probability.  
"""

import numpy as np
import matplotlib.pyplot as plt
import ClusterGrowth as cluster
import ClusterFind as find
    

# warmup: get clusters vs occupation probability for 50x50 lattice
N = 50
iters = 10

# set occupation probability to loop over
OccProb = np.arange(0,1.001,0.001)
Ncluster = []

for i in OccProb:
    # call my module to generate percolation cluster
    Ox,Oy,Px,Py = cluster.Percolation(N,i)
    # call my module to find spanning cluster
    clusters,spanning = find.ClusterFindFast(Ox,Oy,N,iters)
    # get number of clusters for p
    Ncluster.append(clusters)
    
# plot trend
plt.plot(OccProb,Ncluster,'b.', label='50 x 50 Lattice')
plt.xlabel('Site Occupation Probability')
plt.ylabel('Number of Clusters')
plt.legend()