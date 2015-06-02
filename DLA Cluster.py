"""
Created on Sun Oct 27 17:50:45 2013
"""

# DLA Cluster
"""
This generates a DLA particle cluster.
"""

import numpy as np
import matplotlib.pyplot as plt
import ClusterGrowth as cluster
    

# call function for DLA growth model
#-------------------------------------------------------------------
# initial parameters
N = 500
seed = [0,0]
R = 3

# call function
x,y,Px,Py = cluster.DLA(N, seed, R)

# get max point for setting axis limits
mx = [max(np.abs(x)), max(np.abs(y))]
lim = max(mx)+5

# plot DLA growth cluster
plt.plot(x,y,'b.', label='DLA Model')
plt.plot(Px,Py,'m.')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(-lim,lim)
plt.ylim(-lim,lim)
plt.annotate('N = '+str(N),fontsize=14,xy=(0.14,0.84),xycoords='figure fraction')
plt.legend()
#-------------------------------------------------------------------