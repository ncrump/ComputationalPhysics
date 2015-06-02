"""
Created on Sun Oct 27 17:50:45 2013
"""

# Eden Cluster
"""
This generates an Eden particle cluster.
"""

import numpy as np
import matplotlib.pyplot as plt
import ClusterGrowth as cluster


# call function for Eden growth model
#-------------------------------------------------------------------
# initial parameters
N = 500
seed = [0,0]

# call function
x,y = cluster.Eden(N, seed)

# get max point for setting axis limits
mx = [max(np.abs(x)), max(np.abs(y))]
lim = max(mx)+5

# plot Eden growth cluster
plt.plot(x,y,'b.', label='Eden Model')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(-lim,lim)
plt.ylim(-lim,lim)
plt.annotate('N = '+str(N),fontsize=14,xy=(0.14,0.84),xycoords='figure fraction')
plt.legend()
#-------------------------------------------------------------------