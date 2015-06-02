"""
Created on Sun Oct 27 17:50:45 2013
"""

# Percolation Cluster
"""
This generates a percolation cluster lattice.
"""

import numpy as np
import matplotlib.pyplot as plt
import ClusterGrowth as cluster


# call function for Percolation growth model
#-------------------------------------------------------------------
# initial conditions
N = 50
p = 0.6
# call function
Ox,Oy,Px,Py = cluster.Percolation(N, p)

stack = np.column_stack((Ox,Oy))
#np.savetxt('PercolationCluster_'+str(N)+'n'+str(p)+'p.txt',stack,fmt='%i')

# set marker size from derived visual scaling
msize = int(128325*(N**-2.179))

# plot Percolation growth cluster
plt.scatter(Ox,Oy,s=msize,marker='s',c='b',edgecolor='k')
plt.scatter(Px,Py,s=msize,marker='s',c='none',edgecolor='k')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(-1.0,N)
plt.ylim(-1.0,N)
#-------------------------------------------------------------------