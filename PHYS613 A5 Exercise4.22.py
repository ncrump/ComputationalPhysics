"""
Created on Wed Sep 25 23:14:56 2013
PHYS 613, Assignment 5
Nick Crump
"""

# Exercise 4.22
"""
From Computational Physics by Devries
"""

import Integration as Int
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Exercise 4.22
#*******************************************************************
# arrays to store x, y, and electric potential values
x = []
y = []
pG = []
pM = []

# iterate over xp,yp that define coordinate point for potential
for i in range(2,22,2):
    for j in range(2,22,2):
        
        # create function to integrate for each (xp,yp) point
        f = '1/((x-'+str(i)+')**2 + (y-'+str(j)+')**2)**0.5'
        
        # call gaussLegendre2D method from my integration module
        pGauss = Int.gaussLegendre2D(f,-1.0,1.0,-1.0,1.0,10)
        
        # call monteCarlo2D method from my integration module
        # Note: uses 1000 random points to eval integral
        pMonte = Int.monteCarlo2D(f,-1.0,1.0,-1.0,1.0,1000)
        
        # store values
        x.append(i)
        y.append(j)
        pG.append(pGauss)
        pM.append(pMonte)
        
      
# plot potential computed from 2D Gauss-Legendre as 3D scatter plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot([-1,-1,1,1,-1],[-1,1,1,-1,-1],[0,0,0,0,0],label='2D Gauss-Legendre')
cb = ax.scatter(x,y,pG,pG,s=50,c=pG,cmap=cm.jet)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Electrostatic Potential')
plt.xlim(-2,22)
plt.ylim(-2,22)
plt.colorbar(cb)
plt.legend(loc=9)
plt.show()

# plot potential computed from 2D Monte Carlo as 3D scatter plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot([-1,-1,1,1,-1],[-1,1,1,-1,-1],[0,0,0,0,0],label='2D Monte Carlo')
cb = ax.scatter(x,y,pM,pM,s=50,c=pM,cmap=cm.jet)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Electrostatic Potential')
plt.xlim(-2,22)
plt.ylim(-2,22)
plt.colorbar(cb)
plt.legend(loc=9)
plt.show()
#*******************************************************************