"""
Created on Fri Nov 29 21:27:07 2013
PHYS 613, Assignment 11
Nick Crump
"""

# Exercise 7.8
"""
Calculate the temperature profile over a square plate by solving
the 2-D steady state heat equation as a difference equation using
fixed temperature boundary conditions (Dirichlet Boundary). 
"""

import numpy as np
from matplotlib import cm
from math import sqrt,cos,pi
import matplotlib.pyplot as plt


# set iteration parameters
# --------------------------------------------
tol = 1e-5    # desired relative accuracy
Nx = 33       # x-grid size 
Ny = 33       # y-grid size
h = 1         # grid step size
# --------------------------------------------

# set SOR relaxation parameter
# --------------------------------------------
alpha = (4.0/(2+sqrt(4-4*cos(pi/Nx)**2)))-1
# --------------------------------------------


# initialize solution matrix and boundary conditions
# --------------------------------------------
temp = np.zeros((Nx,Ny))

# iterate over xy-grid (i = y-row, j = x-col)
for i in range(Nx):
    for j in range(Ny):
        
        # fix temp at 100C along left and bottom edges of plate
        # leave temp at 0C along right and top edges of plate
        if i == 0 or j == 0:
            temp[i][j] = 100   

# start main loop to calculate values
# --------------------------------------------
done = 'no'
iters = 0

while done == 'no':
    done = 'yes'
    iters = iters + 1    
    
    # iterate over xy-grid (i = y-row, j = x-col)
    for i in range(1,Ny-1):
        for j in range(1,Nx-1):
            temp0 = temp[i][j]
            
            # compute temperature using finite difference            
            temp1 = 0.25*(temp[i][j+1] + temp[i][j-1] + temp[i+1][j] + temp[i-1][j])
                            
            # compute improved SOR approximation and error
            tempSOR = temp1 + alpha*(temp1 - temp0)        
            err = abs(tempSOR-temp1)/tempSOR
            temp[i][j] = tempSOR
            
            # if desired tolerance met then stop
            if err > tol:
                done = 'no'

print '\n','iterations = ',iters
# --------------------------------------------


# set contour levels for plotting
TMax = np.max(temp)
levels = np.arange(0,int(TMax),0.2)
print 'max temp =',TMax
        
# build x,y meshgrids for plotting axes
xpts = np.linspace(0,1,Nx)
ypts = np.linspace(0,1,Ny)
xMesh,yMesh = np.meshgrid(xpts,ypts)

# plot temperature distribution over square plate
plt.figure()
plt.contourf(xMesh,yMesh,temp, cmap=cm.jet,levels=levels)
plt.xlabel('x (m)',fontsize=14)
plt.ylabel('y (m)',fontsize=14)
plt.colorbar()
plt.show()

# plot isotherms - contour lines of constant temperature
plt.figure()
plt.contour(xMesh,yMesh,temp, cmap=cm.jet,levels=range(0,100,10))
plt.xlabel('x (m)',fontsize=14)
plt.ylabel('y (m)',fontsize=14)
plt.annotate('$90$',fontsize=18,xy=(0.19,0.19))
plt.annotate('$80$',fontsize=18,xy=(0.28,0.28))
plt.annotate('$70$',fontsize=18,xy=(0.35,0.35))
plt.annotate('$60$',fontsize=18,xy=(0.42,0.42))
plt.annotate('$50$',fontsize=18,xy=(0.48,0.48))
plt.annotate('$40$',fontsize=18,xy=(0.54,0.54))
plt.annotate('$30$',fontsize=18,xy=(0.60,0.60))
plt.annotate('$20$',fontsize=18,xy=(0.67,0.67))
plt.annotate('$10$',fontsize=18,xy=(0.77,0.77))
plt.show()