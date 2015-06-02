"""
Created on Fri Nov 29 21:27:07 2013
PHYS 613, Assignment 11
Nick Crump
"""

# Exercise 7.11
"""
Calculate the temperature profile over a square plate by solving
the 2-D steady state heat equation as a difference equation using
heat flow through the boundary conditions (Neumann Boundary). 
"""

import numpy as np
from matplotlib import cm
from math import sqrt,cos,pi
import matplotlib.pyplot as plt


# set iteration parameters
# --------------------------------------------
tol = 1e-5    # desired relative accuracy
Nx = 30       # x-grid size 
Ny = 30       # y-grid size
h = 1         # grid step size
# --------------------------------------------

# set SOR relaxation parameter
# --------------------------------------------
alpha = (4.0/(2+sqrt(4-4*cos(pi/Nx)**2)))-1
# --------------------------------------------

# set heat flow boundary conditions
# --------------------------------------------
L = -700.0/Ny   # left side flow (degC/m)
R = -200.0/Ny   # right side flow (degC/m)
B =  400.0/Nx   # bottom side flow (degC/m)
T = -100.0/Nx   # top side flow (degC/m)
Tref = 750.0    # steady state reference temp at lower left corner (degC)
# --------------------------------------------
 

# start main loop to calculate values
# --------------------------------------------
# initialize solution matrix 
temp = np.zeros((Nx,Ny))

done = 'no'
iters = 0

while done == 'no':
    done = 'yes'
    iters = iters + 1    
    
    # iterate over xy-grid (i = y-row, j = x-col)
    for i in range(0,Ny-1):
        for j in range(0,Nx-1):
            temp0 = temp[i][j]
            
            # compute along left edge
            if j == 0:
                temp1 = 0.25*(2*temp[i][1] + temp[i+1][0] + temp[i-1][0] - 2*h*L)
                
            # compute along right edge
            if j == Nx-1:
                temp1 = 0.25*(2*temp[i][Nx-2] + temp[i+1][Nx-1] + temp[i-1][Nx-1] + 2*h*R)
                
            # compute along bottom edge
            if i == 0:
                temp1 = 0.25*(2*temp[1][j] + temp[0][j+1] + temp[0][j-1] - 2*h*B)
                
            # compute along top edge
            if i == Ny-1:
                temp1 = 0.25*(2*temp[Ny-2][j] + temp[Ny-1][j+1] + temp[Ny-1][j-1] + 2*h*T)
                
            # compute at lower left corner
            if i == 0 and j == 0:
                temp1 = 0.5*(temp[0][1] + temp[1][0] - h*B - h*L)
                
            # compute at upper left corner
            if i == Ny-1 and j == 0:
                temp1 = 0.5*(temp[Ny-1][1] + temp[Ny-2][0] + h*T - h*L)
                
            # compute at lower right corner
            if i == 0 and j == Nx-1:
                temp1 = 0.5*(temp[0][Nx-2] + temp[1][Nx-1] + h*R - h*B)
                
            # compute at upper right corner
            if i == Ny-1 and j == Nx-1:
                temp1 = 0.5*(temp[Ny-2][Nx-1] + temp[Ny-1][Nx-2] + h*R + h*T)
                
            # compute everywhere else inside plate
            if 0 < i < Ny-1 and 0 < j < Nx-1: 
                temp1 = 0.25*(temp[i][j+1] + temp[i][j-1] + temp[i+1][j] + temp[i-1][j])
                
            
            # compute improved SOR approximation and error
            tempSOR = temp1 + alpha*(temp1 - temp0)        
            err = abs(tempSOR-temp1)/tempSOR
            temp[i][j] = tempSOR
            
            # shift all temps so that lower left reference temp equals Tref
            diff = Tref - temp[0][0]
            temp = temp + diff
            
            # if desired tolerance met then stop
            if err > tol:
                done = 'no'

print '\n','iterations = ',iters
# --------------------------------------------


# set contour levels for plotting
TMin = np.min(temp)
TMax = np.max(temp)
levels = np.arange(int(TMin),int(TMax),1)
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
plt.contour(xMesh,yMesh,temp, cmap=cm.jet,levels=range(int(TMin),int(TMax),50))
plt.xlabel('x (m)',fontsize=14)
plt.ylabel('y (m)',fontsize=14)
plt.show()