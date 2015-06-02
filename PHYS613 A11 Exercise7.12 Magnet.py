"""
Created on Fri Nov 29 21:27:07 2013
PHYS 613, Assignment 9
Nick Crump
"""

# Exercise 7.12
"""
Calculate the magnetic potential of a ring magnet by 
solving the difference equations for the PDEs at the 
boundary conditions of the magnet. 
"""

import numpy as np
from matplotlib import cm
from math import sqrt, cos, pi
import matplotlib.pyplot as plt


# set magnet parameters
# --------------------------------------------
M = 1e5/1000  # magnetization (amp/mm)
IR = 5        # magnet inner radius (mm)
OR = 14       # magnet outter radius (mm)
H = 3         # magnet height/2 (mm)
# --------------------------------------------

# set iteration parameters
# --------------------------------------------
tol = 1e-5    # desired relative accuracy
Nx = 25       # x-grid size (mm)
Ny = 25       # y-grid size (mm)
h = 1         # grid step size (mm)
# --------------------------------------------

# set SOR relaxation parameter
# --------------------------------------------
alpha = (4.0/(2+sqrt(4-4*cos(pi/Nx)**2)))-1
# --------------------------------------------


# initialize solution matrix and boundary conditions
# --------------------------------------------
phi = np.zeros((Nx,Ny))

# iterate over xy-grid (i = y-row, j = x-col)
for i in range(Ny):
    for j in range(Nx):
            
        # initialize grid boundary far from magnet (magnetic potential Eq. 7.94)
        if j == Nx-1 or i == Ny-1:
            r = sqrt(i**2 + j**2)
            phi[i][j] = (M*i)/(4*pi*r**3)
        
# start main loop to calculate values
# --------------------------------------------
done = 'no'
iters = 0

while done == 'no':
    done = 'yes'
    iters = iters + 1    
    
    # iterate over xy-grid (i = y-row, j = x-col)
    for i in range(1,Ny-1):
        for j in range(0,Nx-1):
            phi0 = phi[i][j]
            
            # compute along z-axis: (Laplacian Eq. 7.102) 
            if j == 0:
                phi1 = (4*phi[i][1] + phi[i+1][0] + phi[i-1][0])/6.0
                
            # compute outside magnet not on boundary (Laplacian Eq. 7.98)
            elif 0 < j < IR and 0 < i < Ny or OR < j < Nx and 0 < i < Ny:
                a = 1 + 0.5*h/j
                b = 1 - 0.5*h/j
                phi1 = 0.25*(a*phi[i][j+1] + b*phi[i][j-1] + phi[i+1][j] + phi[i-1][j])
            
            # compute on side surfaces (Eq. 7.104)
            elif j == IR and 0 < i < H or j == OR and 0 < i < H:
                phi1 = 0.5*(phi[i][j+1] + phi[i][j-1])
            
            # compute on top surface (Eq. 7.107)
            elif i == H and IR < j < OR:
                phi1 = 0.5*(M*h + phi[i+1][j] + phi[i-1][j])
                
            # compute at corners of magnet cross-section (Eq. 7.108)
            elif i == H and j == IR or i == H and j == OR:
                phi1 = 0.25*(M*h + phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1])
                
                            
            # compute improved SOR approximation and error
            phiSOR = phi1 + alpha*(phi1 - phi0)        
            err = abs(phiSOR-phi1)/phiSOR
            phi[i][j] = phiSOR
            
            # if desired tolerance met then stop
            if err > tol:
                done = 'no'

print '\n','iterations = ',iters
# --------------------------------------------


# set contour levels for plotting
phiMax = np.max(phi)
levels = np.arange(0,int(phiMax),0.5)
print 'max potential =',phiMax

# build matrix for plotting full cross-section
# this reflects matrix in all 4 quadrants
z = np.zeros((2*Nx,2*Nx))
for i in range(Nx):
    for j in range(Nx):
        z[i][Nx-1+j] = phi[Nx-1-i][j]
        z[i+Nx-1][j+Nx-1] = phi[i][j]
        z[i+Nx-1][j] = phi[i][Nx-1-j]
        z[Nx-1-i][j] = phi[i][Nx-1-j]
        
# build x,y meshgrids for plotting full cross-section
xpts = range(-Nx+1,Nx+1,1)
ypts = range(-Nx+1,Nx+1,1)
xMesh,yMesh = np.meshgrid(xpts,ypts)


# plot magnetic potential of cross-section in first quadrant
plt.figure()
plt.contourf(phi, cmap=cm.jet,levels=levels)
plt.plot([IR,IR,OR,OR],[0,H,H,0],'k--',linewidth=3)
plt.colorbar()
plt.show()

# plot magnetic potential of full cross-section
plt.figure()
plt.contourf(xMesh,yMesh,z, cmap=cm.jet,levels=levels)
plt.colorbar()
plt.show()