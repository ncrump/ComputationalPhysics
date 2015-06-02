"""
Created on Fri Nov 29 21:27:07 2013
PHYS 613, Assignment 11
Nick Crump
"""

# Exercise 7.10
"""
Calculate the temperature profile over a octagonal plate by solving
the 2-D steady state heat equation as a difference equation using 
fixed temperature boundary conditions (Irregular Boundary). 

NOTE: this one needs to be fixed - not a good method
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
R = 10        # radius of circular plate
# --------------------------------------------

# set SOR relaxation parameter
# --------------------------------------------
alpha = (4.0/(2+sqrt(4-4*cos(pi/Nx)**2)))-1
# --------------------------------------------


# initialize solution matrix and boundary conditions
# --------------------------------------------
temp = np.zeros((Nx,Ny))

# get circle center
x0 = int(0.5*Nx)
y0 = int(0.5*Ny)

# iterate over xy-grid (i = y-row, j = x-col)
for i in range(Nx):
    for j in range(Ny):
        
        # check if points on circumference of circular plate
        if (R-0.15)**2 <= (j-x0)**2 + (i-y0)**2 <= (R+0.6)**2:

            # fix temp at 100C around one quarter of plate
            # leave temp at 0C around rest of plate
            if 0 <= i <= 15 and 0 <= j <= 15:
                temp[i][j] = 100


# set manually since radius check misses them 
temp[5][11] = 100
temp[11][5] = 100

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
            
            # check if points inside circular plate
            if (j-x0)**2 + (i-y0)**2 < (R-0.15)**2:

                temp0 = temp[i][j]
                
                # if point outside radius then use average of two points 
                if sqrt((j+1-x0)**2 + (i-y0)**2) > R:
                    A = 0.5*(temp[i][j] + temp[i][j+1])
                else:
                    A = temp[i][j+1]
                
                if sqrt((j-1-x0)**2 + (i-y0)**2) > R:
                    B = 0.5*(temp[i][j] + temp[i][j-1])                
                else:
                    B = temp[i][j-1]

                if sqrt((j-x0)**2 + (i+1-y0)**2) > R:
                    C = 0.5*(temp[i][j] + temp[i+1][j])
                else:
                    C = temp[i+1][j]
                    
                if sqrt((j-x0)**2 + (i-1-y0)**2) > R:
                    D = 0.5*(temp[i][j] + temp[i-1][j])
                else:
                    D = temp[i-1][j]
                
                # do finite difference
                temp1 = 0.25*(A + B + C + D)
                                
                # compute improved SOR approximation and error
                tempSOR = temp1 + alpha*(temp1 - temp0)        
                err = abs(tempSOR-temp1)/tempSOR
                temp[i][j] = tempSOR
                
                # if desired tolerance met then stop
                if err > tol:
                    done = 'no'
               
# set manually since radius check gets too many
temp[15][5] = 0
temp[14][5] = 0
temp[13][5] = 0
temp[12][5] = 0
temp[11][5] = 0
temp[5][11] = 0
temp[5][12] = 0
temp[5][13] = 0
temp[5][14] = 0
temp[5][15] = 0


print '\n','iterations = ',iters
# --------------------------------------------


# set contour levels for plotting
TMax = np.max(temp)
levels = np.arange(0,int(TMax),0.5)
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
plt.show()