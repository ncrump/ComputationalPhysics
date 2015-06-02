"""
Created on Wed Oct 23 20:55:21 2013
PHYS 613, Assignment 7 
Nick Crump
"""

# Problem 1 (EBB8): Unit Cell Lattice
"""
Generates positions of atoms that form clusters of the following unit cells:
1. Simple Cubic (SC)
2. Body Centered Cubic (BCC)
3. Face Centered Cubic (FCC) 
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Part 1: simple cubic lattice
#*******************************************************************
def SC(Lsides, Ncubes):   
    
    # arrays to store atoms positions
    x = []
    y = []
    z = []
    
    # get x corner locations
    intvl = (Ncubes+1)*Lsides
    xCorners = np.arange(0, intvl, Lsides)    
    
    # loop over line of atoms to build 3D lattice of corners
    for y0 in xCorners:
        for z0 in xCorners:
            # duplicates corner lattice in y,z
            yCNR = np.ones(Ncubes+1)*y0
            zCNR = np.ones(Ncubes+1)*z0
            # write lattice points to arrays
            x = np.append(x,xCorners)
            y = np.append(y,yCNR)
            z = np.append(z,zCNR)
    
    return x,y,z
#*******************************************************************
    
    
# Part 2: body centered cubic lattice
#*******************************************************************
def BCC(Lsides, Ncubes):   
    
    # arrays to store atoms positions
    x = []
    y = []
    z = []
    
    # get x corner locations
    intvl = (Ncubes+1)*Lsides
    xCorners = np.arange(0, intvl, Lsides)    

    # loop over line of atoms to build 3D lattice of corners
    for y0 in xCorners:
        for z0 in xCorners:
            # duplicates corner lattice in y,z
            yCNR = np.ones(Ncubes+1)*y0
            zCNR = np.ones(Ncubes+1)*z0
            # write lattice points to arrays
            x = np.append(x,xCorners)
            y = np.append(y,yCNR)
            z = np.append(z,zCNR)
    
    # ----------------------------------------------        
    # get x center locations
    intvl = Ncubes*Lsides
    xCenters = np.arange(Lsides/2.0, intvl, Lsides)
        
    # loop over line of atoms to build 3D lattice of centers
    for y0 in xCenters:
        for z0 in xCenters:
            # duplicates center lattice in y,z
            yCTR = np.ones(Ncubes)*y0
            zCTR = np.ones(Ncubes)*z0
            # write lattice points to arrays
            x = np.append(x,xCenters)
            y = np.append(y,yCTR)
            z = np.append(z,zCTR)
     
    return x,y,z
#*******************************************************************
    
    
# Part 3: face centered cubic lattice
#*******************************************************************
def FCC(Lsides, Ncubes):   
    
    # arrays to store atoms positions
    x = []
    y = []
    z = []
    
    # get x corner locations
    intvl1 = (Ncubes+1)*Lsides
    xCorners = np.arange(0, intvl1, Lsides)    

    # loop over line of atoms to build 3D lattice of corners
    for y0 in xCorners:
        for z0 in xCorners:
            # duplicates corner lattice in y,z
            yCNR = np.ones(Ncubes+1)*y0
            zCNR = np.ones(Ncubes+1)*z0
            # write lattice points to arrays
            x = np.append(x,xCorners)
            y = np.append(y,yCNR)
            z = np.append(z,zCNR)
    
    # ----------------------------------------------
    # get x face locations
    intvl2 = Ncubes*Lsides
    xFace1 = np.arange(0, intvl1, Lsides)
    xFace2 = np.arange(Lsides/2.0, intvl2, Lsides)
        
    # loop over line of atoms to build 3D lattice of faces in yz-plane
    for y0 in xFace2:
        for z0 in xFace2:
            # duplicates face lattice in y,z
            yFCE = np.ones(Ncubes+1)*y0
            zFCE = np.ones(Ncubes+1)*z0
            # write lattice points to arrays
            x = np.append(x,xFace1)
            y = np.append(y,yFCE)
            z = np.append(z,zFCE)
            
    # loop over line of atoms to build 3D lattice of faces in xy-plane
    for y0 in xFace2:
        for z0 in xFace1:
            # duplicates face lattice in y,z
            yFCE = np.ones(Ncubes)*y0
            zFCE = np.ones(Ncubes)*z0
            # write lattice points to arrays
            x = np.append(x,xFace2)
            y = np.append(y,yFCE)
            z = np.append(z,zFCE)
     
    # loop over line of atoms to build 3D lattice of faces in xz-plane
    for y0 in xFace1:
        for z0 in xFace2:
            # duplicates face lattice in y,z
            yFCE = np.ones(Ncubes)*y0
            zFCE = np.ones(Ncubes)*z0
            # write lattice points to arrays
            x = np.append(x,xFace2)
            y = np.append(y,yFCE)
            z = np.append(z,zFCE)
            
    return x,y,z
#*******************************************************************
# ------------------------------------------------------------------
# Lsides = length of cubic sides
# Ncubes = number of times to replicate unit cube to form lattice
# ------------------------------------------------------------------


x,y,z = FCC(2.2, 6)

# write atom positions to txt file for visualizing lattice in Avogadro
# label with '20' for Calcium 
atoms = np.ones(len(x))*20
stack = np.column_stack((atoms,x,y,z))
np.savetxt('FaceCenteredCubic.txt',stack,fmt=('%i','%f4','%f4','%f4'))

# plot atom positions
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x,y,z,'bo',markersize=8)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
