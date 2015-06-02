"""
Created on Wed Oct 02 18:53:32 2013
PHYS 613, Project 1
Nick Crump
"""

# Project 1: Minimum Energy Molecules
"""
Finding the optimal atomic positions that minimize the interatomic potential 
energy of Potassium molecules (K). 
"""

import numpy as np
from math import sqrt,exp
from scipy import optimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # used for 3D plotting


# Part 2
# define function for interatomic potential Ucoh = Urep - Uel
# returns value of multivariable function of 3N variables for N atoms    
#*******************************************************************
# ------------------------------------------------------------------
# input argument 'pos0' is an array of initial positions 
# 'pos0' is input as [x1,y1,z1,x2,y2,z2 ... xn,yn,zn]
# ------------------------------------------------------------------
def P(pos0):
    # constants for Potassium molecule (K)
    eps = (1.51/1000)*(27.21138386/2.0)  # (eV)
    gam = (19.3/1000)*(27.21138386/2.0)  # (eV)
    r0 = (8.253*0.52917721)              # (angstroms)
    p = 10.58
    q = 1.34
    
    # get number of atoms for loop from input positions
    n = int(len(pos0)/3)
    
    # build array of initial positions of atoms in (x0,y0,z0)
    x0 = [pos0[k] for k in range(0, 3*n, 3)]
    y0 = [pos0[k] for k in range(1, 3*n, 3)]
    z0 = [pos0[k] for k in range(2, 3*n, 3)]
    
    # terms for iteration
    Urep = 0
    Uel = 0
    
    # iterate through double sum of energy function
    for i in range(n):
        # term for iteration
        Uel0 = 0
        # get initial xi positions
        xi,yi,zi = x0[i],y0[i],z0[i]
        
        for j in range(n):
            if j != i:
                # get initial xj positions
                xj,yj,zj = x0[j],y0[j],z0[j]  
                
                # construct 3-D distance equation
                rij = sqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
                
                # calculate repulsive and attractive energy terms
                Urep = Urep + eps*(exp(-p*((rij/r0)-1.0)))
                Uel0 = Uel0 + (gam**2.0)*(exp(-2.0*q*((rij/r0)-1.0)))
        
        # take care of sqrt in sum over j
        Uel = Uel + sqrt(Uel0)
    
    # calculate total cohesive potential energy
    Ucoh = Urep - Uel
    
    # return value for interatomic potential
    return Ucoh
#******************************************************************* 
    

# find [x,y,z] positions of atoms that globally minimize the interatomic potential
# this determines a likely configuration of atoms for a minimum energy molecule
# ------------------------------------------------------------------
# set initial parameters
Natoms = 3                 # number of atoms to initialize
Nrandom = 01               # number of iterations to run

r0 = (8.253*0.52917721)    # atomic bond distance (angstroms)
Pmin = []                  # stores minima of function (eV)
Ppos = []                  # stores atom positions at function minima (angstroms)

# randomly generate N atoms within a cube of side length 2.5*r0
# this initializes N atoms with randomly selected initial integer positions 
for i in range(Nrandom):
    # initial positions of atoms as [x1,y1,z1,x2,y2,z2 ... xn,yn,zn]
    pos0 = np.random.randint(0,int(2.5*r0),3*Natoms)
    
    # use scipy minimization method for nonlinear conjugate gradient
    # returns function minimum and the parameter values that minimize the function
    # calls function to calculate interatomic potential above
    minPos,minFun,fc,gc,wf = optimize.fmin_cg(P,pos0,full_output=1,disp=0)
    
    # keep energies and positions when new minima are found
    if round(minFun,4) not in Pmin:
        Pmin.append(round(minFun,4))
        Ppos.append(minPos)
        
print 'minima = ',Pmin
# ------------------------------------------------------------------


## set up for 3-D plotting of atoms in minimum energy configuration
# also return distances between atoms
# ------------------------------------------------------------------
# break arary of returned position values into x,y,z arrays for plotting
n = int(len(minPos)/3)
x = [minPos[k] for k in range(0, 3*n, 3)]  # gets array of minimized x coords
y = [minPos[k] for k in range(1, 3*n, 3)]  # gets array of minimized y coords
z = [minPos[k] for k in range(2, 3*n, 3)]  # gets array of minimized z coords

# connect the dots
x.append(minPos[0])
y.append(minPos[1])
z.append(minPos[2])

# plot positions of atoms as 3-D scatter plot  
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter(x,y,z,s=200)
ax.plot(x,y,z,label='N='+str(n)+' atoms')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.legend()
plt.show()

# loop through positions to return distances between atoms
for i in range(n):
    for j in range(n):
        if j > i:
            print 'distance = ',sqrt((x[i]-x[j])**2 + (y[i]-y[j])**2 + (z[i]-z[j])**2)
# ------------------------------------------------------------------


# collect all minima found above from N = 3 to N = 7 for plotting
# plot minima along with global minimums for each N
# ------------------------------------------------------------------
# repeatable minima found over many iterations
min3R = [0,-0.4919]
min4R = [0,-0.9922,-0.4919,-1.4746]
min5R = [-1.5734,-0.9922,-1.4841]
min6R = [-1.5734,-2.1479,-2.0652,-0.9922]
min7R = [-3.3628,-2.1479,-2.7891]

# non-repeatable minima found over many iterations
min5N = [0,-1.9510]
min6N = [0,-1.4746,-0.4919]
min7N = [-1.9759,-1.4841,-0.9922]

# global minimums for each N
minGlobal = [-0.9922,-1.5734,-2.1479,-2.7891,-3.394]

plt.figure()
# plot repeatable minima in blue
plt.plot([3,3],min3R,'bo',label='Repeatable Minima')
plt.plot([4,4,4,4],min4R,'bo')
plt.plot([5,5,5],min5R,'bo')
plt.plot([6,6,6,6],min6R,'bo')
plt.plot([7,7,7],min7R,'bo')

# plot non-repeatable minima in red
plt.plot([5,5],min5N,'rx',label='Non-Repeatable Minima')
plt.plot([6,6,6],min6N,'rx')
plt.plot([7,7,7],min7N,'rx')

# plot global minimums in blue with connecting line
plt.plot([3,4,5,6,7],minGlobal,'bo-',markersize=8,label='Global Minimum')
plt.xlabel('N Atoms')
plt.ylabel('Energy Minima')
plt.xlim(2,8)
plt.ylim(-3.6,0.1)
plt.legend(loc=3)