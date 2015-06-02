"""
Generates random uniformly distributed points on the surface of a unit sphere
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#****************************************************
# number of points to generate
n = 1000

# generates randomly distributed theta,phi points at r = 1
r = np.ones(n)
# this prevents clumping at poles for theta,phi
theta0 = np.random.uniform(-1,1,n)    
theta = np.arccos(theta0)
#theta = np.random.uniform(0,2*np.pi,n) 
phi = np.random.uniform(0,2*np.pi,n)

# converts spherical [r,theta,phi] points to cartesian [x,y,z] points
x = r*np.sin(theta)*np.cos(phi)
y = r*np.sin(theta)*np.sin(phi)
z = r*np.cos(theta)

# plots distribution of points
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z,'bo')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
#****************************************************