"""
Generates concentric circles of random uniformly distributed points. 
"""

import numpy as np
import matplotlib.pyplot as plt


#****************************************************
# number of points to generate on each circle
# number of concentric circles to generate
nPoints = 2000
nCircles = 10

for i in np.linspace(0.1,1.0,nCircles):
    # generates randomly distributed theta points at various r
    r = np.ones(nPoints)*i
    theta = np.random.uniform(0,2*np.pi,nPoints) 
    
    # converts polar coords points to cartesian coords points
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    
    # plots distribution of points
    plt.plot([0],[0],'ro',markersize=10)
    plt.plot(x,y,'b.')
    plt.xlim(-1.1,1.1)
    plt.ylim(-1.1,1.1)
#****************************************************