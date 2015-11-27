"""
this generates a Mandelbrot set using a simple slow method
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# input parameters
# ------------------------------
# set xy-domain size
xmin = -1.8
xmax =  0.8
ymin = -1.0
ymax =  1.0

# set xy-grid size
xpts = 200
ypts = 200

# set max iterations
maxn = 20
# ------------------------------

# iterate over Mandelbrot equation
# z = z^2 + c for complex numbers c and z

# start timer
t0 = datetime.now()

# initialize variables
im = np.complex(0,1)
z0 = np.complex(0,0)
x  = []
y  = []

# get grid
xgrid = np.linspace(xmin,xmax,xpts)
ygrid = np.linspace(ymin,ymax,ypts)

# evaluate grid points
for i in xgrid:
    for j in ygrid:
        c = i + im*j
        z = z0
        k = 0
        # if abs(z) stays less than about 2 then point is in set
        while abs(z) < 2 and k < maxn:
            z = z*z + c
            k += 1
        if k >= maxn:
            x.append(i)
            y.append(j)

# make plot
plt.plot(x,y,'.k')

# print elapsed time
t1 = datetime.now()
print '\nelapsed time:'
print t1-t0