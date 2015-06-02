"""
Example of how to use the scipy optimization tools for
finding the global minimum of a 2-dimensional function.
"""

from scipy import optimize
import numpy as np


# define 2-dimensional function to minimize in (u,v)
# this is an elliptic paraboloid
# **************************************************
def f(x, *args):
    u, v = x
    a, b, c, d, e, f = args
    return a*u**2 + b*u*v + c*v**2 + d*u + e*v + f
# **************************************************

# define gradient of function (df/du), (df/dv)
# this step is optional 
# **************************************************
def gradf(x, *args):
    u, v = x
    a, b, c, d, e, f = args
    gu = 2*a*u + b*v + d         # u-component of the gradient
    gv = b*u + 2*c*v + e         # v-component of the gradient
    return np.asarray((gu, gv))  # converts input to an array
# **************************************************
    
# define parameter values
args = (2, 3, 7, 8, 9, 10)
# define array of initial guess
x0 = np.asarray((0, 0))

# call scipy minimization function using a nonlinear conjugate gradient algorithm
res1 = optimize.fmin_cg(f, x0, fprime=gradf, args=args)
# returns an array of the (u,v) coordinates of the function global minimium
print 'res1 = ', res1

