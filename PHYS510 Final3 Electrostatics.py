"""
Created on Fri Jul 26 14:36:45 2013
PHYS 510, Final Exam Problem 3
"""

# Final Exam Problem 3: Electrostatics
"""
Consider a conducting sphere with azimuthally symmetric potential across
the surface in spherical coordinates. The potential is symmetric in azimuth
around the sphere but not in elevation angle off the z-axis. The equation for
the potential V(r,theta) inside the sphere can be written as an infinite series
of Legengre Polynomials with coefficients An. The coefficients can be solved by 
integration to get the potential V(r,theta) inside the spehere. 

1) Use the recurrence relation to write a subroutine for generating
   higher order Legendre Polynomials. 
2) Pick your favorite numerical integrator and write a routine that
   can calculate the coefficients An.
3) Determine the number of coefficients needed so that the boundary
   equation given for the potential over the surface is satisfied to 
   single precision.
4) Use the coefficients obtained to plot the electric potential inside
   the sphere V(r,theta) as a function of r and theta. 
"""


import numpy as np
from matplotlib import cm
from math import pi,exp,sin,cos
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # used for 3D plotting


# Part 1: Generate Legendre Polynomials
#------------------------------------------------------------------
# function that generates higher order Legendre Polynomials P(n,x)
#*******************************************************************
def LegendrePoly(n,x):
    # n = order of polynomial, x = value at which to evaluate P(n,x)
    # uses recurrence relation to generate higher Legendre Polynomials
    
    P = 1.0                  # initial P(0,x) = P(n-1,x) to start recurrence relation
    Pn = x                   # initial P(1,x) = P(n,x) to start recurrence relation
    
    if n == 0: return P      # catch the n=0 case: P(0,x) = 1
    if n == 1: return Pn     # catch the n=1 case: P(1,x) = x
    
    else: 
        nRange = range(1,n)    # range of polynomial orders to iterate over
        
        # iterate recurrence relation up to order n-1 to get P(n,x)
        for i in nRange:
            Pnn = 2*x*Pn - P - (x*Pn - P)/(i+1)  # recurrence relation for P(n+1,x)
            
            P = Pn    # update P(n-1,x)
            Pn = Pnn  # update P(n,x)
            
        return Pnn
#*******************************************************************


# Part 2: Integrate to get coefficients An
#------------------------------------------------------------------

# Function to generate integrand to feed integrator to get An
#*******************************************************************
def Afunc(n,x):
    # n = order of function, x = value at which to evaluate function
    # NOTE: calls Legendre Polynomial function to get P(n,cos(x))
    
    r = 1.0                # radius of conducting sphere
    c = (2*n+1)/(2*r**n)   # constant outside integral
    v = exp(-(cos(x)**2))  # function for potential over surface of sphere
    
    f = c*v*sin(x)*LegendrePoly(n,cos(x))  # integrand function for An
    
    return f
#*******************************************************************
    

# Romberg integration method (integrator)
#*******************************************************************
def Romberg(a,b,n,N):
    # a, b = limits of integration (entered as decimals)
    # n = order of function An, N = step size up to 2**(N-1) points
    # NOTE: calls function 'Afunc(n,x)' above to get value of integrand
    # NOTE: using modified version of Romberg method written for Assignment 4
           
    Rnm = np.zeros((N,N))                  # create Romberg matrix
    
    fn0 = Afunc(n,a) + Afunc(n,b)          # get initial function value for matrix
    Rnm[0][0] = 0.5*(b-a)*fn0              # populate first element of matrix
        
    # this part populates the initial column of the Romberg matrix 
    # using modified Trapezoid rule
    # ---------------------------------------------------
    # iterate over interval step size up to N
    for s in range(1,N):              
        
        h = (b-a)/(2**s)                  # interval size
        pSum = 0                          # partial sum in integral routine
            
        # iterate over interval size up to 2**(N-1)
        for k in range(1, 2**s):
            if k % 2 != 0:                # only odds (k=1,3,5,..2**(N-1))
                fnx = Afunc(n, a + k*h)   # get function value for integration
                pSum = pSum + fnx         # partial sum in integral routine
              
        R = 0.5*Rnm[s-1][0] + h*pSum      # R(N,0) of Romberg method
        Rnm[s][0] = R                     # populate first column of matrix
    # ---------------------------------------------------
    
    # this part uses previous values in Romberg matrix to get new values
    # each new value is closer to the actual answer
    # --------------------------------------------------- 
    # iterate over rows and columns of Romberg matrix to get new values      
    for col in range(1,N):
        for row in range(col,N):
            romb1 = Rnm[row][col-1]
            romb0 = Rnm[row-1][col-1]
            coeff = 1/((4.0**col)-1)
                        
            # R(N,m) of Romberg method
            Rnm[row][col] = romb1 + coeff*(romb1 - romb0)
            
    An = Rnm[N-1][N-1]
    # ---------------------------------------------------
    
    return An  
#*******************************************************************


# Part 3: Determine number of terms needed for single precision 
#------------------------------------------------------------------
# NOTE: want accuracy = 7 decimal places for single precision
# NOTE: here we are approximating the electric potential over the sphere surface

# set initial parameters for integrations
n = 15                            # number of terms in series solution, also order of P(n,x)
a = 0                             # limit of integration
b = pi                            # limit of integration
r = 1.0                           # radius of conducting sphere
pts = 200                         # number of points to evaluate over function domain

# define arrays to store values
An = []                           # stores coefficients An
Vseries = []                      # stores series approx to surface potential
Vactual = []                      # stores actual values of surface potential (V actual)
Verrors = []                      # stores errors between approx potential and actual

theta = np.linspace(a,b,pts)      # points to evaluate potential functions

# loop over order n to get An coefficients for calculating potential V
for i in range(n):
    if i % 2 != 0: Ai = 0         # An=0 if n=odd due to potential function Vactual
    else: Ai = Romberg(a,b,i,10)  # otherwise calculate An if n=even
    An.append(Ai)

# loop over domain of function (0,pi) and order of coefficients (n) 
# this is to see what n is needed to approximate actual potential to single precision
# actual potential here is V=exp(-cos(x)**2)
for i in theta:
    V = 0                            # initiate potential                   
    
    for j in range(n):
        Pn = LegendrePoly(j,cos(i))  # get Legendre Polynomial P(n,x)
        Vpartial = An[j]*(r**j)*Pn   # calculate ith term in series potential
        V = V + Vpartial             # sum terms together in series
        
    Vact = exp(-cos(i)**2)           # calculate Vactual
    error = abs(Vact - V)            # get error between Vapprox and Vactual at each point
        
    Vseries.append(V)                # store approx potential values
    Vactual.append(Vact)             # store actual potential values
    Verrors.append(error)            # store errors at each point
    
# calculate RMS error between approximated and actual potentials over (0,pi)
squares = [i**2 for i in Verrors]
RMS = ((sum(squares))**0.5) / len(Verrors)
print 'n =','%02.0f' % n, '  RMS =','%0.10f' % RMS

# plot series approx to potential with actual potential for given order n
fig = plt.figure()
plt.title('Series Approximation to Known Potential Function')
plt.xlabel('Theta (deg) - Elevation Angle in Spherical Coords')
plt.ylabel('Potential over Spherical Surface')
plt.plot(theta*(180/pi),Vactual, linewidth=2.5, label='Known Potential')
plt.plot(theta*(180/pi),Vseries, marker='o', markerfacecolor='none',label='Order n='+str(n))
plt.ylim(ymax=1.1)
plt.legend()
plt.show()
#------------------------------------------------------------------


# Part 4: Use An's to plot the electric potential inside the sphere
#------------------------------------------------------------------
# NOTE: here we are generating the electric potential inside the sphere
# NOTE: this is the potential as a function of r and theta --> V(r,theta)    
    
radius = np.linspace(0,1,pts)                # define radius points to eval potential
theta = np.linspace(a,b,pts)                 # define theta points to eval potential

# to show a surface plot of potential, have to do 'meshgrid' to turn 1D arrays into 2D arrays
Rmesh, Tmesh = np.meshgrid(radius,theta)     # turns radius and theta arrays into 2D grid
Vinside = np.zeros([pts,pts])                # 2D array to store potential values over grid

# loop over domain of function V(r=0..1,theta=0..pi) to get potential inside sphere
for i in range(pts):                         # loops over rows of (r,theta) grid
    for j in range(pts):                     # loops over cols of (r,theta) grid
        V = 0                                # initiate potential
        Rij = Rmesh[i][j]                    # get radius value
        Tij = Tmesh[i][j]                    # get theta value
        
        # loop over coefficients (n)
        for k in range(n):
            Pn = LegendrePoly(k,cos(Tij))    # get Legendre Polynomial P(n,x)
            Vpartial = An[k]*(Rij**k)*Pn     # calculate term in series potential
            V = V + Vpartial                 # sum terms together in series
            
        Vinside[i][j] = V                    # store 2D array of potential values over grid
       
# generate 3D surface plot of electric potential inside sphere at a fixed azimuth
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.title('Electric Potential Inside Conducting Sphere')
ax.set_xlabel('Radius')
ax.set_ylabel('Theta (deg)')
ax.set_zlabel('Electric Potential')
ax.plot_surface(Rmesh,Tmesh*(180/pi),Vinside, cmap=cm.jet)
plt.show()

# generate contour plot of electric potential inside sphere at a fixed azimuth
fig = plt.figure()
plt.title('Electric Potential Inside Conducting Sphere')
plt.xlabel('Radius')
plt.ylabel('Theta (deg)')
plt.contourf(Rmesh,Tmesh*(180/pi),Vinside, cmap=cm.jet, levels=np.arange(0.30,1.05,0.05))
plt.colorbar()
plt.show()
#------------------------------------------------------------------