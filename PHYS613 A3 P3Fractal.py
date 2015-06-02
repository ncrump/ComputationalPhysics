"""
Created on Mon Sep 16 18:31:35 2013
PHYS 613, Assignment 3
Nick Crump
"""

# Problem 3 (EBB3)
"""
Finding best fit parameters using a non-linear least squares method to 
data with a functional dependence of the fractal dimension. 
"""

from math import log,e
import numpy as np
import matplotlib.pyplot as plt


# import sample data set "P3FractalTable.txt"
c, fd = np.loadtxt('P3FractalTable.txt', skiprows=1, unpack=True)


# Function for fractal dimension f = a + b x**c
#*******************************************************************
def fitfunc(a,b,c,x):
    f = a + b*(x**c)
    dfa = 1
    dfb = x**c
    dfc = b*(x**c)*log(x,e)
    
    return f,dfa,dfb,dfc
#*******************************************************************


# Gauss-Newton fitting method 
#*******************************************************************
def GaussNewtonFit(x,y,a,b,c,n):
    # solves for coefficients using a matrix solve of modified 'normal equations'
    # uses fractal dimension model f = a + b x**c
    # a,b,c are the best fit model parameters to be determined
        
    pts = len(x)                                    
    fdVals = []  # stores fit values
    
    # loop over specified number of iterations to update coefficients a,b,c
    for k in range(n):

        matrixA = np.zeros((pts,3))
        matrixB = np.zeros((pts,1))
   
       # loop to populate arrays
        for i in range(pts):
            f,dfa,dfb,dfc = fitfunc(a,b,c,x[i])
            
            # populate matrix elements 
            matrixB[i][0] = y[i]-f                     
            matrixA[i][0] = dfa                        
            matrixA[i][1] = dfb  
            matrixA[i][2] = dfc
        
        # create normal equation (A^T)Ax=(A^T)b 
        At = np.transpose(matrixA)
        AtA = np.dot(At,matrixA)
        AtB = np.dot(At,matrixB)
        
        coeff = np.linalg.solve(AtA,AtB)
                      
        # get coefficients
        a = a + coeff[0]                           
        b = b + coeff[1]
        c = c + coeff[2]
        
    
    # we now have the coefficients for the fit
    # next plug x values back into fit model and compute new y fit values
    varsum = 0     
    
    # loop to get y fit values and compute variance on fit                           
    for i in range(pts):
        val,dfa,dfb,dfc = fitfunc(a,b,c,x[i])
        
        fdVals.append(val)
                                
        # compute variance from erorr (error=actualY-fitY)
        error = y[i] - val
        
        varsum = varsum + error**2
    
    # using variance as measure for 'goodness of fit'
    # computed as variance = (sum of squares of errors) / (nDataPts-nFitCoefficients-1)
    variance = varsum/(pts-3-1)
    
    print 'coeffs=',a,b,c
    print 'variance=',variance
            
    # plot data set with fractal dimension fit model
    plt.xlabel('Concentration')
    plt.ylabel('Fractal Dimension')
    plt.plot(x,y,'bo--',linewidth=1.8, label='Data Set')
    
    plt.plot(x,fdVals, 'ro-', linewidth=1.5, label='Fractal Model')
    plt.annotate('$d_f = 1.78 + 0.86x^{0.45}$',fontsize=15,xy=(0.16,0.70),xycoords='figure fraction')
    
    plt.legend(loc=2)
#******************************************************************* 
    
    
GaussNewtonFit(c,fd,1.8,1.0,0.6,5)