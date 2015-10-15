"""
Created on Fri Jul 26 14:36:45 2013
PHYS 510, Final Exam Problem 1
"""

# Final Exam Problem 1: Fukushima Radiation Leak (Version 2)
"""
The Fukushima data file contains the observational radiation leak data for the 
Fukushima Power Plant after the earthquake in 2011. The first column is the 
observed time every 10 minutes after the accident started and the second column 
is radiation measurement in micro Gray. Based on your knowledge from physics of 
radiation, find a model that fits the data the best.  
"""
# NOTE: Tries a double exponential model (works very well!)


# NOTE: below fitting algorithm is a modified version of the Gauss-Newton Least Squares
#       method written for Assignment 5: code is from 'PHYS510 A5 Least Squares.py'. 
#       All comments about the code are in the original program named above. 


from math import exp
import numpy as np
import matplotlib.pyplot as plt


# import sample data set "FukushimaDataFile.txt"
t, rad = np.loadtxt('FukushimaDataFile.txt', skiprows=1, unpack=True)


# Function for double exponential model and derivatives for fitting
#*******************************************************************
def fitfunc(a,b,c,d,x):
    f = a*exp(b*x) + c*exp(d*x)
    dfa = exp(b*x)
    dfb = a*x*exp(b*x)
    dfc = exp(d*x)
    dfd = c*x*exp(d*x)
    
    return f,dfa,dfb,dfc,dfd
#*******************************************************************


# Gauss-Newton fitting method 
#*******************************************************************
def GaussNewtonFit(x,y,a,b,c,d,n):
    # solves for coefficients using a matrix solve of modified 'normal equations'
    # uses double exponential model y=a*exp(bt) + c*exp(dt)
    # a,b,c,d are the best fit model parameters to be determined
        
    pts = len(x)                                    
    expval = []  # stores fit values
    
    # loop over specified number of iterations to update coefficients a,b,c,d
    for k in range(n):

        matrixA = np.zeros((pts,4))
        matrixB = np.zeros((pts,1))
   
       # loop to populate arrays
        for i in range(pts):
            f,dfa,dfb,dfc,dfd = fitfunc(a,b,c,d,x[i])
            
            # populate matrix elements 
            matrixB[i][0] = y[i]-f                     
            matrixA[i][0] = dfa                        
            matrixA[i][1] = dfb  
            matrixA[i][2] = dfc
            matrixA[i][3] = dfd
        
        # create normal equation (A^T)Ax=(A^T)b 
        At = np.transpose(matrixA)
        AtA = np.dot(At,matrixA)
        AtB = np.dot(At,matrixB)
        
        coeff = np.linalg.solve(AtA,AtB)
                      
        # get coefficients
        a = a + coeff[0]                           
        b = b + coeff[1]
        c = c + coeff[2]
        d = d + coeff[3]
        
    
    # we now have the coefficients for the fit
    # next plug x values back into fit model and compute new y fit values
    varsum = 0     
    
    # loop to get y fit values and compute variance on fit                           
    for i in range(pts):
        val,dfa,dfb,dfc,dfd = fitfunc(a,b,c,d,x[i])
        
        expval.append(val)
                                
        # compute variance from erorr (error=actualY-fitY)
        error = y[i] - val
        
        varsum = varsum + error**2
    
    # using variance as measure for 'goodness of fit'
    # computed as variance = (sum of squares of errors) / (nDataPts-nFitCoefficients-1)
    variance = varsum/(pts-4-1)
    
    print 'coeffs=',a,b,c,d
    print 'variance=',variance
            
    # plot data set with exponential fit and annotate variance
    plt.title('Fukushima Power Plant Radiation Leak (2011 Earthquake)\n Gauss-Newton Least Squares Fit')
    plt.xlabel('Time (min)')
    plt.ylabel('Absorbed Dose (micro Gray)')
    plt.plot(x,y,'b', linewidth=1.8, label='Observed Data')
    
    plt.plot(x,expval, 'r', linewidth=1.5, label='Exponential Model')
    plt.annotate('$f = 7.86\exp(-0.00027t) + 13.99\exp(-0.00153t)$',fontsize=14,xy=(0.40,0.68),xycoords='figure fraction')
    plt.annotate('$variance=0.24$',fontsize=13,xy=(0.65,0.60),xycoords='figure fraction')
    
    plt.legend()
#******************************************************************* 
    
    
GaussNewtonFit(t,rad,20.0,-0.0001,5.0,-0.001,5)