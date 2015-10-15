"""
Created on Fri Jul 26 14:36:45 2013
PHYS 510, Final Exam Problem 1
"""

# Final Exam Problem 1: Fukushima Radiation Leak (Version 1)
"""
The Fukushima data file contains the observational radiation leak data for the 
Fukushima Power Plant after the earthquake in 2011. The first column is the 
observed time every 10 minutes after the accident started and the second column 
is radiation measurement in micro Gray. Based on your knowledge from physics of 
radiation, find a model that fits the data the best.  
"""
# NOTE: Tries both exponential and logarithmic fit models (not very good fits...)


# NOTE: below fitting algorithm is a modified version of the Gauss-Newton Least Squares
#       method written for Assignment 5: code is from 'PHYS510 A5 Least Squares.py'. 
#       All comments about the code are in the original program named above. 


from math import exp,log
import numpy as np
import matplotlib.pyplot as plt


# import sample data set "FukushimaDataFile.txt"
t, rad = np.loadtxt('FukushimaDataFile.txt', skiprows=1, unpack=True)


# function for exponential model and derivatives for fitting
#*******************************************************************
def fitfuncsExp(c0,c1,x):
    f = c0*exp(c1*x)
    df1 = exp(c1*x)
    df2 = c0*x*exp(c1*x)
    
    return f,df1,df2
#*******************************************************************

    
# function for logarithmic model and derivatives for fitting
#*******************************************************************    
def fitfuncsLog(d0,d1,x):
    g = d0*log(x) + d1
    dg1 = log(x)
    dg2 = 1.0    
    
    return g,dg1,dg2
#*******************************************************************


# Gauss-Newton fitting method 
#*******************************************************************
def GaussNewtonFit(x,y,c0,c1,d0,d1,n):
    # solves for coefficients using a matrix solve of modified 'normal equations'
    # c0, c1 are initial guess of parameters in exponential y=c0*exp(c1*t)
    # d0, d1 are initial guess of parameters in logarithmic y=d0*ln(t)+d1
        
    pts = len(x)                                    
    expval = []  # stores exponential fit values
    logval = []  # stores logarithmic fit values
    
    # loop over specified number of iterations to update coefficients c0 and c1
    for k in range(n):
        # for exponential fit
        matrixAf = np.zeros((pts,2))
        matrixBf = np.zeros((pts,1))
        
        # for logarithmic fit
        matrixAg = np.zeros((pts,2))
        matrixBg = np.zeros((pts,1))
   
       # loop to populate arrays
        for i in range(pts):
            f,df1,df2 = fitfuncsExp(c0,c1,x[i])
            g,dg1,dg2 = fitfuncsLog(d0,d1,x[i])
            
            # for exponential fit
            matrixBf[i][0] = y[i]-f                     
            matrixAf[i][0] = df1                        
            matrixAf[i][1] = df2

            # for logarithmic fit
            matrixBg[i][0] = y[i]-g                    
            matrixAg[i][0] = dg1                        
            matrixAg[i][1] = dg2              
        
        # create normal equation (A^T)Ax=(A^T)b 
        # for exponential fit
        fAt = np.transpose(matrixAf)
        fAtA = np.dot(fAt,matrixAf)
        fAtB = np.dot(fAt,matrixBf)
        
        # for logarithmic fit
        gAt = np.transpose(matrixAg)
        gAtA = np.dot(gAt,matrixAg)
        gAtB = np.dot(gAt,matrixBg)
        
        coeffExp = np.linalg.solve(fAtA,fAtB)
        coeffLog = np.linalg.solve(gAtA,gAtB)
                      
        # exponential fit coefficients
        c0 = c0 + coeffExp[0]                           
        c1 = c1 + coeffExp[1]
        # logarithmic fit coefficients
        d0 = d0 + coeffLog[0]
        d1 = d1 + coeffLog[0]                           
        
    
    # we now have the coefficients for the fit
    # next plug x values back into fit model and compute new y fit values
    varsumExp = 0     
    varsumLog = 0
    # loop to get y fit values and compute variance on fit                           
    for i in range(pts):
        Expval,df1,df2 = fitfuncsExp(c0,c1,x[i])
        Logval,dg1,dg2 = fitfuncsLog(d0,d1,x[i])
        
        expval.append(Expval)
        logval.append(Logval)
                                
        # compute variance from erorr (error=actualY-fitY)
        errorExp = y[i] - Expval
        errorLog = y[i] - Logval
        
        varsumExp = varsumExp + errorExp**2
        varsumLog = varsumLog + errorLog**2
    
    # using variance as measure for 'goodness of fit'
    # computed as variance = (sum of squares of errors) / (nDataPts-nFitCoefficients-1)
    varianceExp = varsumExp/(pts-2-1)
    varianceLog = varsumLog/(pts-2-1)
    
    print c0,c1,d0,d1
    print varianceExp,varianceLog
            
    # plot data set with exponential fit and annotate variance
    plt.title('Fukushima Power Plant Radiation Leak (2011 Earthquake)\n Gauss-Newton Least Squares Fit')
    plt.xlabel('Time (min)')
    plt.ylabel('Absorbed Dose (micro Gray)')
    plt.plot(x,y,'b', linewidth=1.8, label='Observed Data')
    
    plt.plot(x,expval, 'r', linewidth=1.5, label='Exponential Model')
    plt.annotate('$variance=1.12$',fontsize=13,xy=(0.73,0.66),xycoords='figure fraction')
    plt.annotate('$f = 18.12\exp(-0.00056t)$',fontsize=14,xy=(0.40,0.66),xycoords='figure fraction')
    
    plt.plot(x,logval, 'g', linewidth=1.5, label='Logarithmic Model')
    plt.annotate('$variance=0.63$',fontsize=13,xy=(0.73,0.58),xycoords='figure fraction')
    plt.annotate('$f = -4.59\ln(t)+40.41$',fontsize=14,xy=(0.40,0.58),xycoords='figure fraction')
    
    plt.legend()
#******************************************************************* 
    
    
GaussNewtonFit(t,rad,24.0,-0.001,-5.0,40.0,5)