"""
Created on Thu Aug 29 19:14:15 2013
PHYS 613, Assignment 1
Nick Crump
"""

# Root Finding
"""
Design a program that finds the roots of a given function so the user can enter 
root brackets between which a single root is sought for the following methods:
    * Bisection method
    * Hybrid Bisection / Newton-Raphson method
"""

import math
import sympy


# enter root finding algorithm by Bisection method
#*********************************************************************** 
def rootBisection(g, f, xI, xF, Tol, nMax, Mthd):
        # initialize variables    
        error = 1
        n = 1
        xiMid = 0  # initial midpoint value to store the n-1 value
        
        # loop until error is less than input tolerance
        while error > Tol:
            xMid = 0.5*(xI+xF)
            
            # set up main Bisection method:
            # make bracket interval smaller each iteration until root is found
            # check conditions and update bracket points     
            if f(xI)*f(xMid) > 0:
                xI = xMid
                error = abs(xMid - xiMid)  # calculate approx error
                n = n + 1
                xiMid = xMid               # store the n-1 midpoint
                
            elif f(xI)*f(xMid) < 0:
                xF = xMid
                error = abs(xMid - xiMid)  # calculate approx error
                n = n + 1
                xiMid = xMid               # store the n-1 midpoint
        
        # output results to user         
        print 'Root =', xMid
        print 'Iterations = ', n-1
            
# end rootBisection function
#***********************************************************************


# enter root finding algorithm by hybrid Bisection/Newton-Raphson method
#*********************************************************************** 
def rootNewtonRaphson(g, f, xI, xF, Tol, nMax, Mthd):
        # initialize variables    
        error = 1
        n = 1
        
        # get symbolic derivative of input function      
        x = sympy.Symbol('x')        # define x as symbolic variable
        sdf = sympy.diff(g, x)       # get symbolic derivative of input function g
        df = sympy.lambdify(x, sdf)  # turn symbolic derivative into numeric function
        
        # check condition for starting x value
        if f(xI) < f(xF): xn0 = xI
        else: xn0 = xF
        
        # loop until error is less than input tolerance
        while error > Tol:
            
            # set up main Newton-Raphson method:
            # use derivative of function as tangent line to get new x point 
            # calcuate new x value from f(x) and df(x)
            xn1 = xn0 - (f(xn0)/df(xn0))
            
            # if new x point is outside bracket interval then use midpoint instead
            if xn1 < xI or xn1 > xF:
                xn1 = 0.5*(xI+xF)
                
                # check conditions and update bracket points     
                if f(xI)*f(xn1) > 0:
                    xI = xn1
                    
                elif f(xI)*f(xn1) < 0:
                    xF = xn1
                
            error = abs(xn1 - xn0)  # calculate approx error
            n = n + 1
            xn0 = xn1               # store the n-1 value for x
            
        # output results to user
        print 'Root =', xn1
        print 'Iterations =', n-1
            
# end rootNewtonRaphson function
#***********************************************************************  


# main program that calls functions above
#***********************************************************************
#--------------------------------------------------------------------------------------
# g = input function of x inside quotes (Ex: 'cos(x)-x')
# xI = initial x bracket value on one side of the root (entered as decimal)
# xF = initial x bracket value on the other side of the root (entered as decimal)
# Tol = the degree of accuracy in which to compute the root (Ex: 10e-5)
# nMax = max number of iterations to compute if Tolerance has not been met
# Mthd = choice of root finding method to use, entered inside quotes
       # (options are: 'bisection' or 'newton')
#-------------------------------------------------------------------------------------- 
def rootFind(g, xI, xF, Tol, nMax, Mthd):
    f = lambda x: eval(g)  # create function from input
    
    # check input bracket condition to ensure root is between braket points
    # check to make sure function is defined at bracket points too
    # this checks for examples like ln(0) = -inf or ln(0)*ln(4) = nan
    notnum1 = math.isnan(f(xI)*f(xF))
    notnum2 = math.isinf(f(xI)*f(xF))
    if f(xI)*f(xF) >= 0 or notnum1 == True or notnum2 == True: 
        print '\nNO ROOTS BETWEEN THESE POINTS OR FUNCTION IS UNDEFINED.'
        print 'PLEASE ADJUST BRACKET VALUES.\n'
        return
    
    if Mthd == 'bisection':
        # enter root finding algorithm by Bisection method
        rootBisection(g, f, xI, xF, Tol, nMax, Mthd)

    elif Mthd == 'newton':
       # enter root finding algorithm by hybrid Bisection/Newton-Raphson method
       rootNewtonRaphson(g, f, xI, xF, Tol, nMax, Mthd)

    else:
        print '\nINVALID ENTRY. PLEASE ADJUST INPUT AND TRY AGAIN.\n'
#***********************************************************************


#rootFind('x**3 - 2*x*cos(x) - sin(x) - 3', -1.0, 2.0, 10e-5, 30, 'bisection')
#rootFind('12*(x**4) - 44*(x**3) + 8*(x**2) + 29*x - 10', -2.0, 0, 10e-5, 30, 'bisection')
#rootFind('x**3 - 169', 3.0, 6.0, 10e-5, 30, 'newton')
