"""
Created on Wed Jul 17 22:32:39 2013
PHYS 510, Final Project GUI Version
"""

# Simulating the Foucault Pendulum 
"""
This program simulates a Foucault Pendulum using the fifth-order Runge-Kutta-Fehlberg
method (RKF45) with adaptive step size to solve the equations of motion. The equations
of motion are solved in the Earth's rotating coordinate frame so the Coriolis affect is 
seen as a precession of the pendulum motion. The graphical output is displayed in an 
interactive GUI to allow the user to adjust parameters of the pendulum's motion. 
"""

# *****IMPORTANT SOLUTION NOTE*****
"""
It turns out that the ODEs of the Foucault Pendulum are STIFF ODEs that tend to be 
unstable when using certain numerical solutions unless the step size is taken to be 
extremely small. It was noticed here that even when using the RKF45 method with 
adaptive step size, the solution experiences numerical damping when solved over a 
long time interval. For this reason there are special numerical methods for these 
ODEs called STIFF SOLVERS which should work better.
"""


import numpy as np
from math import pi,sin,cos
import matplotlib.pyplot as plt
from matplotlib.widgets import Button,CheckButtons,RadioButtons,Slider

#-------------------------------------------------
# Foucault ODE (Coupled 2nd-Order): 
# x''(t) =  2A y'(t) - B x(t)
# y''(t) = -2A x'(t) - B y(t)
#-------------------------------------------------
# Becomes a system of 1st-Order ODEs:
# u'(t) =  2A v(t) - B x(t)
# v'(t) = -2A u(t) - B y(t)
# x'(t) = u(t)
# y'(t) = v(t)
#-------------------------------------------------
# A = z-component earth's angular velocity (omegaZ)
# B = natural freq of pendulum (omegaP)
# u(t) = x'(t), u'(t) = x''(t)
# v(t) = y'(t), v'(t) = y''(t)
#*************************************************************************************
def ODEfunc(S,omegaZ,omegaP):
    # array of functions for derivative of ODE (vector equation represented as array)
    
    ui = S[0]  # x'(t)
    vi = S[1]  # y'(t)
    xi = S[2]  # x (t)
    yi = S[3]  # y (t)
    
    f = np.array([2*omegaZ*vi - (omegaP**2)*xi, -2*omegaZ*ui - (omegaP**2)*yi, ui, vi])
    
    return f
#*************************************************************************************


# RKF45 method for Foucault Pendulum ODE with adaptive step size
#*************************************************************************************
def RKF45(t0,tf,x0,y0,dx0,dy0,omegaZ,omegaP,h,eps):
    #-----------------------------------------------------------
    # t0, tf = start and end time (sec)
    # x0, y0, dx0, dy0 = initial conditions (m, m, m/sec, m/sec)
    # omegaZ, omegaP = Earth and pendulum angular velocities (rad/sec)
    # h = initial step size
    # eps = max relative error for computing adaptive h
    #-----------------------------------------------------------  

    # create arrays for storing solution values
    t = [t0]
    x = [x0]
    y = [y0]
    
    # set initial values
    ti = t0
    Si = np.array([dx0, dy0, x0, y0])
            
    # RKF45 coefficients for intermediate functions, solution, and error functions
    # intermediate function f1
    a1, b1, = 0.25, 0.25
    # intermediate function f2
    a2, b2, c2 = 3.0/8.0, 3.0/32.0, 9.0/32.0
    # intermediate function f3
    a3, b3, c3, d3 = 12.0/13.0, 1932.0/2197.0, 7200.0/2197.0, 7296.0/2197.0
    # intermediate function f4
    a4, b4, c4, d4, e4 = 1, 439.0/216.0, 8, 3680.0/513.0, 845.0/4104.0
    # intermediate function f5
    a5, b5, c5, d5, e5, f5 = 0.5, 8.0/27.0, 2, 3544.0/2565.0, 1859.0/4104.0, 11.0/40.0
    
    # fifth-order solution function
    Y1, Y2, Y3, Y4, Y5 = 16.0/135.0, 6656.0/12825.0, 28561.0/56430.0, 9.0/50.0, 2.0/55.0
    # error function
    E1, E2, E3, E4, E5 = 1.0/360.0, 128.0/4275.0, 2197.0/75240.0, 1.0/50.0, 2.0/55.0
    
    # RKF45 method to solve ODE over each adaptive step h
    # NOTE: since we have a system of 4 equations, the RKF45 functions are now arrays
    # NOTE: intermediate functions f0,f1,f2,f3,f4,f5 are vectors represented as arrays
    # NOTE: solution and error functions (F and error) are also vector arrays
    while ti <= tf:
                
        # compute initial intermediate function f0
        f0 = ODEfunc(Si,omegaZ,omegaP)
        
        # compute dependent variable arrays and intermediate function arrays f1-f5
        SI = Si + b1*h*f0
        f1 = ODEfunc(SI,omegaZ,omegaP)
        
        SI = Si + b2*h*f0 + c2*h*f1
        f2 = ODEfunc(SI,omegaZ,omegaP)
        
        SI = Si + b3*h*f0 - c3*h*f1 + d3*h*f2
        f3 = ODEfunc(SI,omegaZ,omegaP)
        
        SI = Si + b4*h*f0 - c4*h*f1 + d4*h*f2 - e4*h*f3
        f4 = ODEfunc(SI,omegaZ,omegaP)
        
        SI = Si - b5*h*f0 + c5*h*f1 - d5*h*f2 + e5*h*f3 - f5*h*f4
        f5 = ODEfunc(SI,omegaZ,omegaP)
        
        # compute solution function array
        F = Si + h*(Y1*f0 + Y2*f2 + Y3*f3 - Y4*f4 + Y5*f5)  
        
        
        # compute error function array and determine RMS error
        error = h*(E1*f0 - E2*f2 - E3*f3 + E4*f4 + E5*f5)
        RMS = (sum(error**2) / len(F))**0.5
        
        # set adaptive step size based on RMS error vs specified max acceptable error
        hNew = h*(eps/RMS)**0.2
        
        # if error is small enough, accept solution and continue
        # otherwise iterate again with adjusted smaller step size
        if hNew >= h: 
            ti = ti + h     # update independent variable t
            Si = F          # update dependent variables [dx,dy,x,y]
            t.append(ti)    # store t values to array
            x.append(F[2])  # store x values of solution array to array        
            y.append(F[3])  # store y values of solution array to array
        
        h = 0.9*hNew        # try smaller step size this time
                
    return t,x,y
#*************************************************************************************


# Main program that creates GUI widget and plots solutions
#*************************************************************************************        
# set initial conditions and start parameters
#-------------------------------------------------------
t0, tf = 0, 120                       # start & end time (sec)
x0, y0 = 1.0, 0                       # initial x,y pos (m)
dx0, dy0 = 0, 0                       # initial x,y vel (m/sec)
h, eps = 0.1, 1.0e-5                  # initial steps size and desired accuracy
l0, lf = 1, 100                       # widget bounds on string length (m)
lat0, latf = -90, 90                  # widget bounds on latitude (deg)

# fixed parameters
omegaE = (2*pi)/(24*3600)             # earth angular vel (rad/sec)
g = 9.8                               # accel of gravity  (m/sec**2)

# adjustable initial parameters
ls = 67                               # string length (m)
lats = 48.85                          # latitude      (deg)

# resulting initial parameters
lat = lats*(pi/180)                   # latitude            (rad)
omegaP = (g/ls)**0.5                  # oscillator freq     (rad/sec)
omegaZ = omegaE*sin(lat)              # earth z-angular vel (rad/sec)
rate = omegaZ*(180/pi)*3600           # precession rate     (deg/hr)
pRate = '%0.2f' % rate
#-------------------------------------------------------


# call methods and setup plot window for numerical and analytic solutions
#-------------------------------------------------------
# call ODE solver method for RKF45 numerical solution
t,x,y = RKF45(t0,tf,x0,y0,dx0,dy0,omegaZ,omegaP,h,eps)

# generate analytic solution for overlay plotting with numeric solution
xx = [x0*cos(omegaP*i)*cos(omegaZ*i) for i in t]
yy = [-x0*cos(omegaP*i)*sin(omegaZ*i) for i in t]

# setup plot window
fig = plt.subplot(111)
plt.subplots_adjust(left=0.18, right=0.80, bottom=0.28)
plt.title('Oscillation of Foucault Pendulum')
plt.xlabel('x-East (m)')
plt.ylabel('y-North (m)')
plt.annotate('Precession Rate',xy=(0.82,0.88),xycoords='figure fraction')
plt.annotate('(deg/hr)',xy=(0.86,0.84),xycoords='figure fraction')
pa, = fig.plot(xx,yy,'r', visible=False)
pn, = fig.plot(x,y, 'b')

# setup precession rate output box and text
plt.text(1.085,0.82,'XXXXX',color='lightgrey',fontsize=16,transform=fig.transAxes,bbox=dict(facecolor='lightgrey'))
plt.text(1.1,0.82,pRate,fontsize=16,transform=fig.transAxes)
#-------------------------------------------------------


# create widget buttons and slider bars
#-------------------------------------------------------
# define slider and button positions
buttonChkLoc = plt.axes([0.835, 0.56, 0.13, 0.15])    # data plotting check box
buttonBoxLoc = plt.axes([0.85, 0.38, 0.10, 0.15])     # linestyle plotting check box
buttonResLoc = plt.axes([0.85, 0.30, 0.10, 0.05])     # reset button
sliderLenLoc = plt.axes([0.18, 0.14, 0.62, 0.025])    # string length slider
sliderLatLoc = plt.axes([0.18, 0.09, 0.62, 0.025])    # latitude slider
sliderTimeLoc = plt.axes([0.18, 0.04, 0.62, 0.025])   # time slider

# draw sliders and buttons
buttonChk = CheckButtons(buttonChkLoc, (' Analytic',' Numeric'), (False,True))
buttonBox = RadioButtons(buttonBoxLoc, (' line',' point',' circle'), active=0)
buttonRes = Button(buttonResLoc, 'Reset',color='w',hovercolor='lightgreen')
sliderLen = Slider(sliderLenLoc, 'Length (m)', l0, lf, valinit=ls)
sliderLat = Slider(sliderLatLoc, 'Latitude (deg)', lat0, latf, valinit=lats)
sliderTime = Slider(sliderTimeLoc, 'Time (sec)', t0, tf, valinit=tf)
#-------------------------------------------------------


# create functions to handle widget events from buttons and sliders
#-------------------------------------------------------
# function to update plots when SLIDER BARS change 
def updatefunc(val):
    uLen = sliderLen.val
    uLat = sliderLat.val
    uTime = sliderTime.val
    
    # resulting adjusted parameters from slider bar selection
    ulat = uLat*(pi/180)                   # rad
    uomegaP = (g/uLen)**0.5                # rad/sec
    uomegaZ = omegaE*sin(ulat)             # rad/sec
    
    # update precession rate and redraw box and text when latitude changes
    rate = round(uomegaZ*(180/pi)*3600,2)  # deg/hr
    pRate = '%0.2f' % rate
    plt.text(1.085,0.82,'XXXXX',color='lightgrey',fontsize=16,transform=fig.transAxes,bbox=dict(facecolor='lightgrey'))
    plt.text(1.1,0.82,pRate,fontsize=16,transform=fig.transAxes)
    
    # update solutions from function calls with adjusted parameters
    t,x,y = RKF45(t0,uTime,x0,y0,dx0,dy0,uomegaZ,uomegaP,h,eps)
    xx = [x0*cos(uomegaP*i)*cos(uomegaZ*i) for i in t]
    yy = [-x0*cos(uomegaP*i)*sin(uomegaZ*i) for i in t]
    
    # adjust plot with new solution values
    pn.set_xdata(x)
    pn.set_ydata(y)
    pa.set_xdata(xx)
    pa.set_ydata(yy)
    
sliderLen.on_changed(updatefunc)
sliderLat.on_changed(updatefunc)
sliderTime.on_changed(updatefunc)

# function to handle CHECK BUTTON for turning on/off analytic & numeric overlays
def checkfunc(label):
    if label == ' Analytic': pa.set_visible(not pa.get_visible())
    elif label == ' Numeric': pn.set_visible(not pn.get_visible())
    plt.draw()
buttonChk.on_clicked(checkfunc)

# function to handle plot marker RADIO BUTTONS 
def markerfunc(label):
    if label == ' line': 
        pn.set_linestyle('-')
        pn.set_marker('')
    elif label == ' point': 
        pn.set_linestyle('.')
        pn.set_marker('.')
    elif label == ' circle': 
        pn.set_linestyle('.')
        pn.set_marker('o')
    plt.draw()
buttonBox.on_clicked(markerfunc)

# function to handle RESET BUTTON 
def resetfunc(event):
    sliderLen.reset()
    sliderLat.reset()
    sliderTime.reset()
buttonRes.on_clicked(resetfunc)
#-------------------------------------------------------
#*************************************************************************************