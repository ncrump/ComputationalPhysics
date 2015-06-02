"""
Created on Tue Nov 26 10:34:33 2013
PHYS 613, Project 4 
Nick Crump
"""

# Perturbed Harmonic Oscillator
"""
This solves the Hamiltonian equations of motion for the 
2-D perturbed (chaotic) harmonic oscillator. 
"""

import numpy as np
import ODEsolve as ode
import FFTspectrum as fft
import matplotlib.pyplot as plt


# function for Hamiltonian equations of motion 
# for 2-D perturbed harmonic oscilator 
# ********************************************************
def f(ti,IC): 
    C = 5.2        # coupling parameter (change for chaos)
    L = 1.0        # scaling parameter
    
    xi = IC[0]   # x(t)
    yi = IC[1]   # y(t)
    pxi = IC[2]  # px(t) 
    pyi = IC[3]  # py(t)
    
    x = pxi
    y = pyi
    px = -xi - 4*L*xi**3 - 4*L*C*xi*yi**2
    py = -yi - 4*L*yi**3 - 4*L*C*yi*xi**2
    
    func = np.array([x,y,px,py])
    return func
# ********************************************************


# solve ODE
# -----------------------------
# set initial conditions as x(0), y(0), px(0), py(0)
theta0 = 10*(np.pi/180)
x0 = np.sin(theta0)
y0 = np.cos(theta0)

IC = [x0, y0, 0, 0]
t0 = 0
tf = 200
h = 0.01
tol = 1e-5

# call RKF45 method from my ODE solve module
t,s = ode.RKF45HO(f,IC,t0,tf,h,tol,stp=0)
x = s[0]
y = s[1]
px = s[2]
py = s[3]
# -----------------------------

# get points to construct Poincare section
# -----------------------------
# get points where py = 0
indx = []
N = len(py)
for i in range(N-1):
    if py[i]*py[i+1] < 0:
        indx.append(i)

print '\n','points = ',len(indx)

# Poincare section will plot phase portrait at these points
yy = y[indx]
pxx = px[indx]
# -----------------------------

# take FFT of ODE solution
# -----------------------------
# define freq to sample signal
sfreq = 1.0/h    # Hz

# call FFT from my module to get freq spectrum
fsX,specX = fft.FFTspectrum(t,x,sfreq)
fsY,specY = fft.FFTspectrum(t,y,sfreq)
# convert freq axis in Hz to omega (omega = 2*pi*f)
fsX = fsX*2*np.pi
fsY = fsY*2*np.pi
# convert freq spectrum to power spectrum
psX = specX**2
psY = specY**2
# find freq at strongest peak
indX = np.where(psX == max(psX))[0]
indY = np.where(psY == max(psY))[0]
mX = round(fsX[indX],3)
mY = round(fsY[indY],3)
print 'x peak = ',mX
print 'y peak = ',mY

# plot results
# -----------------------------
# plot position vs time
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,x,'b', label='x')
plt.legend()
plt.subplot(2,1,2)
plt.plot(t,y,'r', label='y')
plt.xlabel('Time')
plt.legend()

# plot phase portrait - y vs x
plt.figure(2)
plt.plot(x,y,'b')
plt.xlabel('x')
plt.ylabel('y')

# plot phase trajectory - px vs y
plt.figure(3)
plt.plot(y,px,'b')
plt.xlabel('y-Position')
plt.ylabel('x-Momentum')

# plot Poincare section - px vs y at py = 0
plt.figure(4)
plt.plot(yy,pxx,'b.')
plt.xlabel('y-Position')
plt.ylabel('x-Momentum')

# plot power spectrum
plt.figure(5)
plt.plot(fsX,psX,'b', label='x')
plt.plot(fsY,psY,'r', label='y')
plt.xlabel('$\omega$',fontsize=18)
plt.ylabel('Power')
plt.xlim(0,4)
plt.legend()