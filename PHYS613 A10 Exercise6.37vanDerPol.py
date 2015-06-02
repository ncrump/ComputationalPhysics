"""
Created on Sat Nov 09 13:53:23 2013
PHYS 613, Assignment 10 
Nick Crump
"""

# Exercise 6.37
"""
Solve the van der Pol oscillator for parameters that exhibit chaotic
behavior and analyze the results in the frequency spectrum.
"""

import numpy as np
import ODEsolve as ode
import FFTspectrum as fft
import matplotlib.pyplot as plt


# function for van der Pol oscillator with driving force
# ********************************************************
def f(ti,IC): 
    eps = 3.0
    
    omg = 1.788                 # driving freq
    f0 = 5.0                    # driving amplitude
    dForce = f0*np.cos(omg*ti)  # driving force
    
    ui = IC[0]  # x'(t)
    xi = IC[1]  # x (t)
    
    func = np.array([-xi-eps*ui*(xi**2 - 1) + dForce, ui])
    return func
# ********************************************************


# first solve ODE
# -----------------------------
# set initial conditions as x'(0), x(0)
IC = [0, 0.5]
t0 = 0
tf = 30*2**np.pi/1.788
h = 0.002
tol = 1e-5

# call RKF45 method from my ODE solve module
t,s = ode.RKF45HO(f,IC,t0,tf,h,tol,stp=0)
pos = s[1]
vel = s[0]


# next do FFT on ODE solution
# -----------------------------
# define freq to sample signal
sfreq = 1.0/h    # Hz

# call FFT from my module to get freq spectrum
fs,spectrum = fft.FFTspectrum(t,pos,sfreq)
# convert freq axis in Hz to omega (omega = 2*pi*f)
fs = fs*2*np.pi
# find freq at strongest peak
indx = np.where(spectrum==max(spectrum))[0]
mx = round(fs[indx],3)
print 'peak = ',mx

# plot results
# -----------------------------
# plot phase space trajectory of chaotic oscillator
plt.figure(1)
plt.plot(pos,vel,'g', label='van der Pol Chaos')
plt.xlabel('Position')
plt.ylabel('Velocity')
plt.legend(loc=2)

# plot freq domain signal
plt.figure(2)
plt.plot(fs,spectrum, label='van der Pol Spectrum')
plt.xlabel('$\omega$',fontsize=18)
plt.ylabel('Normalized Magnitude')
plt.annotate('$peak=1.782$',fontsize=18,xy=(0.63,0.75),xycoords='figure fraction')
plt.xlim(0,10)
plt.legend()