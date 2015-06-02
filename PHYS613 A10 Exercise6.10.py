"""
Created on Mon Nov 11 10:15:15 2013
PHYS 613, Assignment 10 
Nick Crump
"""

# Exercise 6.10
"""
Consider signals periodic in time and analyze sampling rates
using FFT to determine transformed freq spectrum of signal. 
"""

import numpy as np
import FFTspectrum as fft
import matplotlib.pyplot as plt


# define freq to sample signal
sfreq = 128  # Hz
Ttot = 5     # sec 
dt = 1.0/sfreq 

# define time points and signal values
t = np.arange(0,Ttot+dt,dt)
f = np.cos(6*np.pi*t)

# call FFT from my module to get freq spectrum
fs,spectrum = fft.FFTspectrum(t,f,sfreq)

# print max freq found in spectrum
mx = max(spectrum)
indx = np.where(spectrum==mx)[0]
print '\n','max freq =',round(fs[indx],2)

# plot time domain signal
plt.subplot(211)
plt.plot(t,f)
plt.xlabel('Time (sec)')
plt.ylabel('Amplitude')

# plot freq domain signal
plt.subplot(212)
plt.plot(fs,spectrum, label='$128\ Hz\ sample$' +'\n' '$\ \ \ \ 3\ Hz\ freq$')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Normalized Magnitude')
plt.xlim(0,20)
plt.legend()