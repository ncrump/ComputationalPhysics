"""
Created on Fri Jul 26 14:36:45 2013
PHYS 510, Final Exam Problem 2
"""

# Final Exam Problem 2: Fourier Signal Analysis
"""
Use Fourier power spectrum analysis to determine the frequency content
of the data in the Fourier data file. The first column in the data file
corresponds to time in seconds and the second column is the amplitude 
of the signal. Buried in the noise of the data is a signal of the form

      f(t) = sin(2pi*f1*t) + sin(2pi*f2*t) + sin(2pi*f3*t) + ...
      
Determine the frequencies present in the data and reconstruct the signal.
1) Determine the frequencies present in the signal.
2) Reconstruct the original signal using the frequencies found.
3) Add an extra freq (sin term) to the reconstructed signal with 
   fNew = fMax*1.05, where fMax is the max freq component found.
4) Determine the sampling rate required to clearly distinguish 
   the highest two frequencies. Show this by performing another
   power spectrum analysis and plot the two resolvable freqs. 
"""


import numpy as np
from math import pi
import matplotlib.pyplot as plt


# import sample data set "FourierData.txt"
t, amp = np.loadtxt('FourierData.txt', skiprows=1, unpack=True)


# Part 1: FFT function
#------------------------------------------------------------------
# Function that converts time domain data to freq domain using FFT
#*******************************************************************
def FFTspectrum(t,amp):
    # do FFT on time domain signal to get freq domain signal
    # t = time data in seconds, amp = amplitude values of signal
    
    pts = len(t)                # length of signal
    s = abs(np.fft.fft(amp))    # full two-sided spectrum, take abs() since FFT gives complex values
    spectrum = s[range(pts/2)]  # take only one side of spectrum
    
    # get min, max, number of points in freq spectrum
    fmin = min(spectrum)
    fmax = max(spectrum)
    npts = len(spectrum)
    
    # construct frequency data
    freq = np.linspace(fmin,fmax,npts)
    
    # plot time domain signal (amplitude vs time)
    fig = plt.subplot(211)
    plt.title('FFT of Time Signal to Freq Signal')
    plt.xlabel('Time (sec)')
    plt.ylabel('Amplitude')
    plt.plot(t,amp,'b.')
    
    # plot freq domain signal (freq spectrum vs freq)
    fig = plt.subplot(212)
    plt.xlabel('Frequency (kHz)')
    plt.ylabel('Spectrum')
    plt.plot(freq/1000,spectrum, linewidth=1.4)
    plt.xlim(xmin=-5.0)
    
    # make annotations on plots to identify freq components found
    plt.annotate('f1 = 0.10 Hz', fontsize=12, xy=(0.32,0.35), xycoords='figure fraction')
    plt.annotate('f2 = 2400.27 Hz', fontsize=12, xy=(0.32,0.30), xycoords='figure fraction')
    plt.annotate('f3 = 11200.89 Hz', fontsize=12, xy=(0.32,0.25), xycoords='figure fraction')
    plt.annotate('f4 = 18801.43 Hz', fontsize=12, xy=(0.32,0.20), xycoords='figure fraction')
    plt.show(fig)
    
    # print freq components observed in plot to get exact value
    print [round(freq[i],2) for i in range(npts) if spectrum[i] > 15000]
    
    
FFTspectrum(t,amp)
#*******************************************************************


# Parts 2 & 3: Reconstruct original signal from FFT
#------------------------------------------------------------------
# reconstruct origninal signal that was buried in the noise
# frequency components of signal from FFT of time domain signal
# scaling frequency and time data to be in kHz so we see signal and not noise
f1 = 0.10/1000      # kHz
f2 = 2400.27/1000   # kHz
f3 = 11200.89/1000  # kHz
f4 = 18801.43/1000  # kHz
fExtra = f4*1.05    # kHz (add extra freq to be 1.05*maxFreq for second plot)

# generate components of origninal signal f(t) = sin(2pi*fn*t) + ...
sig1 = np.sin(2*pi*f1*t)
sig2 = np.sin(2*pi*f2*t)
sig3 = np.sin(2*pi*f3*t)
sig4 = np.sin(2*pi*f4*t)
sigExtra = np.sin(2*pi*fExtra*t)  # add extra freq term for fExtra plot

# build reconstructed signal
newSig = sig1 + sig2 + sig3 + sig4  # original signal
ExtraSig = newSig + sigExtra        # original signal + fExtra term

# plot original signal on msec scale (t/1000) to see strong signal
fig = plt.subplot(211)
plt.title('Reconstructed Time-Domain Signals from FFT\n $f(t)= sin(2 \pi f_1 t) + sin(2 \pi f_2 t) + sin(2 \pi f_3 t) + sin(2 \pi f_4 t)$')
plt.ylabel('Amplitude')
plt.plot(t/1000,newSig, linewidth=1.2)

# plot original signal + fExtra term on msec scale (t/1000) to see strong signal
fig = plt.subplot(212)
plt.title('$g(t)= f(t) + sin(2 \pi 1.05f_{max} t)$')
plt.xlabel('Time (msec)')
plt.ylabel('Amplitude')
plt.plot(t/1000,ExtraSig, linewidth=1.2)
plt.show(fig)
#------------------------------------------------------------------


# Part 4: Determine required sampling rate to resolve highest freqs
#------------------------------------------------------------------
# do another FFT, now on the reconstructed data to resolve freq components
# mostly a copy/paste of the FFT function above with some parts modified

# do FFT on time domain reconstructed signal to get freq domain signal
# t = time data in seconds, amp = amplitude values of signal
# reconstructing signal again here but with freqs placed back in units of Hz to do FFT
f1Hz = 0.10           # Hz
f2Hz = 2400.27        # Hz
f3Hz = 11200.89       # Hz
f4Hz = 18801.43       # Hz
fExtraHz = f4Hz*1.05  # Hz (add extra freq to be 1.05*maxFreq for second plot)

# generate components of origninal signal again (but in Hz)
sig1Hz = np.sin(2*pi*f1Hz*t)
sig2Hz = np.sin(2*pi*f2Hz*t)
sig3Hz = np.sin(2*pi*f3Hz*t)
sig4Hz = np.sin(2*pi*f4Hz*t)
sigExtraHz = np.sin(2*pi*fExtraHz*t)  # add extra freq term for fExtra plot

# build reconstructed signal (but in Hz)
newSigHz = sig1Hz + sig2Hz + sig3Hz + sig4Hz  # original signal
ExtraSigaHz = newSigHz + sigExtraHz           # original signal + fExtra term
amp = ExtraSigaHz    

# now do FFT again, but on reconstructed signal (in Hz)
pts = len(t)                # length of signal
s = abs(np.fft.fft(amp))    # full two-sided spectrum, take abs() since FFT gives complex values
spectrum = s[range(pts/2)]  # take only one side of spectrum

# get min, max, number of points in freq spectrum
fmin = min(spectrum)
fmax = max(spectrum)
npts = len(spectrum)

# construct frequency data
freq = np.linspace(fmin,fmax,npts)

# plot time domain signal (amplitude vs time) of the kHz data to see strong signal
fig = plt.subplot(211)
plt.title('FFT of Reconstructed Time Signal $g(t)$ to Freq Signal')
plt.xlabel('Time (msec)')
plt.ylabel('Amplitude')
plt.plot(t/1000,ExtraSig,'r', linewidth=1.2)

# plot freq domain signal (freq spectrum vs freq) of the FFT Hz data
fig = plt.subplot(212)
plt.xlabel('Frequency (kHz)')
plt.ylabel('Spectrum')
plt.plot(freq/1000,spectrum, 'r', linewidth=1.4)
plt.xlim(xmin=-5.0)

# make annotations on plots to identify freq components found
plt.annotate('f1 = 0.86 Hz', fontsize=12, xy=(0.50,0.35), xycoords='figure fraction')
plt.annotate('f2 = 5281.36 Hz', fontsize=12, xy=(0.50,0.30), xycoords='figure fraction')
plt.annotate('f3 = 24640.26 Hz', fontsize=12, xy=(0.50,0.25), xycoords='figure fraction')
plt.annotate('f4 = 41359.46 Hz', fontsize=12, xy=(0.50,0.20), xycoords='figure fraction')
plt.annotate('f5 = 43427.22 Hz', fontsize=12, xy=(0.50,0.15), xycoords='figure fraction')
plt.show(fig)

# print freq components observed in plot to get exact value
print [round(freq[i],2) for i in range(npts) if spectrum[i] > 60000]
#------------------------------------------------------------------