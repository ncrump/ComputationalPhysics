"""
Created on Wed Sep 25 23:14:56 2013
PHYS 613, Assignment 5
Nick Crump
"""

# Problem 2 (EBB5) 
"""
Implement random number generating algorithms to obtain uniformly distributed and
normally distributed random number sequences and plot their distributions. 
"""

import RandomNumber as RN
import numpy as np
import matplotlib.pyplot as plt


# generate 10**6 random numbers using the Lehmer modulo generator
# make a histogram to check it is a uniform distribution
#*******************************************************************
# initial seed to start random number generator
seed = 112183.0

# call randUniform method from my random number module
randU = RN.randUniform(seed,10**6)

# generate random numbers using numpy for comparison
numpU = np.random.uniform(0,1,10**6)

# plot histogram of randUniform distribution
plt.figure()
plt.hist(randU,bins=20,label='Lehmer Modulo Generator')
plt.xlabel('Random Number')
plt.ylabel('Frequency')
plt.legend()

# plot histogram of numpy uniform distribution
plt.figure()
plt.hist(numpU,bins=20,label='Numpy Random Generator')
plt.xlabel('Random Number')
plt.ylabel('Frequency')
plt.legend()
#*******************************************************************


# implement Polar method for 10**6 random normal numbers
# implement Box-Muller method for 10**6 random normal numbers
# make a histogram to check it is a normal distribution
#*******************************************************************
# call randNormal1 Polar method from my random number module
randN1 = RN.randNormal1(10**6)

# call randNormal2 Box-Muller method from my random number module
randN2 = RN.randNormal2(10**6)

# generate random numbers using numpy for comparison
numpN = np.random.normal(0,1,10**6)

# plot histogram of randNormal1 (Polar) distribution
plt.figure()
plt.hist(randN1,bins=20,label='Polar Method')
plt.xlabel('Random Number')
plt.ylabel('Frequency')
plt.legend()

# plot histogram of randNormalal2 (Box-Muller) distribution
plt.figure()
plt.hist(randN2,bins=20,label='Box-Muller Method')
plt.xlabel('Random Number')
plt.ylabel('Frequency')
plt.legend()

# plot histogram of numpy normal distribution
plt.figure()
plt.hist(numpN,bins=20,label='Numpy Random Generator')
plt.xlabel('Random Number')
plt.ylabel('Frequency')
plt.legend()
#*******************************************************************