# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 12:00:25 2017

@author: Dominika
"""

import numpy as np
import matplotlib.pyplot as plt

def wavefunc(x):
    # Define particle wavefunction as a function of x
    return (2./np.sqrt(np.pi))*np.exp(-0.5*x**2)*np.exp(15j*x)

def prob(x):
    # Calculates integrand at x
    return abs(wavefunc(x))**2

def prob_normalised(x):
    # Calculates normalised integrand at x
    return abs((1. / c1m)*wavefunc(x))**2    
    
def ymax_sampling(f, x, x1, x2):
    return 1.01*max([f(xi) for xi in x if xi >= x1 and xi <= x2])

def calculate_area(f, x1, x2):
    return (x2 - x1)*ymax_sampling(f, x, x1, x2)

def montecarlo(f, x1, x2, n = 1000):
    count = 0
    for i in range(n):
        # Generate n test coords (x,y)
        xtest = np.random.uniform(x1, x2)
        ytest = np.random.uniform(0., ymax_sampling(f, x, x1, x2))
        ytrue = f(xtest)
        if ytest <= ytrue:
            count += 1
    #print(count)
    return (count / n)*calculate_area(f, x1, x2)
        
    

#---------------------SOME TEST PARAMETERS-----------------------

# Wavepacket wavenumber ?
k = 15

# Full range of x-coordinates
nsteps = 1000
xend = 5
xstart = -xend
x = np.linspace(-xend,xend,nsteps)
#tol = 0.2

# Calculate normalisation constant such that area under probability curve == 1
c1m = (montecarlo(prob, xstart, xend))**0.5

# Desired position range x1 --> x2 to integrate between
x1 = 1
x2 = 2

# Calculate normalised probability
p1m = montecarlo(prob_normalised, x1, x2)
#p1m = montecarlo(prob_normalised, xstart, xend)

# Plot probability distribution
plt.figure()
plt.plot(x,[(1./c1m**2)*abs(wavefunc(x[i]))**2 for i in range(len(x))])
plt.plot(x, -0.48*abs(x) + 0.98)
plt.xlabel('$x$')
plt.ylabel('$P(x)$')