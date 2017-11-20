# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 10:15:59 2017

@author: Dominika
"""

#Write a general routine to perform numerical integration on a user-supplied 1D function using the extended
#trapezoidal rule. Design your program so that it automatically refines the result to a user-specified relative accuracy
#Here ‘user-supplied’ means that you should design your integration function so that one of its arguments
#is a pointer to another function, which will be the integrand. This other function should take a single argument
#(the value of the integration variable) and return the value of the integrand at that value of the integration variabl

import numpy as np

def phase_func(x,k):
    # Define phase of wavepacket at x for a wavenumber k
    return complex(0,k*x)

def wavefunc(x):
    # Define particle wavefunction as a function of x
    return (2./np.sqrt(np.pi))*np.exp(-0.5*x**2)*np.exp(phase_func(x,k))

def prob(x):
    # Calculates integrand at x
    return abs(wavefunc(x))**2

def prob_normalised(x):
    # Calculates normalised integrand at x
    return abs((1./c1)*wavefunc(x))**2

def calculate_step_static(x0,xend,n):
    # Calculates fixed width step h
    return (xend-x0)/n

#-------------------TRAPEZOIDAL RULE---------------------------

def trapezoidal_step(f,x_current,h):
    # Single step of trapezoidal rule 
    return f(x_current)*h

#--------------------SIMPSONS RULE-------------------------------

def simpson_step(f,x_current,h):
    return (4./3.)*trapezoidal_step(f,x_current+h,h) - (1./3.)*trapezoidal_step(f,x_current,h)

#-----------------GENERAL INTEGRATION METHOD--------------------------

def calculate_integral(f,x,x1,x2,method):
    # Integrate between x1 <= x <= x2 using method = trapezoidal/simpson
    # Calculate static stepsize h:
    h = calculate_step_static(-xend,xend,n)
    
    # Decide integration method:
    if method == 's':
        # Sum of f(x) evaluated inside integration range using simpson
        pts_mid = sum([(4./3.)*trapezoidal_step(f,xi+h,h) - (1./3.)*trapezoidal_step(f,xi,h) for xi in x if xi > x1 and xi < x2])
    
        # Sum of f(x) evaluated at end points of integration range
        pts_end = sum([0.5*trapezoidal_step(f,xi,h) for xi in x if xi == x1 or xi == x2])
        
    if method == 't':
        # Sum of f(x) evaluated inside integration range using trapezoidal
        pts_mid = sum([trapezoidal_step(f,xi,h) for xi in x if xi > x1 and xi < x2])
    
        # Sum of f(x) evaluated at end points of integration range
        pts_end = sum([0.5*trapezoidal_step(f,xi,h) for xi in x if xi == x1 or xi == x2])
    
    return (pts_mid + pts_end)


#---------------------SOME TEST PARAMETERS-----------------------

# Full range of x-coordinates
n = 10000
xend = 5
x = np.linspace(-xend,xend,n)

# Calculate normalisation constant such that area under probability curve == 1
c1 = (calculate_integral(prob,x,x[0],x[-1],'s'))**0.5

# Desired position range x1 --> x2 to integrate between
x1 = -5
x2 = 5

# Wavepacket wavenumber ?
k = 15

#Calculate normalised probability
p1 = calculate_integral(prob_normalised,x,x1,x2,'s')

import matplotlib.pyplot as plt

# Plot wavefunction
plt.figure()
plt.plot(x,[(1./c1)*wavefunc(x[i]) for i in range(len(x))])

# Plot probability distribution
plt.figure()
plt.plot(x,[(1./c1**2)*abs(wavefunc(x[i]))**2 for i in range(len(x))])
#plt.plot(x,[abs(pdf(x[i],k))**2 for i in range(len(x))],'k--')



    