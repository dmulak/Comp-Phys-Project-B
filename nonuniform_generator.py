# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 11:20:08 2017

@author: Dominika
"""
import numpy as np


# Computes the inverse of the cumulative distribution function of pdf(x)
def inverse_cdf(y):
    ac = -0.48
    bc = 0.98
    return (1./ac)*(-bc + np.sqrt(bc**2 + 2*ac*y))
    
# Computes n pseudo-random numbers from a seed in 
def non_uniform_gen():
    nrand = np.random.uniform() # generate seed in range [0,1]
    a = 1664525; m = 2**32; c = 2013904223 # recommended parameters from Numerical Recipes
    y = ((a*nrand+c)%m)/float(m)
    x = inverse_cdf(y)
    return x