#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 18:21:41 2022

@author: johanna
"""

from scipy.stats import nbinom
import matplotlib.pyplot as plt
import numpy as np


n = 5
p = 0.01

nes = np.arange(1, 10, 1)
for n in nes : 
    fig, ax = plt.subplots(1,1)

    x = np.arange(nbinom.ppf(0.001, n, p),
                  nbinom.ppf(0.99, n, p))
    
    ax.plot(x, nbinom.pmf(x, n, p), 'bo', ms=8, label='nbinom pmf')
    ax.vlines(x, 0, nbinom.pmf(x,n,p), colors = 'b', lw=5, alpha=0.5)



