#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 18:21:41 2022

@author: johanna
"""
import os
import scipy
from scipy.stats import nbinom
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm


# check means in our real data
dataloc = os.path.expanduser("~/proj_LDHnet202204/data4py/")
emx = pd.read_csv(dataloc + "p4py_rawq.tsv", sep="\t",
                   index_col=0) # expression

means_obs = []
for kc in emx.columns:
    tmp = np.mean(emx[kc])
    means_obs.append(tmp)
    
    
plt.plot(means_obs, "o", c="darkgreen")
plt.title("observed mean of raw counts across subjects")
plt.show()

# as E(X) = k / p , save p and k for which the likelihood of observed distr is max
p_k_ = []

ares





# p (X = n ) = (n-1 k-1) (1-p )**n-k * p **k  is the counts distribution law
# E(X) = k / p
# V(X) = k(1-p) / p**2


initialN = 15
subjmu_ = [150, 200]

# means by subject from a uniform distribution
means_subs_sim =  np.random.uniform(low=subjmu_[0], high=subjmu_[1], size=initialN)


plt.plot(means_subs_sim, "o")
plt.title("synthetic mean expression across subjects, n="+str(initialN))
plt.show()

for k in range(len(means_sub_sim)):
    

# once the matrices are ready, calculate CPM, across columns : 
# store pickle dico :  with raw counts and cpm

################ ====
def initialworkaround():
    n = 5
    p = 0.01
    nes = np.arange(1,3, 1)
    for n in nes : 
        fig, ax = plt.subplots(1,1)
    
        x = np.arange(nbinom.ppf(0.001, n, p),
                      nbinom.ppf(0.99, n, p))
        
        ax.plot(x, nbinom.pmf(x, n, p), 'bo', ms=8, label='nbinom pmf')
        ax.vlines(x, 0, nbinom.pmf(x,n,p), colors = 'b', lw=5, alpha=0.5)
    
    n = 300
    he = nbinom.rvs(1,p, size = 100)
    fig, ax = plt.subplots(1,1)
    plt.plot(he)
