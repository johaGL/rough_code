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
import seaborn as sns

#import statsmodels.api as sm
from scipy.special import gammaln
from scipy.special import psi
from scipy.special import factorial
from scipy.optimize import fmin_l_bfgs_b as optim


os.chdir(os.path.expanduser("~/rough_code/"))

from fit_nbinom import *  # authored: Gokcen Eraslan


# check means in our real data
dataloc = os.path.expanduser("~/proj_LDHnet202204/data4py/")
rawqx = pd.read_csv(dataloc + "p4py_rawq.tsv", sep="\t",
                   index_col=0) # expression

rawqx.columns

means_obs = []
for kc in rawqx.columns:
    tmp = np.mean(rawqx[kc])
    means_obs.append(tmp)
    
    
#plt.plot(means_obs, "o", c="darkgreen")
#plt.title("observed mean of raw counts across subjects")
#plt.show()

conditions = [i.split("-")[0]+"-"+i.split("-")[1] for i in rawqx.columns]
oxygenlev = [i.split("-")[1] for i in rawqx.columns]
x = np.arange(1, len(conditions)+1,1)
xdf = pd.DataFrame(data={'subjects' : x, 
                         'mu_counts' : means_obs,
                         'condition' : oxygenlev})
xdf = xdf.sort_values('mu_counts')
xdf['subjects'] = x

sns.catplot(x="subjects",y="mu_counts", data=xdf, hue='condition')

i = 0

print("\nreal (mean, variance) for first observed subject")
print((np.mean(rawqx[rawqx.columns[i]]), np.var(rawqx[rawqx.columns[i]])))

np.sum(rawqx[rawqx.columns[i]])

# as E(X) = k / p , save p and k for which the likelihood of observed distr is max
# using Gokcen's fit_nbinom code:
    

ngenes = 30000
kp = fit_nbinom(rawqx[rawqx.columns[i]])
print(kp)
#print(kp['prob'])
#print(kp['size'])

singlesim_v = nbinom.rvs( kp['size'], kp['prob'], size = ngenes) # 30000 genes

print("\nestimated (mean, variance) by maximum likelihood ")
print((np.mean(singlesim_v), np.var(singlesim_v)))

#plt.plot([1,1])
plt.hist( rawqx[rawqx.columns[i]],  color = 'blue');plt.title("obs")
plt.show()
plt.hist(  singlesim_v,  color = "purple"); plt.title("simul")
plt.show()

def getvals(alist, maxv):
    ol_ = [i for i in alist if i <= maxv]
    return ol_
    
cutoffc = 10 # only for viz purposes
plt.hist( getvals(rawqx[rawqx.columns[i]], cutoffc),  color = 'blue', bins=30)
plt.title(f"obs, raw counts <= {cutoffc}")
plt.show()
plt.hist(  getvals(singlesim_v, cutoffc),  color = "purple", bins=30)
plt.title(f"simul, counts <= {cutoffc}")
plt.show()

cutoffc = 10 # only for viz purposes
plt.hist( getvals(rawqx[rawqx.columns[i]], cutoffc),  color = 'blue', bins=30)
plt.title(f"obs, raw counts <= {cutoffc}")
plt.show()
plt.hist(  getvals(singlesim_v, cutoffc),  color = "purple", bins=30)
plt.title(f"simul, counts <= {cutoffc}")
plt.show()

probs_ = np.arange(kp['prob'], 0., 0.01 )  #0.0277
nn_ = np.arange(kp['size'],0.05, 0.001 ) # 0.00001

nn_ = []
print(len(nn_))
print(nn_)
print(len(probs_))
plt.plot([1,1])
for v in nn_:
    for w in probs_:
        try:
            singlesim_v = nbinom.rvs(n = v, p = w, size = 300)
            plt.hist( singlesim_v,  color = "purple", bins=30)
            plt.title(f"simul, counts <= {cutoffc}, p = {w}, n={v}")
            plt.show()
        except ValueError as e:
            print(e)
            
    

# p (X = n ) = (n-1 k-1) (1-p )**n-k * p **k  is the counts distribution law
# E(X) = k / p
# V(X) = k(1-p) / p**2

"""
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
"""
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

#dist = getattr(scipy.stats, "gamma")