#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 18:21:41 2022
Simulating counts for desired number of subjects
(negative binomial distributions)
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

# NEGATIVE BINOMIAL : 
# p (X = n ) = (n-1 k-1) (1-p )**n-k * p **k  is the counts distribution law
# E(X) = k / p
# V(X) = k(1-p) / p**2

NGENESFINAL = 300
NSUBJECTSFINAL = 15

# check means in our real data
dataloc = os.path.expanduser("~/proj_LDHnet202204/data4py/")
rawqx = pd.read_csv(dataloc + "p4py_rawq.tsv", sep="\t",
                   index_col=0) # expression

rawqx.columns

means_obs = []
vars_obs = []
for kc in rawqx.columns:
    tmp = np.mean(rawqx[kc])
    means_obs.append(tmp)
    vars_obs.append(np.var(rawqx[kc]))
    
    
#plt.plot(means_obs, "o", c="darkgreen")
#plt.title("observed mean of raw counts across subjects")
#plt.show()

conditions = [i.split("-")[0]+"-"+i.split("-")[1] for i in rawqx.columns]
oxygenlev = [i.split("-")[1] for i in rawqx.columns]
indivs = np.arange(1, len(conditions)+1, 1)
xdf = pd.DataFrame(data={'subjects' : indivs, 
                         'mu_counts' : means_obs,
                         'variance' : vars_obs,
                         'condition' : oxygenlev})
xdf = xdf.sort_values('mu_counts')
xdf['subjects'] = indivs

#sns.catplot(x="subjects",y="mu_counts", data=xdf, hue='condition')
sns.scatterplot(data=xdf, x="mu_counts", y="variance", hue="condition")
plt.title("Observed mean and variances across subjects\n")
plt.show()
i = 0

"""
estimating p and n parameters with MLE for negative binom
"""
# as E(X) = k / p , save p and k for which the likelihood of observed distr is max
# using Gokcen's fit_nbinom code:
 
ngenes = 30000 # only to test fit function
kp = fit_nbinom(rawqx[rawqx.columns[i]])
print(kp)
#print(kp['prob'])
#print(kp['size'])

fitcheck = nbinom.rvs( kp['size'], kp['prob'], size = ngenes) 

print("\nreal (mean, variance) for first observed subject")
print((np.mean(rawqx[rawqx.columns[i]]), np.var(rawqx[rawqx.columns[i]])))

np.sum(rawqx[rawqx.columns[i]])
print("\nestimated (mean, variance) by maximum likelihood ")
print((np.mean(fitcheck), np.var(fitcheck)))
print("\n")

#plt.plot([1,1])
plt.hist( rawqx[rawqx.columns[i]],  color = 'blue');plt.title("obs")
plt.show()
plt.hist(  fitcheck,  color = "purple"); 
plt.title(f"simul,\
   nbinom.rvs({np.round(kp['size'],3)},{np.round(kp['prob'],5)},size={ngenes})")
plt.show()


def truncateat(alist, minv,maxv):
    ol_ = [i for i in alist if (i <= maxv and i >= minv)]
    return ol_
    
cutoffc = 50 # only for viz purposes
plt.hist( truncateat(rawqx[rawqx.columns[i]], 1, cutoffc),  color = 'blue', bins=30)
plt.title(f"obs, raw counts truncated (<= {cutoffc})")
plt.show()
plt.hist(  truncateat(fitcheck, 1,  cutoffc), color = "purple", bins=30)
plt.title(f"simul, counts truncated (<= {cutoffc})")
plt.show()


print("CAUTION: with MLE, n was not an integer as expected.\
      Nevertheless, the distribution histogram is still coherent")


"""
testing other p, and bigger n
"""
probs_ = [0.01, 0.015, 0.020, 0.025] 
nn_ = [1, 2,3,4,5  ] 

mean_ = []
vari_ = []
pro_str = []
nparam = []
for v in nn_:
    for w in probs_:
        try:
            singlesim_v = nbinom.rvs(n = v, p = w, size = NGENESFINAL)
            #plt.hist( singlesim_v,  color = "orange", bins=30)
            #plt.title(f"simul, counts <= {cutoffc}, p = {w}, n={v}")
            #plt.show()
            mean_.append(np.mean(singlesim_v))
            vari_.append(np.var(singlesim_v))
            pro_str.append(str(w))
            nparam.append(v)
        except ValueError as e:
            print(e)
            
aodf = pd.DataFrame(data={"mean":mean_,
                          "variance":vari_,
                          "p_param": pro_str,
                          "n_param":nparam})

sns.scatterplot(data=aodf, x="mean", y="variance", 
                hue='p_param', size='n_param')
plt.title("Simulated distributions")
plt.show()

"""
select interval of mean and variance
(these intervals will serve to pick n and p)
"""
print("\nSelecting estimated parameters")

pulledpar = aodf[(aodf["mean"] >= 130)]
pulledpar = pulledpar[(pulledpar["mean"] <= 500)]
print(pulledpar)

#MINP = min([float(i) for i in pulledpar["p_param"]])
#MAXP = max([float(i) for i in pulledpar["p_param"] if float(i) <= 0.025])
MINP = 0.018
MAXP = 0.02

MIN_N = 3
MAX_N = 3
#MAX_N = max(pulledpar["n_param"])
#if max(pulledpar["n_param"]) > 5:
#    MAX_N = 5  #

print((MINP, MAXP))
print(MIN_N, MAX_N)

print("\n************* sampling parameters  **************")

print(f"desired number of simulated subjects is :{NSUBJECTSFINAL}")
nsubjsplus = NSUBJECTSFINAL + 0  # at the end retain desired
n_ = np.random.choice(np.arange(MIN_N,MAX_N+1,1), 
                      size=nsubjsplus, replace=True)
# to simplify, pull p from a uniform distribution
p_ = np.random.uniform(low=MINP, high=MAXP, size=nsubjsplus)

print(MIN_N)

p_ = np.sort(p_)
n_ = np.sort(n_)

print(p_)
print(n_)

outdf = pd.DataFrame(data={"p_param": p_,
                           "n_param": n_})

print("*************** simulating gene counts by indiv***********************")
arrfill = np.zeros([NGENESFINAL,nsubjsplus])
for s in range(nsubjsplus): # fill columns
    tmp = nbinom.rvs(n = n_[s], p = p_[s], size=NGENESFINAL)
    arrfill[:,s] = np.sort(tmp)
    #plt.hist(arrfill[:,s], color="gray")
    #plt.show()

meanarr = []
vararr = []
for s in range(nsubjsplus):
    meanarr.append(np.mean(arrfill[:,s]))
    vararr.append(np.var(arrfill[:,s]))
    
plt.plot(meanarr,vararr, "o")
plt.plot()

#picksubset with highest m
thevars = arrfill.var(axis=1)


print("normalizing (CPM)")



"""

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
# def initialworkaround():
#     n = 5
#     p = 0.01
#     nes = np.arange(1,3, 1)
#     for n in nes : 
#         fig, ax = plt.subplots(1,1)
    
#         x = np.arange(nbinom.ppf(0.001, n, p),
#                       nbinom.ppf(0.99, n, p))
        
#         ax.plot(x, nbinom.pmf(x, n, p), 'bo', ms=8, label='nbinom pmf')
#         ax.vlines(x, 0, nbinom.pmf(x,n,p), colors = 'b', lw=5, alpha=0.5)
    
#     n = 300
#     he = nbinom.rvs(1,p, size = 100)
#     fig, ax = plt.subplots(1,1)
#     plt.plot(he)

#dist = getattr(scipy.stats, "gamma")
# need seed : https://stackoverflow.com/questions/16016959/scipy-stats-seed


# def printisto(v, w):
#     plt.hist(x=nbinom.rvs(n=v, p=w, size=300),
#          color ="orange", bins=30)
#     plt.title(f"n ={v}   p= {w}")
#     plt.show()
    
# def printdots(v, w):
#     dotos = nbinom.rvs(n=v, p=w, size=300)
#     plt.plot(dotos,dotos, "o" ) #, color ="darkgreen", s=3)
#     plt.title(f"n ={v}   p= {w}")
#     plt.show()
# for v in [0.03, 0.3,0.4,0.5, 0.6]:
#     for w in [ 0.000189, 0.002]:
#         printisto(v,w)
#         printdots(v,w)


