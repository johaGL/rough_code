#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 15:04:34 2022

@author: johanna
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 15:34:38 2022

@author: johanna
"""

import os
import scipy
from scipy import stats
from scipy.stats import nbinom
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from numpy.random import seed
seed(42)

# Understand one single feature(gene) across many samples :
# this team https://www.biorxiv.org/content/10.1101/572693v1.full.pdf
# fitted all TGCA rna-seq genes to lognlizal, Cauchy and Gamma, Pareto:
# they found that 31% genes were well fitted by Gamma. All the others
# were likely to follow the other laws, each law in smaller percentage.
# this other team : https://www.biorxiv.org/content/10.1101/2022.02.14.480329v1.full
# models expression as a gaussian mixture
# an example of shapes simmulation (Australia): https://www.mdpi.com/2073-4425/11/10/1231/htm
 # an example of true shapes (logTPM) : https://rpubs.com/marsluo1127/expDistColon   
    
# check several genes at choice in: GlioVis
# for simplicity, simulate counts as negative binomial
saz = 300 # nb of subjects
n_ = [ 5, 10 ]
pr = 0.1 #np.arange(0.3, 1, 0.2)

icD = dict()

#two uncorrelated

fig, axs = plt.subplots(len(n_) + 1 ,1, figsize=(7,10))
genename = 0
for v in range(len(n_)):
    simul = []
    simul = np.random.negative_binomial(n=n_[v], p = pr, size=saz)
    print(len(simul))
    axs[v].hist(simul, bins=40)
    axs[v].set_title(f"n = {n_[v]}, p = {np.round(pr,4)}")
    icD[genename] = {"counts" : simul,
                     "nliz" : np.array(simul) / sum(simul) , 
                     "n" : n_[v],
                     "p" : pr }
    genename += 1

s1 = np.var(icD[0]["counts"])
s2 = np.var(icD[1]["counts"])
cov = np.sqrt(s1+s2) / (s1*s2)
cov
myrho = stats.spearmanr(icD[0]["nliz"],
                       icD[1]["nliz"], axis = 0)
axs[ len(n_) ].text( x= 0, y= 0.5, s=f'covariance: {cov}')
axs[ len(n_) ].text(x=0, y = 0.3, s = f"rho : {myrho}")
fig.subplots_adjust(hspace=0.2)
fig.suptitle("two uncorrelated ", fontsize=  14)
fig.show()
plt.show()
foo = "counts"
plt.scatter(x = icD[0][foo], y = icD[1][foo])
plt.show()

# Transform both genes to be correlated:
# https://oscarnieves100.medium.com/simulating-correlated-random-variables-in-python-c3947f2dbb10
# we will use oscar nieves explanation, even if intended for normally distributed variables. 
# So for our data, desired rho and obtained rho will differ (because not normal):

def newX(oldX, sigma, mu):
    vec = mu + (sigma * oldX)
    return vec

def newY(oldY, oldX, sigma, mu, rho):
    factor = (rho * oldX) + ( np.sqrt(1-rho**2) )  * oldY 
    vec = mu + (sigma * factor)
    return vec

rho = 0.9
np.mean(icD[0]["counts"])
fii = "counts"
X = newX(np.array(icD[0][fii]), 1.2, 40)
Y = newY(np.array(icD[1][fii]), np.array(icD[0][fii]), 3, 1000, rho)

print(stats.spearmanr(X,Y))
plt.scatter(X,Y)
