#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:10:40 2022

@author: johanna
"""
import os
import scipy
from scipy.stats import nbinom
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from numpy.random import seed
seed(42)


# viz some distributions :
saz = 1000
n_ = [10, 15, 30, 50, 100]
p_ = np.arange(0.2, 1, 0.2)
fig, axs = plt.subplots(len(n_),len(p_), figsize=(13,17))
for v in range(len(n_)):
    for w in range(len(p_)):
        hypo = np.random.binomial(n=n_[v], p=p_[w], size=saz)
        axs[v, w].hist(hypo, color="gray", bins=40)
        axs[v,w].set_title(f"n = {n_[v]}, p = {np.round(p_[w],2)}")
#for ax in axs.flat:
#    ax.set()
#for ax in axs.flat:
#    ax.label_outer()
fig.subplots_adjust(hspace=0.6)
fig.suptitle("BINOMIAL DISTRIBUTION ", fontsize=  14)


lambdas = [ 5,  25, 50,  500]
fig, axs = plt.subplots(1,len(lambdas), figsize=(14,4))
for w in range(len(lambdas)):
    hypo = np.random.poisson(lam=lambdas[w], size=saz)
    axs[w].hist(hypo, color="gray", bins=40)
    axs[w].set_title(f" lambda = {lambdas[w]}")
for ax in axs.flat:
    ax.set()
for ax in axs.flat:
    ax.label_outer()
fig.suptitle("POISSON DISTRIBUTION\n", fontsize=  14)

# negative binomial
n_ = [15, 50, 100, 10000, 10000]
p_ = [0.3, 0.5, 0.7, 0.9]#np.arange(0.3, 1, 0.2)
fig, axs = plt.subplots(len(n_),len(p_), figsize=(13,17))

for v in range(len(n_)):
    for w in range(len(p_)):
        hypo = np.random.negative_binomial(n=n_[v], p=p_[w], size=saz)
        axs[v, w].hist(hypo, color="gray", bins=40)
        axs[v,w].set_title(f"n = {n_[v]}, p = {np.round(p_[w],4)}")

fig.subplots_adjust(hspace=0.6)
fig.suptitle("NEGATIVE BINOMIAL DISTRIBUTION ", fontsize=  14)


k_ = [0.5, 5, 10, 15,50]
t_ = [1, 5, 10]
fig, axs = plt.subplots(len(k_),len(t_), figsize=(13,17))
for v in range(len(k_)):
    for w in range(len(t_)):
        hypo = np.random.gamma(shape=k_[v], scale=t_[w], size=saz)
        axs[v, w].hist(hypo, color="gray", bins=40)
        axs[v,w].set_title(f"kappa = {k_[v]}, theta = {np.round(t_[w],2)}")
#for ax in axs.flat:
#    ax.set()
#for ax in axs.flat:
#    ax.label_outer()
fig.subplots_adjust(hspace=0.6)
fig.suptitle("GAMMA DISTRIBUTION ", fontsize=  14)
