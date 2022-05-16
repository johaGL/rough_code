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
# fitted all TGCA rna-seq genes to logNormal, Cauchy and Gamma, Pareto:
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
p_ = [0.2,   0.9]#np.arange(0.3, 1, 0.2)

icD = dict()


fig, axs = plt.subplots(len(n_),len(p_), figsize=(13,17))
genename = 0
for v in range(len(n_)):
    for w in range(len(p_)):
        hypo = np.random.negative_binomial(n=n_[v], p=p_[w], size=saz)
        axs[v, w].hist(hypo, color="gray", bins=40)
        axs[v,w].set_title(f"n = {n_[v]}, p = {np.round(p_[w],4)}")
        icD[genename] = {"counts" : hypo,
                         "n" : n_[v],
                         "p" : p_[w] }
        genename += 1

fig.subplots_adjust(hspace=0.6)
fig.suptitle("NEGATIVE BINOMIAL DISTRIBUTION ", fontsize=  14)
fig.show()
plt.show()
print("ended demo distribution")

newicD = dict()
bigene = [0,1,2]
for k in icD.keys():
    newicD[bigene[0]] = icD[k]["counts"]
    randomfactors = np.random.choice(np.arange(0.2,2,0.2), len(newicD[bigene[0]]))
    artifi1 = [int(i/j) for i,j in zip(newicD[bigene[0]], randomfactors)]
    randomfactorsB = np.random.choice(np.arange(100,200, 1), len(newicD[bigene[0]]))
    artifi2 = [int(i+j) for i,j in zip(newicD[bigene[0]], randomfactorsB)]
    newicD[bigene[1]] = np.array(artifi1)
    newicD[bigene[2]] = np.array(artifi2)
    bigene = [i+3 for i in bigene]

# ## artificiall thing : 
# randomfactors = np.random.choice(np.arange(0.2,0.3, 0.05), scou.shape[1])
# to_append = [ i*j for i,j in zip(scou.iloc[0], randomfactors) ] 
# a_series = pd.Series(to_append)
# print(a_series.T.shape)
# scou = pd.concat([scou.T, a_series], axis=1).T
# scou.index = range(scou.shape[0])

mean = []
vari = [] 

for k in newicD.keys():
    mean.append( np.mean(newicD[k]))
    vari.append( np.var(newicD[k]))


plt.plot(mean, vari, "o", color="darkgreen")
plt.xlabel("mean")
plt.ylabel("variance")
plt.title("Simulated genes (a dot == a gene), {} values per gene".format(saz) )
plt.show()


# fill scou
scou = np.zeros( (len(newicD.keys()), saz) )
for gk in newicD.keys():
    scou[gk,] = newicD[gk]
    print(sum(newicD[gk]))
scou = pd.DataFrame(scou)
scou.index = newicD.keys()




cpmdf = scou.copy(deep=False)
# "normalize" by individual (across columns)
for i in cpmdf.columns:
    denomin = sum(cpmdf[i])
    cpm = (np.array(cpmdf[i]) / denomin  ) * 1e6
    cpmdf[i] = cpm
 
  
def meltcustom(adf, newcolnames):
    """
    adf : all rows to be piled up
    newcolnames : eighter ["gene", "counts"] or ["gene","cpm"]
    """
    columnA = []
    columnB = []
    for i, row in adf.iterrows():
        for k in row:
            columnB.append(k)
            columnA.append(i)
    apile = pd.DataFrame(np.array([columnA, columnB]).T)
    apile.columns = newcolnames
    return apile


countmelt = meltcustom(scou, ["gene", "counts"])
sns.violinplot(x="gene", y="counts", data=countmelt)
plt.show()


countmelt["logcounts"] = np.log10(countmelt["counts"] + 1)
sns.violinplot(x="gene", y="logcounts", data=countmelt)
plt.show()

meltedcpm = meltcustom(cpmdf, ["gene", "cpm"])
meltedcpm["logcpm"] = np.log10(meltedcpm["cpm"] + 1)
sns.violinplot(x= "gene", y="logcpm", data=meltedcpm)
plt.show()


def calc_spearman(adf):
    """
    adf : genes in rows, subjects in columns
    """
    matcorr = np.zeros((adf.shape[0], adf.shape[0]))
    pvalue = np.zeros((adf.shape[0], adf.shape[0]))
    pvalue[:] = np.nan
    for i in range(adf.shape[0]):
        for j in range(adf.shape[0]):
            geneAx = adf.iloc[i,].tolist()
            geneBx = adf.iloc[j,].tolist()
            resv = stats.spearmanr(geneAx, geneBx , axis = 0)
            matcorr[i,j] = resv.correlation
            pvalue[i,j] = resv.pvalue
    return matcorr, pvalue



corrcpm, pvalcpm = calc_spearman(cpmdf)
corrco, pvalco = calc_spearman(scou)

corrva = corrcpm[corrcpm < 0]
print(min(set(corrva)))
# ## artificiall thing (add arbritrary scalar (60) to a vector ):
#     # plot was done, yes it says they are correlated : 1
#     # no need to do this anymore
# to_append = [ i + 60 for i in scou.iloc[0] ] 
# a_series = pd.Series(to_append)
# print(a_series.T.shape)
# scou = pd.concat([scou.T, a_series], axis=1).T
# cpmdf = scou.copy(deep=True)

# tips = sns.load_dataset("tips")
# sns.set_style('ticks')
# g = sns.FacetGrid(tips, col = 'time')
# g = g.map(plt.hist, "tip")
# g.savefig( "hehehe.pdf", format="pdf")


