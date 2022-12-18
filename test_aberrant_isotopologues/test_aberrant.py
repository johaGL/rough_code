#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 15:27:17 2022

@author: johanna
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# def orderby_negative_values(df):
#     df = df.assign(max_neg=df.min(axis=1))
#     df.sort_values(by='max_neg', ascending = False, inplace=True) # to make bad be bottom
#     ordering_isos = df.index
#     ord_k = 
#     ordering_metabolites =
#     return df
    


def add_metabolite_column(df):
    theindex = df.index
    themetabolites = [i.split("_m+")[0] for i in theindex]
    df = df.assign(metabolite=themetabolites)
    
    return df


def add_isotopologue_type_column(df):
    theindex = df.index
    preisotopologue_type = [i.split("_m+")[1] for i in theindex]
    theisotopologue_type = [int(i) for i in preisotopologue_type]
    df = df.assign(isotopologue_type=theisotopologue_type)
    
    return df



def givelevels(melted):
    another = melted.copy()
    another = another.groupby('metabolite').min()
    another = another.sort_values(by='value', ascending = False)
    levelsmetabolites = another.index
    tmp = melted['metabolite'] 
    melted['metabolite'] = pd.Categorical(tmp, categories=levelsmetabolites)

    return melted



def table_minimalbymet(melted, fileout):
    another = melted.copy()
    another = another.groupby('metabolite').min()
    another = another.sort_values(by='value', ascending = False)
    another.to_csv(fileout,sep='\t', header=True)
    
    
    

# open files
cell = "correctedIsotopologues_myelin_cell.tsv"  
supernatant = "correctedIsotopologues_myelin_sn.tsv"

cell_df = pd.read_csv(cell ,sep='\t', header=0, index_col=0)
sn_df = pd.read_csv(supernatant ,sep='\t', header=0, index_col=0)


# cell plot: 
    
#cell_df = orderby_negative_values(cell_df)
#print(cell_df['max_neg'])



cell_df = add_metabolite_column(cell_df)
cell_df = add_isotopologue_type_column(cell_df)
print(cell_df.shape)
cellmelt = pd.melt(cell_df, id_vars=['metabolite', 'isotopologue_type'])


cellmelt = givelevels(cellmelt)

table_minimalbymet(cellmelt, "minvalIso_cell.tsv")

#cellmelt = cellmelt.iloc[0:100]

fig, ax = plt.subplots(1,1, figsize=(16, 10))
sns.stripplot(ax=ax, data=cellmelt, x="value", y="metabolite", jitter=False,
              hue="isotopologue_type", size=4, palette="tab20")

plt.axvline(x=0, 
          ymin=0, 
          ymax= 1,
          linestyle="--", color="gray")

plt.axvline(x=1, 
          ymin=0, 
          ymax= 1,
          linestyle="--", color="gray")
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
plt.title("intracellular compartment")
plt.xlabel("fraction")
plt.savefig("cell.png")




## supernatant plot: 

#sn_df = orderby_negative_values(sn_df)
#print(sn_df['max_neg'])    

sn_df = add_metabolite_column(sn_df)
sn_df = add_isotopologue_type_column(sn_df)
print(sn_df.shape)
snmelt = pd.melt(sn_df, id_vars=['metabolite', 'isotopologue_type'])

snmelt = givelevels(snmelt)

fig, ax = plt.subplots(1,1, figsize=(15, 10))
sns.stripplot(ax=ax, data=snmelt, x="value", y="metabolite", jitter=False,
              hue="isotopologue_type", size=4, palette="tab20")

plt.axvline(x=0, 
          ymin=0, 
          ymax= 1,
          linestyle="--", color="gray")

plt.axvline(x=1, 
          ymin=0, 
          ymax= 1,
          linestyle="--", color="gray")
plt.title("supernatant")
plt.xlabel("fraction")
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
plt.savefig("sn.png")



table_minimalbymet(snmelt, "minvalIso_sn.tsv")