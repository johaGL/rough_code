#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 12:29:46 2022

@author: johanna
"""

import os
import pandas as pd
import numpy as np
#import seaborn as sn
#import plotly.express as px
import matplotlib.pyplot as plt

dataloc = "/home/johanna/proj_LDHnet202204/data4py/"
gmx = pd.read_csv(dataloc + "p4py_geomdf.tsv", sep="\t",
                   index_col=0)
emx = pd.read_csv(dataloc + "p4py_cpm.tsv", sep="\t",
                   index_col=0)
mtd = pd.read_csv(dataloc + "p4py_metadata.tsv", 
                      sep="\t")


print(3 * len(set(mtd['condition'])))
print(mtd.shape)

myconds = ["sgControl-Hypox", "sgLDHAB-Hypox" ]

selsmp = mtd.loc[mtd['condition'].isin(myconds), "sample"].tolist()
emxsub = emx[selsmp]
print(emxsub.head)
print(emxsub.columns)

gmxsub = gmx[myconds]
print(gmxsub.head)
print(gmxsub.columns)
# plotting the geometric means of these stuff
loda = np.log10(np.array(gmxsub))
loda = pd.DataFrame(loda, columns = myconds)
plt.scatter(myconds[0],
            myconds[1], data=loda)
plt.close()
plt.scatter(myconds[0], myconds[1],
            data=gmxsub)
plt.xlim([-1,100])
plt.ylim([-1,100])
plt.show()


