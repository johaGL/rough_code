#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 16:19:28 2022

@author: johanna
"""

import numpy as np
import pandas as pd

# example Random Walk ("classical")

A = np.array([[0,1,0,0],
              [1,0,1,1],
              [0,1,0,1],
              [0,1,1,0]])
A
# transition matrix
# by row normalization:
R = A / A.sum(axis=1) [:,np.newaxis]   
# by column normalization:
M = A / A.sum(axis=0)  

# lets continue using column normalized transition matrix
print(M)
# probability distribution at time zero
p_0 = np.array([1,0,0,0])


# iterate and put values in table:
iters = 12
p_t = p_0.copy()
res = np.zeros((iters,M.shape[0]))
res[0,:] = p_t
for t in range(iters):
    tmp = M.dot(p_t.T)
    p_t = tmp.copy()
    res[t,:] = p_t
    del(tmp)
    
adfplo = pd.DataFrame( res, columns = ['1','2','3','4'])
adfplo.index.name = "time steps "
adfplo.columns.name = "vertex"
adfplo.plot(title = f"Toy example, {iters} iterations",
            alpha = 0.5, linewidth=3)

susu = adfplo["suma"] = adfplo.sum(axis=1) 

for i in  ['1','2','3','4']:
    print(np.sum(adfplo[i].to_numpy()))
    
adfplo['1'].sum()
