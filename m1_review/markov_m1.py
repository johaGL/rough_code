#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

doing tiny exercices related to Markov processes :
 A. one modified from 
  https://towardsdatascience.com/markov-chain-analysis-and-simulation-using-python-4507cee0b06e
   using https://github.com/FredericBourgeon/python-exercises/blob/master/Markov%20chain/markov-chain-1.py
 
 B. the one  shown in google colab:
  https://colab.research.google.com/drive/1d1wCFQe0mMu4io6G7E-xS4UEYaHFU3rY#scrollTo=G1cAK5Mu5fI_

**
see also wonderful explanation, equations, etc:
    https://personal.math.ubc.ca/~tbjw/ila/dds.html
    https://www.math.drexel.edu/~jwd25/LM_SPRING_07/lectures/Markov.html
**

intro to graph deep learning : https://ericmjl.github.io/essays-on-data-science/machine-learning/graph-nets/
with 'message passing' concept (also intro with Markov models)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

""" 
A part
"""
# here an adjacency matrix, unweighted, just boolean if edge exists
# (undirected)
adjM = np.array([[0,1,0,1],[1,0,1,1], [0, 1, 0 ,0], [1,1,0,0]]) 
Z = adjM / adjM.sum(axis=0)
print(Z)

# note that in this example columns sum up 1, 
# in counterpart, slide 12 by Laurent Tichit has rows summing up 1
# anyway let's continue:
    
#  initial probabilities (T0) distribution, let's say :
z0 = np.array([[0.1, 0.2, 0.1, 0.6]] )

# iterations
niter = 20
adfplo = pd.DataFrame(np.zeros([niter+1, Z.shape[0]]))

adfplo.columns = ['a','b','c','d']
iprobvec = z0.copy() # initialization i prob vector, when t==0
print("\nHERE WE GOOOOO!!!!!!!!!!!!!")
for t in range(niter+1):
    # stock vector at t
    adfplo.iloc[t,] = iprobvec
    # calculate for t+1
    iprobvec = np.dot(iprobvec, Z)
    print(iprobvec)

print("ended iterations. ==>", iprobvec/iprobvec.sum())
# with plt not necessary to melt or stack
#dfplotend = adfplo.reset_index().melt(id_vars = 'index')
#dfplotend.columns = ['iter', 'vertex', 'p']

adfplo.index.name = "time steps "
adfplo.columns.name = "vertex"
adfplo.plot(title = f"Toy example, {niter} iterations (time is discrete)")

## eigenvector and eigen values calc, using FredericBourgeon's code:
w, v = np.linalg.eig(Z.T) # w is the eigenvector OF THE TRANSPOSED Z here
j_stationary = np.argmin(abs(w - 1.0))
p_stationary = v[:,j_stationary].real
p_stationary /= p_stationary.sum()
p_stationary
print("linear algebra solution ==>" , p_stationary)

##########################################
# some nice reminders:
##########################################    
# diagonalization :https://bvanderlei.github.io/jupyter-guide-to-linear-algebra/Diagonalization.html
# A is a nxn matrix. 
# so D is a diagonal matrix with the eigenvalues of A 
# and S is the nxn matrix with the eigenvectors of A as its columns
# A = SDSinv
# and AS = SD
# let's check:
    
A = Z.copy()
eve, eva = np.linalg.eig(A) # DO NOT TRANSPOSE HERE
def dodiagmat(avec):
    foo = np.zeros((len(avec),len(avec)))
    for i in range(len(avec)):
        foo[i,i] = avec[i]
    return foo
            
D = dodiagmat(eve)
S = eva.copy()

print(A.dot(S))
print(S.dot(D))

try:
    Sinv = np.linalg.inv(S)
except:
    print("Inverse not possible")
Sinv    
SD = S.dot(D)
SDSinv = SD.dot(Sinv)
print(SDSinv)

SDSinver = np.around(SDSinv, decimals = 2)
print("\nhere our calculated 'A' from SDSinv equation")
print(SDSinver)
print("\nhere real A")
print(A)
# ha = np.dot(eigres[0],eigres[1])
# print(ha)


"""
B part
"""
# transition matrix, columns must sum up 1
# P = np.array([[0,0.25,0.25],[0.5,0.5,0.25],[0.5,0.25,0.5]])
# print(P)
    
