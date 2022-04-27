#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 14:50:36 2022

@author: johanna

Resolves problem of dgTMatrix not being imported with barcodes and gene names
as happens with saunders data

runs inside bash : ./gz2BCgenes.sh 

"""
import os 
import sys


def frommultiline2list(fullfile):
    """
    input :  .txt file, has multiple lines : 
    %%GENES	Zfp583	Zfp59	Zfp592	Zfp593
    or
    %%CELL_BARCODES	P60FCRep2P2_AGAGACCTGATA
    
    output :  list
    """
    oli = []
    with open(fullfile, "r") as f:
        databc = f.readlines()   
    for k in databc:
        tmp = k.strip().split("\t")
        oli += tmp[1:] #ignores the 0 index (%%GENES or %%CELLBARCODES)
    return oli

def fromlist2txt(mylist, fullfile):
    txtout = ""
    for k in mylist:
        txtout += k+"\n"
    with open(fullfile, "w") as f:
        f.write(txtout)
        

if __name__ == "__main__" :   
    print(sys.argv[1])
    datadir = os.path.expanduser(sys.argv[1])
    thefile = sys.argv[2]
    fullfile = datadir+thefile
    fromlist2txt(frommultiline2list(fullfile), fullfile)
    
    # #simple test:
    # datadir = os.path.expanduser("~/dataref_scRNA/saunders/")
    # thefile = "frontalbarcodes.txt"
    # fullfile = datadir+thefile   
    # theli = frommultiline2list(fullfile)
    # fromlist2txt(theli, datadir+"hehehe")
    print("ok")

         

