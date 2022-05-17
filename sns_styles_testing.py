#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 16:17:19 2022

@author: johanna
"""
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# https://github.com/mwaskom/seaborn/issues/672  == styles
# https://seaborn.pydata.org/tutorial/aesthetics.html
#sns.set(font_scale=1.5)

data = np.random.normal(size=(20, 6)) + np.arange(6) / 2

"""
f, ax = plt.subplots()
sns.violinplot(data=data)
sns.set_style({'font.family': 'serif', 'font.serif':['Times New Roman']})
sns.despine(offset=0, trim=True)
plt.show()
"""

f, ax = plt.subplots()
sns.violinplot(data=data)
sns.set_theme({ 'font.family' : 'sans-serif', 
               'font.serif' : ['Nimbus Roman No9 L'],
               })
sns.despine(offset=0, trim=True)
plt.show()


f, ax = plt.subplots()
sns.set_style({ 'font.family' : 'sans-serif', 
               'font.serif' : ['Nimbus Roman No9 L'],
               })



l_ =  ['Bitstream Vera Serif',
                         'DejaVu Serif',
                         'New Century Schoolbook',
                         'Century Schoolbook L',
                         'Utopia',
                         'ITC Bookman',
                         'Bookman',
                         'Nimbus Roman No9 L',
                         'Times New Roman',
                         'Times',
                         'Palatino',
                         'Charter',
                         'serif']
for thefunk in l_ :
    afo = [thefunk]
    f, ax = plt.subplots()
    sns.violinplot(data=data).set(title = thefunk)
    sns.set_style({ 'font.family': 'serif', 
                   'font.serif': afo, 
                   'font.size' : 11.0})
    sns.despine(offset=0, trim=True)

l_ =  []
for thefunk in l_ :
    afo = [thefunk]
    f, ax = plt.subplots()
    sns.violinplot(data=data).set(title = thefunk)
    sns.set_style({ 'font.family': 'sans-serif', 
                   'font.serif' :
                       , 'font.size' : 11.0})    
    sns.despine(offset=0, trim=True)

# other to see : 
    #https://stackoverflow.com/questions/41325160/seaborn-plots-in-a-loop
    # https://people.duke.edu/~ccc14/sta-663-2017/06_Graphics.html
    # https://stackoverflow.com/questions/36220829/fine-control-over-the-font-size-in-seaborn-plots-for-academic-papers