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

sns.reset_orig()
"""
f, ax = plt.subplots()
sns.violinplot(data=data)
sns.set_style({'font.family': 'serif', 'font.serif':['Times New Roman']})
sns.despine(offset=0, trim=True)
plt.show()
"""

# f, ax = plt.subplots()
# sns.violinplot(data=data)
# sns.set_theme({ 'font.family' : 'sans-serif', 
#                'font.serif' : ['Nimbus Roman No9 L'],
#                })
# sns.despine(offset=0, trim=True)
# plt.show()


# f, ax = plt.subplots()
#sns.set_style({ 'font.family' : 'sans-serif', 
 #              'font.serif' : ['Nimbus Roman No9 L'],
  #             })

f, ax = plt.subplots()
sns.violinplot(data=data).set(title = "nochange no Serif choice")
sns.despine(offset=0, trim=True)

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
                   'font.serif': afo})
    sns.despine(offset=0, trim=True)

l_ =  ['Arial',
  'Liberation Sans',
  'Bitstream Vera Sans',
  'sans-serif']
for thefunk in l_ :
    afo = [thefunk]
    f, ax = plt.subplots()
    sns.violinplot(data=data).set(title = thefunk)
    sns.set_style({ 'font.family': 'sans-serif', 
                   'font.sans-serif' : afo
         })    
    sns.despine(offset=0, trim=True)
    plt.savefig("p_"+afo[0]+".pdf", format = "pdf")


# other to see : 
    #https://stackoverflow.com/questions/41325160/seaborn-plots-in-a-loop
    # https://people.duke.edu/~ccc14/sta-663-2017/06_Graphics.html
    # https://stackoverflow.com/questions/36220829/fine-control-over-the-font-size-in-seaborn-plots-for-academic-papers
    
    
# First create some toy data:
x = np.linspace(0, 2*np.pi, 400)
y = np.sin(x**2)

# Creates just a figure and only one subplot


# Creates two subplots and unpacks the output array immediately
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1.plot(x, y)
ax1.set_title('Sharing Y axis')
ax2.scatter(x, y)

# using facet grid : i dont like the titles it sets
import pandas as pd
iris = sns.load_dataset("iris")
iris_long = pd.melt(iris, "species", var_name="measurement")
g = sns.FacetGrid(iris_long, hue="species", col="measurement", col_wrap=2, sharex=False)
g.map(plt.hist, "value", alpha=.4)


def complextest():
    sns.set_style({ 'font.family': 'sans-serif', 
                     'font.sans-serif' : 'Liberation Sans'   }) 

    tips = sns.load_dataset('tips')
    agg_tips = tips.groupby(['day', 'sex'])['tip'].sum().unstack().fillna(0)
    totals = [ 100/(i + j)  for (i,j) in zip(agg_tips['Male'], agg_tips['Female']) ]
    copyy = agg_tips.copy()
    copyy['Male'] = [ i*j for (i,j) in zip(copyy["Male"], totals) ]
    copyy['Female'] = [ i*j for (i,j) in zip(copyy["Female"], totals) ]

    verify = copyy['Male'].map(float) + copyy['Female'].map(float)

    agg_tips = copyy

    def dobarplot( agg_tips) : 
        fig, ax = plt.subplots()
        
        colors = ['#24b1d1', '#ae24d1']
        bottom = np.zeros(len(agg_tips))
        for i, col in enumerate(agg_tips.columns):
          ax.bar(
              agg_tips.index, agg_tips[col], bottom=bottom, label=col, color=colors[i])
          bottom += np.array(agg_tips[col])
        
        totals = agg_tips.sum(axis=1)
        y_offset = 4
          
        # Let's put the annotations inside the bars themselves by using a
        # negative offset.
        y_offset = -15
        # For each patch (basically each rectangle within the bar), add a label.
        for bar in ax.patches:
          ax.text(
              # Put the text in the middle of each bar. get_x returns the start
              # so we add half the width to get to the middle.
              bar.get_x() + bar.get_width() / 2,
              # Vertically, add the height of the bar to the start of the bar,
              # along with the offset.
              bar.get_height() + bar.get_y() + y_offset,
              # This is actual value we'll show.
              round(bar.get_height()),
              # Center the labels and style them a bit.
              ha='center',
              color='w',
              weight='bold',
              size=8
          )
        ax.set_title('Tips by Day and Gender')
        # ax.legend() no legend because goes into grid
        return ax

    #hey = dobarplot(agg_tips)

    #fig.savefig("youpi", format="pdf")

    f, axs = plt.subplots(1, 2, sharey=True)
    meh = 0
    for ax in axs: 
        #dobarplot(agg_tips=agg_tips)
        #fig, ax = plt.subplots()
        
        colors = ['#24b1d1', '#ae24d1']
        bottom = np.zeros(len(agg_tips))
        for i, col in enumerate(agg_tips.columns):
          ax.bar(
              agg_tips.index, agg_tips[col], bottom=bottom, label=col, color=colors[i])
          bottom += np.array(agg_tips[col])
        
        totals = agg_tips.sum(axis=1)
        y_offset = 4
          
        # Let's put the annotations inside the bars themselves by using a
        # negative offset.
        y_offset = -15
        # For each patch (basically each rectangle within the bar), add a label.
        for bar in ax.patches:
          ax.text(
              # Put the text in the middle of each bar. get_x returns the start
              # so we add half the width to get to the middle.
              bar.get_x() + bar.get_width() / 2,
              # Vertically, add the height of the bar to the start of the bar,
              # along with the offset.
              bar.get_height() + bar.get_y() + y_offset,
              # This is actual value we'll show.
              round(bar.get_height()),
              # Center the labels and style them a bit.
              ha='center',
              color='w',
              weight='bold',
              size=8
          )
        ax.set_title(str(meh))

        # ax.legend() no legend because goes into grid
        meh += 1
    f.suptitle("CELL CONTROL", fontsize=  14)
    return "ended complex test"


complextest()

##################################################
# grid seaborn
gridsel = ["L-Lactic_acid" , "Citric_acid" , "Oxoglutaric_acid",       
                "Succinic_acid",           
                "L-Malic_acid", 
                "L-Alanine",
                "Fumaric_acid",
                "Pyruvic_acid"  , "L-Glutamic_acid", "L-Glutamine", 
                "L-Aspartic_acid", 
                "L-Asparagine"]
df4plot
# neee  = df4plot.loc[df4plot["metabolite"] == "L-Malic_acid"]
# print(neee.dtypes)
# ["Isotopologue Contribution (%)"]
# #neee["ouf" ] = pd.to_numeric(neee["Isotopologue Contribution (%)"])
# neee.dtypes
# #naaa = neee.groupby(["condition", "metabolite",  "m+x", "timepoint"] ).mean(numeric_only = True)

# naaa = neee.groupby(["condition", "metabolite",  "m+x", "timepoint"] ).mean()
# naaa = naaa.reset_index()
# naaa['shit'] = [str(i) for i in naaa["timepoint"].to_list()]


griddata = df4plot.loc[df4plot["metabolite"].isin(gridsel),]
griddata = griddata.groupby(["condition","metabolite", "m+x", "timepoint"]).mean()
griddata = griddata.reset_index()

grid = sns.FacetGrid(griddata, row="condition", col="metabolite",  
                     margin_titles=True)
sns.catplot(
         data = griddata,
         x= "timepoint", y="Isotopologue Contribution (%)",
         hue = "m+x",
         kind= "point",
         row="condition", col="metabolite",  
                             margin_titles=True)

grid = sns.FacetGrid(griddata, row="condition", col="metabolite",  
                     margin_titles=True)

grid.map(sns.histplot,   data=griddata,     
         x= "timepoint", weights="Isotopologue Contribution (%)",
         hue = "m+x",
         multiple = "stack",
         legend = True
         )
         
#order = list(set(griddata["timepoint"])) 