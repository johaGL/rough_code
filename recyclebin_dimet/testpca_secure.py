import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.decomposition import PCA
import seaborn as sns
import math
from scipy import stats
from matplotlib.patches import Ellipse


# calculate Ellipses , many thanks to :
# https://rayblick.gitbooks.io/my-python-scrapbook/content/analysis/plotting/scatterplot_ellipse.html

def make_ellipse(mean, cov,  level=0.95, color=None):
    """Support function for scatter_ellipse."""
    from matplotlib.patches import Ellipse

    v, w = np.linalg.eigh(cov)
    u = w[0] / np.linalg.norm(w[0])
    angle = np.arctan(u[1]/u[0])
    angle = 180 * angle / np.pi # convert to degrees
    v = 2 * np.sqrt(v * stats.chi2.ppf(level, 2)) #get size corresponding to level
    ell = Ellipse(mean[:2], v[0], v[1], 180 + angle, facecolor='none',
                  edgecolor="gray", alpha=0.4,
                  #ls='dashed',  #for debugging
                  lw=1.5)
    return ell

def eigsorted(cov):
    """
    used for calculating ellipses
    many thanks to :
https://rayblick.gitbooks.io/my-python-scrapbook/content/analysis/plotting/scatterplot_ellipse.html
    """
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]


def calcPCAand2Dplot(mymat, metadata, col1, col2, pointlabels, title, odir, *ellipsesparams):
    """
    ellipsesparams is a tuple with info to draw the ellipses, examples:
          ("timepoint")    or   ( "condition")
    """
    print(ellipsesparams)
    col_ellipses = None
    if len(ellipsesparams) > 0 :
        col_ellipses = ellipsesparams[0]

    desdim = 4
    X = np.transpose(np.array(mymat))
    pca = PCA(n_components=desdim)
    pc = pca.fit_transform(X)

    pc_df = pd.DataFrame(data=pc,
                         columns=['PC' + str(i) for i in range(1, desdim + 1)])
    pc_df = pc_df.assign(sample = mymat.columns)

    pc_df = pd.merge(pc_df, metadata, on= "sample")
    print(pca.explained_variance_)
    print(pca.explained_variance_ratio_ * 100)
    dfvare = pd.DataFrame({
        # 'var': pca.explained_variance_ratio_,
        'var': pca.explained_variance_ratio_ * 100,
        'PC': ['PC' + str(i) for i in range(1, desdim + 1)]})
    print(dfvare.head())
    # barplot
    plt.figure()
    sns.barplot(x='PC', y="var", data=dfvare, color="cadetblue")
    plt.title("Principal components: explained variance")
    plt.ylabel("variance (%)")
    plt.savefig(odir+"variance_pca_"+title+".pdf")
    plt.close()


    # scatterplot

    # sns.set_style("whitegrid")
    fig = plt.figure()

    sns.scatterplot(x="PC1", y="PC2",
                    data=pc_df,
                    hue=col1,
                    style=col2,
                    legend=True,
                    s=80, zorder=3)
    plt.axhline(0, ls="--", color="gray",zorder = 1)
    plt.axvline(0, ls="--", color="gray", zorder = 1)

    yesnolabel = "no"
    if pointlabels != "":
        yesnolabel = "yes"
        for l in range(0, pc_df.shape[0]):
            plt.text(pc_df.PC1[l] + 0.5, pc_df.PC2[l]+0.1, pc_df[pointlabels].tolist()[l],
                                    horizontalalignment='center', size='x-small')
    # end if
    plt.xlabel(f"{dfvare.iloc[0, :]['PC']} {round(dfvare.iloc[0, :]['var'], 2) } %")
    plt.ylabel(f'{dfvare.iloc[1, :]["PC"]} {round(dfvare.iloc[1, :]["var"], 2) } %')
    plt.title(title)

    # ellipses, if true
    if col_ellipses is not None:
        myellipsesnames = pc_df[col_ellipses].unique()
        for lab in myellipsesnames:
            xdata = pc_df.loc[pc_df[col_ellipses] == lab, 'PC1']
            ydata = pc_df.loc[pc_df[col_ellipses] == lab, 'PC2']


            # get values to build the ellipse
            cov = np.cov(xdata, ydata)
            ell = make_ellipse((0,0),cov,
                               level= 0.3, color=None)
            # vals, vecs = eigsorted(cov)
            # theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
            # print(vals, vecs, theta)
            # w, h = 2 * 2 * np.sqrt(vals)
            #
            # # create the ellipse
            # ell = Ellipse(xy=(np.mean(xdata), np.mean(ydata)),
            #               width=w, height=h,
            #               angle=theta, color='gray', alpha=0.2)

            # ell.set_facecolor() # reference the colour for each factor
            fig.add_artist(ell)

    # end  if ellipses
    plt.savefig(f"{odir}pca_{title}_label{yesnolabel}.pdf", format="pdf")
    return 0

metadata = pd.read_csv("data/metadata_myelin.csv", header = 0)




tmpdir = "tmp/"
compartments_names = {"cells" : "cell",
                      "supernatant" : "sn"}
suffix = "myelin"
odir = "results/plots/pca/"
tps = ["T0h", "T4h", "T24h"]
# todo : how to define in yaml file the needed stuff for doing

domytest = False
if domytest == True:
    for k in compartments_names.keys():
        co = compartments_names[k]
        df = pd.read_csv(f"{tmpdir}meanEnrich_{suffix}_{co}.tsv", sep='\t', header = 0, index_col = 0)
        metadatasub = metadata.loc[metadata['short_comp'] == co,:]
        df = df[metadatasub['sample']]
        # reduce by std each metabolite :
        #df = df + 1
        df = df.div(df.std(axis=1, ddof=0), axis=0) # reduce rows
        #df = df.div(df.sum(axis = 0 )) # normalize columns
        df = df.dropna(axis = 0 )

        calcPCAand2Dplot(df, metadatasub, col1 = "timepoint" , col2="condition",
                         pointlabels="",
                         title = f'{suffix}-{k}', odir=odir)

    print("now parsing time-points to yield separate time-point dataframe for plot")

    for k in compartments_names.keys():
        co = compartments_names[k]
        for tp in tps:
            df = pd.read_csv(f"{tmpdir}meanEnrich_{suffix}_{co}.tsv", sep='\t', header = 0, index_col = 0)
            metadatasub = metadata.loc[(metadata['short_comp'] == co) & (metadata['timepoint'] == tp),:]
            df = df[metadatasub['sample']]
            print(df.shape)
            print(metadatasub.shape)
            # reduce by std each metabolite :
            #df = df + 1
            df = df.div(df.std(axis=1, ddof=0), axis=0) # reduce rows
            #df = df.div(df.sum(axis = 0 )) # normalize columns
            df = df.dropna(axis = 0 )

            calcPCAand2Dplot(df, metadatasub, col1 = "condition" , col2="condition",
                             pointlabels = "sample_descrip",
                             title = f'{suffix}-{k}-{tp}', odir=odir)

iris = sns.load_dataset("iris")
sns.relplot(data=iris, x = "sepal_width", y="petal_width",
             hue="species")
iris = iris.assign(sample = [str(i) for i in iris.index])
fakemeta = iris[["sample", "species"]]
iris = iris.drop(columns=["sample", "species"])
fakedf = iris.T  # variables rows, samples columns
fakedf = fakedf.div(fakedf.std(axis=1, ddof=0), axis=0)

fakedf.columns = [str(i) for i in iris.index]
calcPCAand2Dplot(fakedf, fakemeta, "species", "species", "",
                 "dontcare", "/home/johanna/myelin_metabo/", ("species"))




## cell compartment
# oho = pd.read_csv("tmp/rawAbundances_myelin_cell.tsv", sep='\t', header = 0, index_col = 0)
# #oho = oho + 1
# oho = oho.div(oho.sd(axis=1))
# #oho = oho.div(oho.sum(axis=0))
#
#
# metadatasub = metadata.loc[metadata['short_comp' ] == "cell",:]
#
# calcPCAand2Dplot(oho, metadatasub, col1 = "timepoint" , col2="condition", title = "cells")
#
# ## supernantant compartment
#
# super = pd.read_csv("tmp/rawAbundances_myelin_sn.tsv", sep='\t', header = 0, index_col = 0)
# print(super.head())
# super = super.div(super.sd(axis = 1))
# #super = super.div(super.sum(axis=0))
# metadata = pd.read_csv("data/metadata_myelin.csv", header = 0)
# metadatasub = metadata.loc[metadata['short_comp' ] == "sn",:]
# print(metadatasub.head())

#calcPCAand2Dplot(super, metadatasub, col1 = "timepoint" , col2="condition", title = "supernatant")
