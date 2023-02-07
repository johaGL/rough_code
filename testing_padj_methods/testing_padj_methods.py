import numpy

import pandas as pd
import numpy as np
import statsmodels.stats.multitest as smm
import statsmodels

import matplotlib.pyplot as plt

import seaborn as sns
from scipy import stats


print("version of statsmodels in this conda env (dimet)" ,
      statsmodels.__version__) # conda most stable version of statsmodels is not the last one

"""
Created these functions to calculate Benjamini-Hochberg 
 without statsmodels package
"""
def compute_ranks(values: 'numpy.array'):
    original_df = pd.DataFrame({'pvalue': values})
    uniq_df = pd.DataFrame({'pvalue': values})
    uniq_df = uniq_df.drop_duplicates()
    uniq_df = uniq_df.sort_values(by="pvalue", ascending=True)
    uniq_df['rank'] = range(1,uniq_df.shape[0]+1)
    original_df = pd.merge(original_df, uniq_df, on='pvalue',  how='left')
    return np.array(original_df['rank'])

def give_rank_column(df):
    pvalues = np.array(df['pvalue'])
    df['rank'] = compute_ranks(pvalues)
    return df

def compute_pvalue_version3(df):
    copydf = df.copy()
    result_padjs = list()
    for i, row in copydf.iterrows():
        padj_res = (row['pvalue'] * df.shape[0]) / row['rank']
        if padj_res > 1:
            padj_res = 1
        result_padjs.append(padj_res)
    copydf['padj'] = result_padjs
    return copydf


print('small test with fake pvalues')
pvalues = [0.001,0.02,0.045,0.1]
df = pd.DataFrame({'pvalue': pvalues})
df.index = ['a', 'b', 'c', 'd']

df = give_rank_column(df)
df2 = compute_pvalue_version3(df)
oh = 1

print("test several padj calculations, using statsmodels")
tsv = "/home/johanna/padj_addmethod.tsv"


df = pd.read_csv(tsv,sep='\t', header=0, index_col=None)
print(df.head())
x = [4000,4000,4000]
y = [10585,11000,10900]

yes = stats.kruskal(x, y)
print(yes)


sns.histplot(df, x="pvalue")
plt.show()

df = df.sort_values(by="pvalue")

order = np.array(df['pvalue']).argsort()

for method in ['bonferroni', 'sidak',
               'holm-sidak',  'holm', 'simes-hochberg',
               'hommel',  'fdr_bh', 'fdr_by',
                'fdr_tsbh', 'fdr_tsbky']:
    print(method)
    (sgs, corrP, _, _) = smm.multipletests(df["pvalue"], alpha=float(0.05),
                                           method=method)
    df[method] = corrP
print("adding function made without statsmodels package")
df = give_rank_column(df)
dfcopy = compute_pvalue_version3(df)
dfcopy['NOstatsmodels BH'] = dfcopy['padj']
df = pd.merge(df, dfcopy['NOstatsmodels BH'], left_index=True, right_index=True, how='left')

df.to_csv(tsv,sep='\t' )




