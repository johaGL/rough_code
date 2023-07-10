from __future__ import division
import sys

import numpy as np
import pandas as pd
import statsmodels.stats.multitest as smm

def compute_p_adjusted(df: pd.DataFrame, correction_method: str) -> pd.DataFrame:
    rej, pval_corr = smm.multipletests(df['pvalue'].values, alpha=float('0.05'), method=correction_method)[:2]
    df['padj'] = pval_corr
    return df

print("using a dummy array of values")
#kn = [0.04, 0.004, 0.1, 0.5, 0.9]
#kn = [0.24, 0.5, 0.5, 1, 0.0001]
kn = [1,  0.290034529038187, 0.320932951534461, 0.2065564914483, 0.174633900731984]


rej_bool, pval_arr = smm.multipletests(kn, alpha=0.05, method="fdr_bh")[:2]
print(kn)
print(pval_arr)
# or using the function
df = pd.DataFrame({"rownb": [i for i in range(len(kn))],
                   "pvalue" : kn})
df = compute_p_adjusted(df, "fdr_bh")
print(df)


#####
# from https://stackoverflow.com/questions/25185205/calculating-adjusted-p-values-in-python
def fdr(pvals):
    tmp = list()
    from scipy.stats import rankdata
    ranked_p_v = rankdata(pvals)
    print(ranked_p_v)
    #fdr = ( pvals * len(pvals)) / ranked_p_v
    for i in range(len(pvals)):
        tmp.append((pvals[i] * len(pvals)) / ranked_p_v[i])
    fdr = np.array(tmp)
    fdr[fdr > 1] = 1
    return fdr


ehe = fdr(np.array(kn))
print(ehe,  "<======== by hand, 1")

# --
# from : https://rosettacode.org/wiki/P-value_correction#Python

def pminf(array):
    x = 1
    pmin_list = []
    N = len(array)
    for index in range(N):
        if array[index] < x:
            pmin_list.insert(index, array[index])
        else:
            pmin_list.insert(index, x)
    return pmin_list
#end function


def cumminf(array):
    cummin = []
    cumulative_min = array[0]
    for p in array:
        if p < cumulative_min:
            cumulative_min = p
        cummin.append(cumulative_min)
    return cummin
#end

def cummaxf(array):
    cummax = []
    cumulative_max = array[0]
    for e in array:
        if e > cumulative_max:
            cumulative_max = e
        cummax.append(cumulative_max)
    return cummax
#end

def order(*args):
    if len(args) > 1:
        if args[1].lower() == 'false':# if ($string1 eq $string2) {
            return sorted(range(len(args[0])), key = lambda k: args[0][k])
        elif list(args[1].lower()) == list('true'):
            return sorted(range(len(args[0])), key = lambda k: args[0][k], reverse = True)
        else:
            print("bad")
            sys.exit()
    elif len(args) == 1:
        return sorted(range(len(args[0])), key = lambda k: args[0][k])

o = order(kn, 'TRUE')
cummin_input = []
for index in range(len(kn)):
    cummin_input.insert(index, (index + 1) * kn[o[index]])
cummin = cumminf(cummin_input)
pmin = pminf(cummin)
ro = order(o)
qvalues = [pmin[i] for i in ro]

print(qvalues, " <================ by hand, 2")


# --

### ####
print("\nUsing some real values")
prekn = [0.586519014724443 , 0.609126394687426 , 0.490666135722899, 0.425655265804451, \
    0.239121611248683, 0.516810415015153, 1, 0.585822774747661 ]

kn = np.array(prekn).round(6)
print(kn)
rej_bool, pval_arr = smm.multipletests(kn, alpha=0.05, method="fdr_bh")[:2]
print(pval_arr)
print("using 'is_sorted'")
# now function but ordered
def compute_p_adjustedV2(df: pd.DataFrame, correction_method: str, is_sorted: bool) -> pd.DataFrame:
    rej, pval_corr = smm.multipletests(df['pvalue'].values,
                                       alpha=float('0.05'),
                                       method=correction_method,
                                       is_sorted=is_sorted)[:2]
    df['padj'] = pval_corr
    return df
df = pd.DataFrame({"rownb": [i for i in range(len(kn))],
                   "pvalue" : kn})
df2 = df.sort_values(by='pvalue', ascending=True)
df2 = compute_p_adjustedV2(df2, "fdr_bh", is_sorted=True)
print(df)
print("useless is_sorted\n")

print("using fdrcorrection function")

rej, pval_corr = smm.fdrcorrection(df['pvalue'].values, method='p',alpha= 0.05)
print(rej)
print(pval_corr)


## ## ##
print("\nUsing several method in pure multipletests function")
print("pvalues :", df['pvalue'])
for imethod in ['bonferroni', 'sidak', 'homl-sidak', 'holm','hommel',
               'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky']:
    try:
        print(imethod)
        rejtmp, pvaltmp = smm.multipletests(df['pvalue'].values, alpha=float('0.05'), method=imethod)[:2]
        print(pvaltmp)
    except:
        pass
