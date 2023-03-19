import scipy
import scipy.stats as stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

geom_meansA = np.array([94,26,58,75,31,42,63,87,62,59,37,68,20,67,29,43,61,55,77,94,14])
print(len(geom_meansA))
geom_meansB = 0
ratios1 = np.array([ 0.000002, 0.002, 0.1, 0.95, 1.02, 0.6, 0.99, 0.94, 1,1.15, 1.22,
              1.4, 1.35, 1.05, 1.09, 0.99, 0.04, 1, 3, 5, 7])

A = {'met1':[7,9, 2], 'met2':[9,np.nan,2], 'met3':[100,50,69]}
B = {'met1':[9,12, 5], 'met2':[11,15,16], 'met3':[10,np.nan,2]}

mets = list(A.keys())
manyratios = list()

def combis_separated_groups(arr1, arr2, grsize):
    ratiosall = []
    res_obj1 = itertools.combinations(arr1,grsize)
    gr1 = np.array(list(res_obj1))
    res_obj2 = itertools.combinations(arr2, grsize)
    gr2 = np.array(list(res_obj2))
    for tup1 in gr1:
        for tup2 in gr2:
            ratiosall.append(
               np.nanmean(np.array(tup1)) / np.nanmean(np.array(tup2))
            )
    return(ratiosall)


#for m in mets:
m = 'met1'
maxgroupzize = max(len(A[m]), len(B[m]))
mingroupsize_admitted = 2
for grsize in range(mingroupsize_admitted, maxgroupzize+1):
    reso = combis_separated_groups(A[m], B[m], grsize)


def combis_asone_group(arr_united, grsize) -> list:
    ratiosall = []
    res_obj = itertools.combinations(arr_united, grsize)
    gr = np.array(list(res_obj))
    i = 0
    while i < len(gr):
        tup1 = gr[i]
        tmp_gr = np.delete(gr, i, axis=0)
        for tup2 in tmp_gr:
            ratiosall.append(
               np.nanmean(np.array(tup1)) / np.nanmean(np.array(tup2))
            )
        i += 1
    out_ratios = [i for i in ratiosall if i is not np.nan]
    return out_ratios


def scale_to_interval(array, targ_max, targ_min):
    mi = array.min()
    ma = array.max()
    i01 = ( array - mi ) / (ma - mi )
    outarr = (i01 * (targ_max - targ_min) ) + targ_min
    return outarr


m = 'met1'
resofinal = []
for m in mets:
    maxgroupzize = max(len(A[m]), len(B[m]))
    mingroupsize_admitted = 2
    uniarra = np.concatenate((np.array(A[m]), np.array(B[m])))
    for grsize in range(mingroupsize_admitted, maxgroupzize):
        reso = combis_asone_group( uniarra , grsize)
        resofinal += reso

print(min(resofinal))
# plt.hist(resofinal)
# plt.title(f"all 'ratios' resulting from a total of {len(resofinal)} permutations")
# plt.annotate(f"min={min(resofinal)}\nmax={max(resofinal)}", (max(resofinal),100))
#
#

#plt.annotate(f"min={min(resofinal)}\nmax={max(resofinal)}", (10,10))
myratios = []
for m in mets:
    mydot = np.nanmean(np.array(A[m])) / np.nanmean(np.array(B[m]))
    myratios.append(mydot)
    # plt.vlines(mydot, 0, 100, colors='pink')
# plt.show()

print(min(myratios))
print(min(resofinal))
print(max(myratios))
print(max(resofinal))

scaleddots = scale_to_interval(np.array(myratios), 1,0)
resofinal_0_1 = scale_to_interval(np.array(resofinal), 1, 0)
# plt.hist(resofinal_0_1)
# for i in scaleddots:
#     plt.vlines(i, 0, 10, colors='pink')
# plt.title(f"scaled (0 to 1), n={len(resofinal)} ")
# plt.show()

# r_df = pd.DataFrame({'res_perm': resofinal})
# r_df['scaled'] = scale_to_interval(np.array(resofinal), 1, 0)
# r_df['yeah'] = "yeah"
#
# print(r_df['scaled'].mean())
# sns.boxplot(data=r_df, x='yeah' ,y='scaled')
# plt.show()
# print("**")
# print(r_df.loc[r_df['scaled'] <= 0.05, :])
# print(np.array(resofinal).mean(), ">----")
# meem = np.array(resofinal).mean()
# print(r_df.loc[r_df['res_perm'] < 1, :])

r_df = pd.DataFrame({'res_perm': resofinal})
print(r_df['res_perm'].to_numpy().mean(), " !!")
print(r_df['res_perm'].median())
# flip right tail to left, "center" is ratio = 1

def perm_results_scaling(df):
    df['inum'] = df.index
    df = df.assign(single_tail_disposal= df['res_perm'])
    df['rel_to_1_origi_locati'] = 'right'
    df = df.sort_values('res_perm')
    tmp_inf = df.loc[df['res_perm'] <= 1, :]
    inf_mirror = (tmp_inf['res_perm'].max() - tmp_inf['res_perm'].to_numpy())
    inf_scal = scale_to_interval(inf_mirror, targ_max=max(df['res_perm']),
                                 targ_min=1)
    tmp_inf = tmp_inf.assign(single_tail_disposal=list(inf_scal))

    for i in tmp_inf['inum'].tolist():
        df.loc[df['inum'] == i, 'single_tail_disposal'] = tmp_inf.loc[tmp_inf['inum'] == i, 'single_tail_disposal']
        df.loc[df['inum'] == i, 'rel_to_1_origi_locati'] = 'left'

    df = df.assign(scaled=scale_to_interval(df['single_tail_disposal'].to_numpy(),
                                            targ_max=1,
                                            targ_min=0))

    df = df.sort_values('inum')
    df['pvalues'] = 1 - df['scaled']
    return df

print(r_df)

r_df = perm_results_scaling(r_df)
plt.hist(r_df['scaled'].tolist())
plt.show()
sns.histplot(r_df, x="scaled", hue="rel_to_1_origi_locati", alpha=0.5)
plt.show()
#k = 9
print()
cool = r_df.loc[r_df['pvalues'] <= 0.05,:]
print(cool.sort_values('pvalues'))