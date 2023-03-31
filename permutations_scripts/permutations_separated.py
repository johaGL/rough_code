import scipy.stats as stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools


def compute_gmean_nonan(anarray):
    """ do not recopy in project, already in functions general"""
    anarray = np.array(anarray, dtype=float)
    anarray = anarray[~np.isnan(anarray)]
    if sum(anarray) == 0:  # replicates all zero
        outval = 0
    else:
        outval = stats.gmean(anarray)
    return outval


def compute_ratios_p_permutations(arr_united, p_size) -> list:
    """
    output : the ratios of all possible combinations of size r of the array
    """
    ratiosall = []
    res_obj = itertools.combinations(arr_united, p_size)
    gr = np.array(list(res_obj))
    i = 0
    while i < len(gr):
        tup1 = gr[i]
        tmp_gr = np.delete(gr, i, axis=0)
        for tup2 in tmp_gr:
            ratiosall.append(
               compute_gmean_nonan(np.array(tup1)) / compute_gmean_nonan(np.array(tup2))
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


def perm_results_2_pvalues(df, column_name):
    df = df.assign(single_tail_disposal= df[column_name]) # initialize
    df['rel_to_1_origi_locati'] = 'right'
    df = df.sort_values(column_name)
    tmp_inf = df.loc[df[column_name] <= 1, :]
    inf_mirror = (tmp_inf[column_name].max() - tmp_inf[column_name].to_numpy())
    inf_scal = scale_to_interval(inf_mirror, targ_max=max(df[column_name]),
                                 targ_min=1)
    tmp_inf = tmp_inf.assign(single_tail_disposal=list(inf_scal))

    for i in tmp_inf['id_unique'].tolist():
        df.loc[df['id_unique'] == i, 'single_tail_disposal'] = tmp_inf.loc[tmp_inf['id_unique'] == i, 'single_tail_disposal']
        df.loc[df['id_unique'] == i, 'rel_to_1_origi_locati'] = 'left'

    df = df.assign(scaled=scale_to_interval(df['single_tail_disposal'].to_numpy(),
                                            targ_max=1,
                                            targ_min=0))
    df = df.sort_values('id_unique')
    df['pvalues'] = 1 - df['scaled']
    return df


def perm_results_2_pvalues(df, column_name):
    df = df.assign(single_tail_disposal= df[column_name]) # initialize
    df['rel_to_1_origi_locati'] = 'right'
    df = df.sort_values(column_name)
    tmp_inf = df.loc[df[column_name] <= 1, :]
    inf_mirror = (tmp_inf[column_name].max() - tmp_inf[column_name].to_numpy())
    inf_scal = scale_to_interval(inf_mirror, targ_max=max(df[column_name]),
                                 targ_min=1)
    tmp_inf = tmp_inf.assign(single_tail_disposal=list(inf_scal))

    for i in tmp_inf['id_unique'].tolist():
        df.loc[df['id_unique'] == i, 'single_tail_disposal'] = tmp_inf.loc[tmp_inf['id_unique'] == i, 'single_tail_disposal']
        df.loc[df['id_unique'] == i, 'rel_to_1_origi_locati'] = 'left'

    df = df.assign(scaled=scale_to_interval(df['single_tail_disposal'].to_numpy(),
                                            targ_max=1,
                                            targ_min=0))
    df = df.sort_values('id_unique')
    df['pvalues'] = 1 - df['scaled']
    return df

##


input_df = pd.DataFrame.from_dict({'met1': [7,9, 2, 9,12, 5],
                 'met2':[14,np.nan,2, 11,15,16] ,
                 'met3' : [100,50,69, 10,np.nan,2],
                 'met4': [240, 310, 200, 80,88,89]})
input_df = input_df.T
input_df.columns = ['a1','a2','a3', 'b1', 'b2', 'b3']

groupA = ['a1','a2','a3']
groupB = ['b1', 'b2', 'b3']
maxgroupsize = max(len(groupA), len(groupB))
dico_permus_res = dict()
df = input_df.copy()
mingroupsize_admitted = 2

for i, row in df.iterrows():
    tmp_arr = df.loc[i,:].to_numpy()
    for p_size in range(mingroupsize_admitted, maxgroupsize + 1):
        reso = compute_ratios_p_permutations( tmp_arr , p_size)
    dico_permus_res[i] = reso

print(dico_permus_res)
# visualization
# for k in dico_permus_res.keys():
#     ratios_permuts = dico_permus_res[k]
#     r_df = pd.DataFrame({'ratios_permut': ratios_permuts})
#     r_df = r_df.assign(id_unique=["permut%" + str(i) for i in r_df.index])
#
#     print("Obtained ratios from", r_df.shape[0], "permutations. Metabolite:", k)
#     print("mean:", r_df['ratios_permut'].to_numpy().mean())
#     print("median:", r_df['ratios_permut'].median())
#
#     r_df = perm_results_2_pvalues(r_df, 'ratios_permut')
#     plt.hist(r_df['ratios_permut'].tolist(), bins=30)
#     plt.title(f"ratios resulting from {r_df.shape[0]} permutations. {k}")
#     plt.show()
#     myvar_plot = "pvalues"  # "scaled"
#     sns.histplot(r_df, x=myvar_plot, hue="rel_to_1_origi_locati", alpha=0.5, bins=30)
#     plt.title(f"computed {myvar_plot}, from permutations. {k}")
#     plt.xlabel("values")
#     plt.show()

# continue
for k in dico_permus_res.keys():
    ratios_permuts = dico_permus_res[k]
    r_df = pd.DataFrame({'ratios': ratios_permuts})
    r_df = r_df.assign(id_unique=["permut%" + str(i) for i in r_df.index])
    row = input_df.loc[k, :]  # row of the metabolite in observations
    interest = row[groupB]
    baseline = row[groupA]
    print(interest, baseline)
    observed_ratio = compute_gmean_nonan(np.array(interest)) / compute_gmean_nonan(np.array(baseline))
    new_row = pd.DataFrame({'ratios': [observed_ratio], 'id_unique': [k]})
    r_df = pd.concat([r_df, new_row]).reset_index()
    r_df = perm_results_2_pvalues(r_df, 'ratios')
    print(r_df.loc[r_df['id_unique']==k, : ])