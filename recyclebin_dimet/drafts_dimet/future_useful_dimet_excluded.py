
import functions_general as fg
import os
import prepare as prp


def plot_overlap_hist(df_overls, colname_symetric, colname_assymetric, fileout):
    import seaborn as sns
    import matplotlib as plt
    """just for debugging or other tests"""
    values_sym = df_overls[colname_symetric]
    a = pd.DataFrame({'value' : values_sym,
                      'type_overlap': ["symm" for i in range(len(values_sym))] })
    vasym = df_overls[colname_assymetric]
    b = pd.DataFrame({'value': vasym,
                      'type_overlap': ["assym" for i in range(len(vasym))]})
    dfplotov = pd.concat([a,b], ignore_index=True, axis=0)

    with sns.axes_style("darkgrid"):
        sns.displot(data=dfplotov, x = 'value', hue='type_overlap',
                       kde=False)
        plt.savefig(fileout)
    plt.close()
    return 0


# # do pytest in tests/test_prepare, needs pytest and pytest-cov and mypy
def test_excelsheets2frames_dic():
    configfile = "examples/toy1/a001/config-1-001.yml"
    confidic = fg.open_config_file(configfile)
    targetedMetabo_path = "examples/toy1/data/results_toy1.xlsx"
    assert prp.excelsheets2frames_dic(targetedMetabo_path, confidic) is not None


# end pytest 
#######################"


def verify_metadata_coherence_further(metadata_df) -> None:
    #use some elements of this function to 
    # create the function to yield good metadata :)
    # then get rid of this function XD
    txt_errors = ""
    poor_sample_names = list()
    for i, r in metadata_df.iterrows():
        if not ((r['condition'] in r['sample']) and
                (str(r['timenum']) in r['sample']) and
                (str(r['timenum']) in r['timepoint']) and
                (r['timepoint'] in r['sample'])):
            #suggestion =  # or timenum ?
            poor_sample_names.append([r['sample'], r['condition'],
                                       str(r['timenum']), r['timepoint']])

    if len(poor_sample_names) > 0:
        txt_errors += "\n-> discrepancies:"
        txt_errors += "\nsample\tcondition\ttimenum\ttimepoint\n"
        for tup in poor_sample_names:
            txt_errors += '\t'.join(tup)
            txt_errors += '\n'
        warnings.warn("your sample names, to be improved: ")
        
        
def wdir_configpaths_validate(wdir, config) -> None:
    aproval = [True,True, True]
    if not os.path.isdir(wdir):
        print(f"not a directory: {wdir}")
        aproval[0] = False
    if not os.path.isfile(config):
        print(f"not a configuration file: {config}")
        aproval[1] = False
    if not os.path.exists(wdir + "data/"):
        print(f"data/ folder is missing in {wdir}")
    if not (aproval[0] and aproval[1] and aproval[2]):
        raise ValueError("\nDid you inverted the order of directory and config file?")


