def red_cv_g_bytime_dffull(names_compartments, dirtmpdata, tableAbund,
                       namesuffix, metadata, levelstimepoints_):
    def give_geommeans(df_red, red_meta, t):
        condilist = red_meta["condition"].unique()
        tmpdico = dict()
        for condi in condilist:
            samplesthiscondi = red_meta.loc[red_meta['condition'] == condi, "sample"]
            subdf = df_red[samplesthiscondi]
            subdf = subdf.assign(geomean=subdf.apply(gmean, axis=1))
            # print(subdf.head())
            tmpdico[f"{condi}_{t}_geomean"] = subdf.geomean.tolist()

        dfout = pd.DataFrame.from_dict(tmpdico)
        dfout.index = df_red.index

        return dfout
    def give_coefvar_bycond(df_red, red_meta, t):
        condilist = red_meta["condition"].unique()
        tmpdico = dict()
        for condi in condilist:
            samplesthiscondi = red_meta.loc[red_meta['condition'] == condi, "sample"]
            subdf = df_red[samplesthiscondi]
            subdf = subdf.assign(CV=subdf.apply(compute_cv, axis=1))
            # print(subdf.head())
            tmpdico[f"{condi}_{t}_CV"] = subdf.CV.tolist()

        dfout = pd.DataFrame.from_dict(tmpdico)
        dfout.index = df_red.index
        return dfout

    # calculate and save : reduced data, coefficient of variation, splitting by timepoint, here only T0h test
    ddof = 0
    for co in names_compartments.values():
        df = pd.read_csv(f"{dirtmpdata}{tableAbund}_{namesuffix}_{co}.tsv", sep='\t', header=0, index_col=0)

        metada_sel = metadata.loc[metadata['short_comp']==co, :]
        #get reduced rows , cv and geommeans,
        for t in ['T0h']: # for t in levelstimepoints_
            print(t)
            samples_t = metada_sel.loc[metada_sel['timepoint'] == t, "sample"]
            samples_t = sorted(samples_t)
            df_t = df[samples_t]
            rownames = df_t.index
            df_t.index = range(len(rownames))#  index must be numeric because compute reduction accepted
            df_t_red = compute_reduction(df_t, ddof)
            df_t_red.index = rownames

            #outfilereduced = f"{dirtmpdata}abund_reduced_{t}_{co}.tsv"
            # df_t_red.to_csv(outfilereduced, header=True, sep='\t')
            # add coefficient of variation, by condition
            red_meta = metada_sel.loc[metada_sel["sample"].isin(df_t_red.columns) ,:]

            df_cv = give_coefvar_bycond(df_t_red, red_meta, t )
            df_t_red_cv = pd.merge(df_t_red, df_cv , left_index=True, right_index=True)
            #outfi_coefvar = f"{dirtmpdata}abund_reduced_coefvar_{t}_{co}.tsv"
            #df_t_red_cv.to_csv(outfi_coefvar)

            # save intervals overlap # TODO?... but here we have 3 conditions (in this timepoint)

            # save geometric means table, and the ratio
            df_t_red_geomean = give_geommeans(df_t_red, red_meta, t)
            print(df_t_red_geomean.head())
            dfo1 = pd.merge(df_t_red_cv, df_t_red_geomean, left_index=True, right_index=True)
            outfi_geomean = f"{dirtmpdata}abund_reduced_geomean_{t}_{co}.tsv"
            dfo1.to_csv(outfi_geomean, header=True)
    return 0


# def getspecificmk(prop_df, isotosrowdata, selmk):
#     """
#     obtain proportion of CI for specified label selmk
#     example : selmk = "m+0"
#     output: columns are [samples], rownames metabolites
#     NOTE: isotosrowdata has been obtained from prop_df (%) via yieldrowdataB (fun_fm)
#     """
#     tmp = prop_df.copy()
#
#     tmp = tmp.assign(isotopolgFull=tmp.index)
#     tmp = pd.merge(tmp, isotosrowdata, on="isotopolgFull", how="inner")
#     tmp = tmp.loc[tmp["m+x"] == selmk, :]
#     tmp = tmp.set_index("metabolite")  # make metabolite to be the index
#     tmp = tmp.drop(["isotopolgFull", "m+x"], axis=1)
#     return tmp


# def yieldmarkedabu(prop_mx, abund, *totalmarked):
#     """
#     * calculates abundance using proportion specific m+x dataframe *
#     NOTE: if option totalmarked is 'totalmarked' (not case sentitive):
#     * * calculates by substracting mzero abundance from total abundance * *
#     Parameters
#     ----------
#     prop_mx : pandas
#         example : prop_mx = getspecificmk(prop_df, isotosrowdata, "m+0")
#           so prop_mx :    sampleConditionA-1   sampleConditionA-2 ....
#                    metabolites           0                      0
#
#     abund : pandas
#
#     Returns
#     -------
#     abu_mx : pandas
#     """
#     total_opt = [o for o in totalmarked]
#     # make sure same order in columns
#     ordercol = prop_mx.columns
#     abund = abund[ordercol]
#     markedout = dict()
#     # calculate by rows, across mx marked proportions
#     # note: if prop_zero as input, row will be zero marked proportions
#     for i, row in prop_mx.iterrows():
#         nameindex = row.name  # the metabolite in row
#         totabu = abund.loc[nameindex, :]
#         if len(total_opt) == 0:
#             abumk = (totabu * row) / 100
#         else:
#             if total_opt[0].lower() == "totalmarked":
#                 zeroabu = (totabu * row) / 100  # zero marked abundance,
#                 abumk = totabu - zeroabu  # total - zero marked abundance
#             else:
#                 # wrote string not recognized (total only accepted)
#                 print("3dr arg not recognized, must be 'totalmarked' or empty")
#                 return 1
#         markedout[nameindex] = abumk
#
#     abu_mx = pd.DataFrame(markedout).T
#     return abu_mx


# def callfuns_perc2abu(odir, dicos, co, tableAbund, tableIC, selmk, *totalmarked):
#     """
#     Parameters
#     ----------
#     dicos :
#     co : compartment
#     tableAbund : pandas
#     tableIC : pandas
#     selmk : "m+3" for example
#         DESCRIPTION.
#     *totalmarked : string
#         "totalmarked" (if willing totalMARKED abundance, set selmk = m+0)
#
#     """
#     total_opt = [o for o in totalmarked]
#     abund = dicos[co][tableAbund]
#     prop_df = dicos[co][tableIC]
#     prop_df = correction_prop_df(prop_df)
#     isotosrowdata = yieldrowdataB(prop_df)  # this funciton is from fun_fm
#
#     prop_mx = getspecificmk(prop_df, isotosrowdata, selmk)
#
#     if len(total_opt) == 0:
#         abu_mx = yieldmarkedabu(prop_mx, abund)
#         nameout = f"{odir}abux_byProp_{selmk}_{co}.tsv"
#         saveabundance(abu_mx, nameout)
#     else:
#         if total_opt[0].lower() == "totalmarked":
#             abu_mkall = yieldmarkedabu(prop_mx, abund, "totalmarked")
#             nameout = f"{odir}abux_byProp_totmk_{co}.tsv"
#             saveabundance(abu_mkall, nameout)
#     return 0


# def saveabundfrompercentagesIC(
#     datadi,
#     tableAbund,
#     tableIC,
#     metadata,
#     names_compartments,
#     namesuffix,
#     odir,
#     max_m_species
# ):
#     if not os.path.exists(odir):
#         os.makedirs(odir)
#
#     for co in names_compartments.values():
#         abun = pd.read_csv(
#             datadi + tableAbund + "_" + namesuffix + "_" + co + ".tsv",
#             sep="\t",
#             index_col=0,
#         )
#         # note that pandas automatically transform any 99.9% in decimal 0.999
#
#         propc = pd.read_csv(
#             datadi + tableIC + "_" + namesuffix + "_" + co + ".tsv",
#             sep="\t",
#             index_col=0,
#         )
#         dicos = dict()
#         dicos[co] = {}
#         dicos[co]["metadata"] = metadata.loc[metadata.short_comp == co]
#         selecols = dicos[co]["metadata"]["sample"]
#         dicos[co][tableAbund] = abun[selecols]
#         dicos[co][tableIC] = propc[selecols]
#
#         # by default calculates total marked (taking away m+0 unmarked)
#         callfuns_perc2abu(odir, dicos, co, tableAbund, tableIC, "m+0", "totalmarked")
#         callfuns_perc2abu(odir, dicos, co, tableAbund, tableIC, "m+0")  # true unmarked
#
#         for k in range(1, max_m_species + 1):
#             callfuns_perc2abu(odir, dicos, co, tableAbund, tableIC, f"m+{k}")
#
#     return 0
