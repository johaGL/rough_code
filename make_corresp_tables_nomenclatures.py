"""
there are discrepancies among 3 different nomenclatures
that I have observed across the targeted metabolomics topic:
- the one of VIB: very clean, follow 'Compound Name' international standard, for most metabolites
- the one of MetaToul : abreviatures are not folowing nor KEGG, nor HMDB ID nor Compound Name standards
- the one of our published figure: abreviations are different to those of MetaToul,
        but the required standards are not met either

Note 1: I am using the compound name as in Human Metabolic Atlas HMA
(caution: for some metabolites, HMA compound name is different to MetaboAnalyst one !!!)
Note 2 : here I am not dealing with VIB nomenclature, only metatoul and HMA ones.

@johaGL
"""
import yaml
import pandas as pd
# take the yaml that I created, and put all metabolites in a table
# with 3 columns : 'Compound_Name' , 'metatoul_name', 'figure_name'

# copied example of output when using metatoul abundances table
aMetaboDiff_file = "cell_TOTAL_AB_T48_Cont_T48_Tt.tsv"
metatoulDiffAbAnalyOut = pd.read_csv(aMetaboDiff_file,sep= '\t')

metatoulnames = metatoulDiffAbAnalyOut['mets']

# Attention: this part is hard-coded, correspondance among metatoul and 'Compound name' official:
dicometatoul_hma = {
    "Fru6P" : "fructose-6-phosphate",
    "Gly3P" : "sn-glycerol-3-phosphate",
    "Fum" : "fumarate",
    "P5P" : "ribose-5-phosphate",
    "Pyr" : "pyruvate",
    "Rib1P" : "ribose-1-phosphate",
    "Suc" : "succinate",
    "ATP" : "ATP",
    "PRPP" : "PRPP",
    "AMP" : "AMP",
    "Cit" : "citrate",
    "a-KG" : "AKG",
    "ADP" : "ADP",
    "Glc6P" : "glucose-6-phosphate",
    "Glycine" : "glycine",
    "Fru1P" : "fructose-1-phosphate",
    "Serine" : "serine",
    "6-PG" : "glucono-1,5-lactone-6-phosphate",
    "FruBP" : "fructose-1,6-bisphosphate",
    "Glutamine" : "glutamine",
    "PEP" : "PEP",
    "Alanine" : "alanine",
    "2_3-PG" : "2,3-bisphospho-D-glycerate",
    "Aspartate" : "aspartate",
    "GTP" : "GTP",
    "IsoCit" : "isocitrate",
    "Glutamate" : "glutamate",
    "Mal" : "malate",
    "2-OHGlu" : "2-hydroxyglutarate",
    "Sed7P" : "sedoheptulose-7-phosphate"
                }


mtouldf = pd.Series(dicometatoul_hma, name = 'compound_name')
mtouldf.index.name = "metatoul_name"
mtouldf = mtouldf.reset_index()
print(mtouldf.head())



jofile = "net_figure.yaml"
with open(jofile, 'r') as f:
    jorisfig = yaml.load(f, Loader=yaml.FullLoader)

dicofig_hma = dict()
for i in jorisfig.keys():
    for j in jorisfig[i].keys():
        tuli = j.split("!")
        hmaname = tuli[0]
        figname = tuli[1]
        dicofig_hma[figname] = hmaname

figdf = pd.Series(dicofig_hma, name = "compound_name")
figdf.index.name = "figure_name"
figdf = figdf.reset_index()
print(figdf.head())

finaldf = mtouldf.merge(figdf, on="compound_name", how="left")
outfile = "nomenclatures_metatoul.csv"
print("saving", outfile)
finaldf.to_csv(outfile, header = True)
