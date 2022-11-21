# note: This script MUST use conda environment called 'scst' (single cell and spatial transcriptome)
# from shell run:
# $ conda activate scst 
# $ rstudio &
# TODO : export new scst.yml because other many pks insts: conda install -c bioconda bioconductor-cardinal and conda-forge r-hdf5r
# TODO continuation : r-irlba r-MALDIquant r-MALDIquantForeign r-viridis
 
library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

# in the files : 
# C stands for Cortex, 
# T stands for tumor, 
# TC stands for Tumor Core, 
# TI stands for Tumor Infiltration.
# - 
# ST: stands for to Spatial Transcriptomics.
mywdir = "~/"
setwd(mywdir)
wdadir = "visium_coauthorsBMC/10XVisium_2/10XVisium_2/"
pxdir = "#UKF241_C_ST/"
px_st = list()
px_st[["one"]] = list()
px_st[["one"]]$thedir = paste0(wdadir, pxdir, "outs/")
#px_st[["one"]]$theh5 = paste0(wdadir, pxdir, "outs/", "filtered_feature_bc_matrix.h5")
px_st[["one"]]$theh5 = "filtered_feature_bc_matrix.h5"
px_st[["one"]]$slice = "detected_tissue_image"
# raw_feature_bc_matrix.h5
# molecule_info.h5
# metrics_summary.csv

# images:
"outs/spatial/" # aligned_fiducials.jpg      scalefactors_json.json  tissue_lowres_image.png
# detected_tissue_image.jpg  tissue_hires_image.png  tissue_positions_list.csv

print(px_st[["one"]]$theh5 )
print(file.exists(px_st[["one"]]$theh5 ))

list.dirs("visium_coauthorsBMC")
print(px_st)
px_one = Seurat::Load10X_Spatial(
  data.dir = px_st[["one"]]$thedir,
  filename = px_st[["one"]]$theh5,
  slice = px_st[["one"]]$slice ,
  assay="one",
  filter.matrix = TRUE, 
  to.upper = FALSE
  
)

# note : we "could" Seurat::merge patients by C T TC TI. But in first time do not
px_one <-  PercentageFeatureSet(px_one, "^MT-", col.name = "percent_mito")
px_one <- PercentageFeatureSet(px_one, "^HB.*-", col.name = "percent_hb")
VlnPlot(px_one, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito",
                            "percent_hb"), pt.size = 0.1, ncol = 2) + NoLegend()
SpatialFeaturePlot(px_one, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito",
                                       "percent_hb"))
SpatialFeaturePlot(px_one, features = c("GFAP", "TUBB3", "LDHA"))



get_genes_inseurobj <- function(seurobj){
  rownames(px_one@assays[[1]])
}



#######################################################
################ end spatial transcripto
# use cardinal to visualize spatial metabolomics : 
# and use also https://paoloinglese.github.io/tut_maldiquant_msi/
library(Cardinal)
library(MALDIquant)
library(MALDIquantForeign)
library(irlba)
library(viridis)


f = "visium_coauthorsBMC/doi_10.5061_dryad.h70rxwdmj__v11/MALDI_1/MALDI_1/raw/20201029scilslab_ncfr_glia_combine_root_mean_square"
print(file.exists(paste0(f,".imzML")))
smetabo = Cardinal::readImzML(file=paste0(f, ".imzML"), name=f) # demands ibd file! 



# note that authors merged the peaks and metadata in the same imzML

folde="visium_coauthorsBMC/doi_10.5061_dryad.h70rxwdmj__v11/MALDI_1/MALDI_1/raw/"
datamz = readImzML(folder=folde, name="20201029scilslab_ncfr_glia_combine_root_mean_square")


pek = MALDIquantForeign::importImzMl(path=paste0(rawmef, ".imzML"),
                                centroided=TRUE, verbose=FALSE)
                 
                 
                 
     

# other packages I would see in future projects: 
# vissE : think about this pathway enrichment in ST
# https://www.biorxiv.org/content/10.1101/2022.03.06.483195v1.full
