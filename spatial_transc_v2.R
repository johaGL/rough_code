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
library(tidyverse)

# in the files : 
# C stands for Cortex, 
# T stands for tumor, 
# TC stands for Tumor Core, 
# TI stands for Tumor Infiltration.
# - 
# ST: stands for Spatial Transcriptomics.
mywdir = "~/"
setwd(mywdir)
wdadir = "visium_coauthorsBMC/10XVisium_2/10XVisium_2/"
ODIR="~/spatial_hei/"
odirqc = paste0(ODIR,"QC/")
odirclus = paste0(ODIR,"clus/")
dir.create(odirqc)
dir.create(odirclus)
#################
# better order patients by tissue available:
################
#TODO 

prep_seu <- function(thedir, slicename, patient, tissue){
  sobj <- Seurat::Load10X_Spatial(
    data.dir = thedir,
    filename = "filtered_feature_bc_matrix.h5",
    slice = slicename,
    assay="Spatial",
    filter.matrix = TRUE
  )
  #metadata
  sobj@meta.data$patient = patient
  sobj@meta.data$tissue = tissue
  return(sobj)
}


QC_figures_function <- function(rds_file, odirqc){
  # thanks to D Challopin
  
  spseu <- readRDS(rds_file)
  
  name_seu_object = paste0(spseu@meta.data$patient[1],spseu@meta.data$tissue[1])
  
  pdf(paste0(odirqc, "plot_nCount_",name_seu_object,".pdf"),width = 20, height = 10)
  plot1 <- VlnPlot(spseu, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(spseu, features = "nCount_Spatial") + theme(legend.position = "right")
  print(wrap_plots(plot1, plot2))
  dev.off()
  
  pdf(paste0(odirqc, "plot_nFeature_",name_seu_object,".pdf"),width = 20, height = 10)
  plot1 <- VlnPlot(spseu, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(spseu, features = "nFeature_Spatial") + theme(legend.position = "right")
  print(wrap_plots(plot1, plot2))
  dev.off()
  
  # Compare normalization methods
  
  # rerun normalization to store sctransform residuals for all genes
  spseu <- SCTransform(spseu, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
  # also run standard log normalization for comparison
  spseu <- NormalizeData(spseu, verbose = FALSE, assay = "Spatial")
  
  # Computes the correlation of the log normalized data and sctransform residuals with the
  # number of UMIs
  spseu <- GroupCorrelation(spseu, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
  spseu <- GroupCorrelation(spseu, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)
  
  pdf(paste0(odirqc, "plot_normalization_comparison_",name_seu_object,".pdf"),width = 20, height = 10)
  p1 <- GroupCorrelationPlot(spseu, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
    theme(plot.title = element_text(hjust = 0.5))
  p2 <- GroupCorrelationPlot(spseu, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
    theme(plot.title = element_text(hjust = 0.5))
  print(p1 + p2)
  dev.off()
  
  # Gene expression visualization
  pdf(paste0(odirqc,"spatialFeaturePlot_",name_seu_object,".pdf"),width = 10, height = 10)
  print(SpatialFeaturePlot(spseu, features = c( "GFAP", "OSP", "TUBB3")))
  dev.off()
  
  return(spseu)
}



clustering_function <- function(rds_file, odirclus){
  # thanks to D Challopin
  
  spseu <- readRDS(rds_file, odirclus)
  
  name_seu_object = paste0(spseu@meta.data$patient[1],spseu@meta.data$tissue[1])
  
  spseu <- SCTransform(spseu, assay = "Spatial", verbose = FALSE)
  
  spseu <- RunPCA(spseu, assay = "SCT", verbose = FALSE)
  pct <- spseu[["pca"]]@stdev / sum(spseu[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

  spseu <- FindNeighbors(spseu, reduction = "pca", dims = 1:co2)
  spseu <- FindClusters(spseu, verbose = FALSE)
  spseu <- RunUMAP(spseu, reduction = "pca", dims = 1:co2)
  
  pdf(paste0(odirclus,"umap_spatial_",name_seu_object,".pdf"),width = 20, height = 20)
  p1 <- DimPlot(spseu, reduction = "umap", label = TRUE)
  p2 <- SpatialDimPlot(spseu, label = TRUE, label.size = 3)
  p3 <- SpatialFeaturePlot(spseu, features = "nCount_Spatial") + theme(legend.position = "right")
  p4 <- SpatialFeaturePlot(spseu, features = "nFeature_Spatial") + theme(legend.position = "right")
  print(p1 + p2 + p3 + p4)
  dev.off()
  
  return(spseu)
}


# check expression of some markers
FeaturePlot_function <- function(rds_file){
  
  spseu <- readRDS(rds_file)
  
  name_seu_object = paste0(spseu@meta.data$patient[1],spseu@meta.data$tissue[1])
  
  list_markers = c("GFAP", "OLIG1", "TUBB3")
  
  pdf(paste0("FeaturePlots_markers_",name_seu_object,".pdf"),width = 20, height = 15)
  print(SpatialFeaturePlot(spseu, features = list_markers))
  dev.off()
  
  return()
}


## start


# spatially variable features ? ==> Todo?


###################
# several patients : 7
###################                
patients_li = list("UKF242" = c("C","T"), "UKF248"=c("C", "T"),
                     "UKF256"= c("C", "TC", "TI") ,
                   "UKF259"= c("C","T"),
                   "UKF265"=c("C","T"), "UKF313"=c("C","T"),  "UKF334"=c("C","T")
                   )    


for (px in names(patients_li)){
  tissues = patients_li[[px]]
  for (tiss in tissues){
    
   
    px_tiss = paste0(px,"_", tiss)
    
    pxdir = paste0("#",px_tiss,"_ST/")

    sobj = prep_seu(thedir=paste0(wdadir, pxdir, "outs/"),
                    slicename=str_replace(px_tiss, "_", ""),
                    patient=px, tissue=tiss)
    
    saveRDS(sobj, paste0(ODIR, px_tiss))
    
    #print(Sys.glob(paste0(wdadir, "/*/" )))
    
    
    sobj = QC_figures_function(paste0(ODIR, px_tiss), odirqc)
    saveRDS(sobj, paste0(ODIR,px_tiss,"_qc"))
    
    
    
    sobj = clustering_function(paste0(ODIR, px_tiss, "_qc") , odirclus)
    saveRDS(sobj, paste0(ODIR,px_tiss,"_clus"))
    
    
    FeaturePlot_function(paste0(ODIR,px_tiss,"_clus"))
    
  }
}



# other packages I would see in future projects: 
# vissE : think about this pathway enrichment in ST
# https://www.biorxiv.org/content/10.1101/2022.03.06.483195v1.full
