#!/usr/bin/env Rscript
#############################################################
#   examples :                                              #
#  $ mywcdir=/scratch/CBIB/jgalvis/recast_scRNApub/saunders #      
#  $ Rscript --vanilla 2_doseurat.R $mycwdir                #  
#                                                           #
#  $ mywcdir=~/recast_scRNApub/saunders                     #      
#  $ Rscript --vanilla 2_doseurat.R $mycwdir                #  
#############################################################
# johaGL

# Note: if it fails,  even with multiprocess and maxSize changed, go to scanpy

library(tidyverse)
library(future)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(sctransform) # not used as mitochondrial reads not present from Sauders

args = commandArgs(trailingOnly=TRUE)
NBsubsau = 2000 # max barcodes by cluster allowed ( ~10000 barcodes "mini" seu).
print(args[1])
setwd(args[1])
if (!is.null(args[2])){
  NBsubsau = args[2] 
}

OUTFILE.frontal <- "FRONTAL_neu-ast-oli_mx.rds"
redofirtsteps <- FALSE
doinitialplots <- FALSE
protofigmarkers <- FALSE
redosubsau <- FALSE


doscale <- function(){
  plan("multicore", workers = 4)
  options(future.globals.maxSize = 10000 * 1024^2)
  future.seed=TRUE
  
  print("loading list (matrix, annots)")
  frontal_l <- readRDS(OUTFILE.frontal)
  frontal.mx <- frontal_l$matrix
  frontal.ann <-frontal_l$annots
  rm(frontal_l) ; gc()  # clear
  print("seurat obj")
  sau <- CreateSeuratObject(frontal.mx, 
                            project="Saunders", 
                            min.cells=3, min.features=200)
  
  colsformd <- frontal.ann[match(rownames(sau@meta.data), 
                                 frontal.ann$bc) , ] 
  rm(frontal.mx, frontal.ann); gc()  # clear
  # as colsformd respects order in sau@meta.data, directly add:
  sau@meta.data$cluster.said <- colsformd$cluster
  sau@meta.data$class <- colsformd$class
  sau@meta.data$full_name <- colsformd$full_name
  sau@meta.data$bc <- colsformd$bc
  sau@meta.data$orig.ident <- "FC"
  sau <- SetIdent(sau, value = sau@meta.data$orig.ident)
  sau[["percent.mt"]] <- PercentageFeatureSet(sau, pattern="^Mt-")
  print("preparing plots")
  # #mitochondrial features have been already removed by Saunders team
  # VlnPlot(sau, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
  # a <- FeatureScatter(sau, feature1="nCount_RNA", feature2="nFeature_RNA") +
  #   labs(title="before subset") + theme_bw()
  # sau <- subset(sau, subset= nFeature_RNA >200 & nFeature_RNA < 6000)
  # b <- FeatureScatter(sau, feature1="nCount_RNA", feature2="nFeature_RNA") +
  #   labs(title="after subset") + theme_bw()
  # c <- VlnPlot(sau, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
  #              ncol=3, pt.size = .05, cols="cadetblue") 
  
  # pdf("plots_feat_counts.v2.png")
  # ( ( a | b ) / c ) + plot_layout(heights = c(3,2), widths=c(1,1))
  # dev.off()
  print("normalize")
  # normalize :  CPM, see https://satijalab.org/seurat/reference/normalizedata
  sau <- NormalizeData(sau, 
                       normalization.method="RC",  
                       scale.factor = 1e6) 
  
  allgenes <- rownames(sau)
  sau <- FindVariableFeatures(sau, selection.method = "vst", nfeatures = 2000)
  sau <- ScaleData(sau, features=allgenes)
  saveRDS(sau,"norm_scal.rds")
}

do.pca.plots.neigh.clus <- function(){
  plan("multicore", workers = 4)
  options(future.globals.maxSize = 10000 * 1024^2)
  future.seed=TRUE
  
  print("loading seurat object ...")
  sau <- readRDS("norm_scal.rds")
  print("run PCA")
  sau <- RunPCA(sau, features=VariableFeatures(object = sau))
  print("doplots")
  pcaplo <- DimPlot(sau, reduction="pca") 
  elbowpl <- ElbowPlot(sau)
  png("plots/sau_pca.png", width=800); (pcaplo | elbowpl) ;dev.off()
  
  pdf("plots/sau_dims.pdf")
  DimHeatmap(sau, dims = 1:12, cells = 500, balanced = TRUE, combine=TRUE)
  dev.off()
  
  print("run neigh")
  sau <- FindNeighbors(sau, dims = 1:5)
  print("run find clus")
  sau <- FindClusters(sau, resolution=0.6)
  saveRDS(sau,"normscalpca_clus.rds")
}

do.tsne.umap <- function(){
  print("loading seurat object")
  sau <- readRDS("normscalpca_clus.rds")
  print("running tsne and umap")
  sau <- RunTSNE(sau, dims=1:10)
  sau <- RunUMAP(sau, dims=1:10)
  saveRDS(sau, "sau-tsne-umap.rds")
}

massageplots <- function(alistofplots){
  outli <- list()
  alegend <- cowplot::get_legend(alistofplots[[1]])
  l <- 1; for (k in alistofplots){
    outli[[l]] <- k + theme_bw() + theme(legend.position = "none")
    l <- l + 1 
  }
  outli[[l]] <- NULL
  outli[[l+1]] <- alegend
  return(outli)
}

# plots and more : 
dotheplots <- function(sau, diu, MYCAT){
  RESCLUSTERS = max(as.numeric(levels(sau@meta.data$seurat_clusters))) + 1
  colall <- colorRampPalette(brewer.pal(n=12,"Set3"))(RESCLUSTERS)
  #names(colall) <- seq(0,20)
  
  print(DimPlot(sau, reduction=diu, label=TRUE) + theme_bw() +
          scale_color_manual(values=colall) + theme(legend.position="none") )
  print( DimPlot(sau, reduction=diu, group.by=MYCAT, 
                 label = TRUE,
                 pt.size=0.001) +  theme_bw() +
           scale_color_brewer(palette="Dark2")  + 
           theme(legend.position="bottom") )
  combi <- FeaturePlot(sau, 
                       features=c("Rbfox2","Tubb3", "Gad1", "Gad2", "Reln", "Syp",
                                  "Gfap","Aldh1l1",  "Opalin", "Flt1", "Cx3cr1", "C1qc",  "Itgam" ),
                       order=TRUE, # order true brings high expressions front
                       keep.scale = "all", 
                       reduction = diu,
                       combine = FALSE, 
                       cols=c("gray95", "red4")) 
  print(cowplot::plot_grid(plotlist=massageplots(combi)))
  sau@meta.data$acla <- sau@meta.data$class
}


# # ------------ START: step by step : 
gc()
if (redofirtsteps){
  doscale()
  do.pca.plots.neigh.clus()
  do.tsne.umap()
}
# -----------------------------------------------------------------------
gc()
if (doinitialplots){
  sau <- readRDS("sau-tsne-umap.rds")
  gc()
  dotheplots(sau, "tsne", "class")
  dotheplots(sau, "umap", "class")
} 

# too many undefined clusters at the first plot, clear from small undefined ones: 
# TO KEEP : 
# 6,2, 17, 9  : Tubb3
# 3, 1, 8 :  endothelial
# 13,5,4,14 : oligo
# 10, 0 : astro  (exclude 7 because + both Tubb3 and Gfap)
# 16 : mural : include 
# 11 : microglia
# 
# in whole : 1,2,3,4,5,6,8,9,10,11,13,14, 16, 17
if (!file.exists("sau_tsne_umap_whole.rds")){
  tokeepclus <- c(0,1,2,3,4,5,6,8,9,10,11,13,14, 16, 17)
  sau.sel <- subset(x = sau, idents = tokeepclus)

  newids_ <- list("Neurons"=c(6,2,17,9),
                "Endothelial"=c(3,1,8),
                "Oligodendrocytes"=c(13,5,4,14),
                "Astrocytes"=c(10,0),
                "Microglia"=c(11),
                "Mural"=c(16))
  tmpdf1 <- data.frame(clus_calc=tokeepclus)  
  tmpdf1$newclass <- NA
  for (k in names(newids_)){
    tmpdf1[tmpdf1$clus_calc %in% newids_[[k]],]$newclass <- k
  }  
  
  tmpdf2 <- sau.sel@meta.data
  tmpdf2$celltype <- tmpdf1[match(tmpdf2$seurat_clusters, tmpdf1$clus_calc),]$newclass
  sau.sel@meta.data <- tmpdf2
  dotheplots(sau.sel, "tsne", "celltype")
  gc()
  saveRDS(sau.sel, "sau_tsne_umap_whole.rds")
}

# -- chunk markers heatmap in full-----------------------------
if (protofigmarkers){
  sau.sel <- readRDS("sau_tsne_umap_whole.rds")
  Idents(sau.sel) <- sau.sel@meta.data$celltype
  allmarks <- FindAllMarkers(sau.sel, only.pos=TRUE, min.diff.pct = 0.1)
  saveRDS(allmarks, "allmarkers_whole.rds")
  # cluster has changed from numbers to  celltype (== cluster).
  topgene.a <- allmarks %>% 
    filter(cluster %in% c("Endothelial","Microglia","Mural", "Oligodendrocytes")) %>% 
    group_by(cluster) %>%
    slice_max(order_by=avg_log2FC, n=5)
  topgene.b <-  allmarks %>% filter(cluster == "Astrocytes") %>% 
    group_by(cluster) %>%
    slice_max(order_by=avg_log2FC, n=15)
  topgene.c <-  allmarks %>% filter(cluster == "Neurons") %>% 
    group_by(cluster) %>%
    slice_max(order_by=avg_log2FC, n=15)  # 330 will include Tubb3
  
  topgene.minisau <- rbind(topgene.a, topgene.b, topgene.c)
  toporder <- intersect(allmarks$gene,topgene.minisau$gene)
  
  
  pdf("plots/sau_markers_whole.pdf", height = 22, width=14)
  DoHeatmap(sau.sel, features= toporder,
            size=3, angle=0, hjust="center")  
  dev.off()
  print("Heatmap sau_markers_whole : ggplot failed ! no image generated")
  print("reason: huge amount of barcodes by cluster")
}


### --mini dataset to print the markers heatmap and do deconvolution ----------
if (redosubsau){
  sau.sel <- readRDS( "sau_tsne_umap_whole.rds")
  gc()
  table(sau.sel@meta.data$celltype)
  # minisau : randomly choose, for each celltype, NBsubsau cells
  table(sau.sel@meta.data$celltype)
  krando <- list() ; for (k in unique(sau.sel@meta.data$celltype)){
    print(k)
    tmpk <- sau.sel@meta.data %>% filter(celltype==k)
    if (nrow(tmpk) > NBsubsau){
      krando[[k]] <- sample(tmpk$bc, NBsubsau, replace=FALSE)
    }else{ krando[[k]] <- tmpk$bc } }
  krando.v <- unlist(krando)
  sau.mini <- subset(sau.sel,cells=krando.v)
  rm(sau.sel);gc()
  
  # chunk  save mini rds seurat object -------------------------------
  saveRDS(sau.mini, "sau_tsne_umap_mini.rds")
  
  # chunk heatmap markers with this mini version  ---------------------
  Idents(sau.mini) <- sau.mini@meta.data$celltype
  allma.minisau <- FindAllMarkers(sau.mini, only.pos=TRUE, min.diff.pct = 0.1)
  saveRDS(allma.minisau,  "allmarkers_mini.rds")   ## save
  # cluster has changed from numbers to  celltype (== cluster).
  # pick different number of markers to show in the markers heatmap: 
  topgene.a <- allma.minisau %>% 
    filter(cluster %in% c("Endothelial","Microglia","Mural", "Oligodendrocytes")) %>% 
    group_by(cluster) %>%
    slice_max(order_by=avg_log2FC, n=5)
  topgene.b <-  allma.minisau %>% filter(cluster == "Astrocytes") %>% 
    group_by(cluster) %>%
    slice_max(order_by=avg_log2FC, n=25)
  topgene.c <-  allma.minisau %>% filter(cluster == "Neurons") %>% 
    group_by(cluster) %>%
    slice_max(order_by=avg_log2FC, n=25)  # 330 will include Tubb3
  
  topgene.minisau <- rbind(topgene.a, topgene.b, topgene.c)
  toporder <- intersect(allma.minisau$gene,topgene.minisau$gene)
  pdf("plots/sau_markers_miniV2.pdf", height = 13, width=10) ## write
  DoHeatmap(sau.mini, features= toporder,
            size=3, angle=0, hjust="center") 
  dev.off()
}

###### ------------- end: list of plots --------------------------------
print(list.files(path="plots/"))
cat("END")

## END



# plan("multicore", workers = 4)
# # Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead, explicitly 
# # specify either 'multisession' or 'multicore'. In the current R session, 'multiprocess' equals 'multicore'.
# 
# # Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) : 
# #   The total size of the 15 globals exported 
# # for future expression (‘FUN()’) is 1.77 GiB.. This exceeds the maximum allowed 
# # size of 500.00 MiB (option 'future.globals.maxSize'). The three largest globals 
# # are ‘object’ (1.77 GiB of class ‘S4’), ‘as’ (228.67 KiB of class ‘function’) and
# # ‘.asCoerceMethod’ (80.58 KiB of class ‘function’)
# 
# options(future.globals.maxSize = 10000 * 1024^2)
# # 850 mib example: 850*1024^2; 

## future.seed=TRUE
# Running Louvain algorithm...
# Maximum modularity in 10 random starts: 0.8642
# Number of communities: 19
# Elapsed time: 5 seconds
# Warning message:
#   UNRELIABLE VALUE: One of the ‘future.apply’ iterations (‘future_lapply-1’) 
# unexpectedly generated random numbers without declaring so. There is a risk that those 
# random numbers are not statistically sound and the overall results might be invalid. 
# To fix this, specify 'future.seed=TRUE'. This ensures that proper, parallel-safe random 
# numbers are produced via the L'Ecuyer-CMRG method. To disable this check, use 'future.seed = NULL', or set option 'future.rng.onMisuse' to "ignore". 
