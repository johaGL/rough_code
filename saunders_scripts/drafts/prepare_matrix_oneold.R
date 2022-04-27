setwd("~/dataref_scRNA/saunders")
library(tidyverse)
OUTFILEV2.frontal <- "FRONTALV2_neu-ast-oli_mx.rds"
myconda <- "/home/johanna/programs_cmd/mambaforge/bin/conda"
mycondaname <- "scanpy"  # name as told when created the conda env

OUTFILE.frontal <- "FRONTAL_neu-ast-oli_mx.rds"

q <- readRDS("annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
# metacells : 
#w <- readRDS("metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
print(dim(w)) # 32307 genes x 565 metacells
print(dim(q)) # the metadata for the 565 metacells
# metacell is the UMI aggregates for each subcluster, see: http://dropviz.org/#tab-2497-4
table(q$tissue)

#CB ENT  FC  GP  HC  PC  SN STR  TH 
#25  35  81  69 103  93  59  52  48 
# meaning : CB:cerebellum, ENT:entopeduncular, FC:frontal cortex
# GP: Globus pallidus  HC:hippocampus   PC: posterior cortex
# SN: substantia nigra STR:striatum TH:thalamus

# reserve q and w, can be of utility some day.

# # POSTERIOR CORTEX : do not open , for now
# dge_poste <- loadSparseDge("F_GRCm38.81.P60Cortex_noRep5_POSTERIORonly.raw.dge.txt.gz")
# clus_poste <- readRDS("F_GRCm38.81.P60Cortex_noRep5_POSTERIORonly.cluster.assign.RDS")

#install.packages("~/dataref_scRNA/saunders/DropSeq.util_2.0.tar.gz",repos=NULL)
library(DropSeq.util)

doviolinprep <- function(genesym, dge_frontal, genes_frontal, mydfviz ){
  geneindex <- which(genes_frontal==genesym)
  gcountsvec <- dge_frontal[geneindex,]
  # pass to mydfviz (ok, order respected)
  mydfviz$rawcnts <- gcountsvec
  mydfviz$lograw <- log10(mydfviz$rawcnts+1)
  gghere <- ggplot(data=mydfviz, aes(factor(cluster),lograw )) +
    geom_violin() + labs(title=genesym)
  print(gghere)
  return(0)
}

if( !file.exists(OUTFILE.frontal)){
  # dgTMatrix object  : dge_frontal
  # FRONTAL CORTEX : clusters info: barcodes and to which cluster they belong to:
  dge_frontal <- loadSparseDge("F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt.gz")
  clus_frontal <- readRDS("F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.cluster.assign.RDS")
  length(clus_frontal) # 156167
  ncol(dge_frontal) #  194027 
  print("caution, dge_frontal ncol not equal to clus_frontal length, use names")
  
  genes_frontal <- readLines("frontalgenes.txt")
  bc_frontal <- readLines("frontalbarcodes.txt")
  print(length(genes_frontal)); print(length(bc_frontal))
  print(dim(dge_frontal))
  print("ok my txt genes and barcodes match with matrix dimensions")
  ## ---------- get only neurons and astrocytes:  FAILED !! --------------
  desiredclass <- c("NEURON","ASTROCYTE", "OLIGODENDROCYTE")
  selq <- q %>% filter(tissue == "FC") %>% filter(class %in% desiredclass)
  clusterskept <- unname(sapply(selq$subcluster, function(x) unlist(str_split(x,"-"))[1]))
  clusterskept <- sort(as.integer(unique(as.vector(clusterskept))))
  ( clusterskept )
  mybool <- sapply(clus_frontal, function(x) ifelse(x %in% clusterskept, TRUE, FALSE))
  # mybool
  # FALSE   TRUE 
  # 15600 140567  
  # this kepts 140 000 cells, makes no big difference if I take all frontal barcodes
  rm(mybool, selq, clusterskept)
  # end FAILED
  
  # ------- good choice : clusters selection by main genes ----------------------
    # mydfviz <- data.frame("bc","cluster","log2rawcount")
  mydfviz <- data.frame("bc"=bc_frontal)
  clus_frodf <- data.frame("bc"=names(clus_frontal), "cluster"=clus_frontal)
  mydfviz <- left_join(mydfviz, clus_frodf, by="bc"); rm(clus_frodf)
    # get index where Tubb3 is
  print(doviolinprep("Tubb3",dge_frontal, genes_frontal, mydfviz))
  # same for Gfap
  print(doviolinprep("Gfap",dge_frontal, genes_frontal, mydfviz))
  # from what I saw in violinplots, pick clusters 6,8,9 (neu neu astro oligo)
  SELECLUS <- c(6,7,8,9)
  pickbcdf <- mydfviz %>% filter(cluster %in% SELECLUS)
  dim(pickbcdf)  # 
    gc() # free unused r memory , garbage collection
  colnames(dge_frontal) <- bc_frontal
  rownames(dge_frontal) <- genes_frontal
  frontal.mx <- dge_frontal[,pickbcdf$bc]
  
  qq <- q %>% filter(tissue=="FC") %>% select(class, full_name,subcluster)
  qq$cluster <- sapply(qq$subcluster, function(x){
    unlist(str_split(x, "-"))[1]
  })
  qq <- qq %>% select(-subcluster) %>% filter(cluster %in% SELECLUS) %>%
    distinct()
  cols2add <- c("class", "full_name")
  frontal.ann = cbind(pickbcdf, qq[match(pickbcdf$cluster, qq$cluster), cols2add])
  frontal.ann$cluster <- factor(frontal.ann$cluster, 
                                levels=sort(unique(frontal.ann$cluster)))
  
  
  if (ncol(frontal.mx)==nrow(frontal.ann)){
    # A. full version .rds
    saveRDS(list("matrix"=frontal.mx, "annots"=frontal.ann), 
            OUTFILE.frontal)
    
    # B.  interoperability sparse matrix R and python, reticulate
    library(reticulate)
    use_condaenv(condaenv=mycondaname,conda=myconda)
    py_config()
    scipy <- import("scipy")
    scipy$sparse$save_npz("FRONTAL_nfull_npz.npz", frontal.mx)
    #  in python open : mymx = scipy.sparse.load_npz("wow.npz")
    # end reticulate  
    # metadata to csv file: 
    write.table(frontal.ann, "FRONTAL_nfullANN.csv",col.names=TRUE,
                sep=",", row.names=FALSE)
    
    }else{ print("discordant ncol .mx and nrow .ann, saveRDS failed") }
  
  nbn = 20000; if ( table(frontal.ann$class)["NEURON"] > (nbn + 1) ){
    #  C. reduced verssion .rds if excessive neurons ! ------
    table(frontal.ann$class)
    # ASTROCYTE          NEURON OLIGODENDROCYTE 
    # 10394           88903            8718 
    # as seen, excessive nb of neurons, pick the nbn neuron expressing highest Tubb3
    #  ( better : reject (ntotalneu - nbn) neuron barcodes showing lowest Tubb3  )
    
    tmp.ann <- frontal.ann
    tmp.v.tubb3 <- frontal.mx["Tubb3",] + 0
    tmp.ann$tubb3 <- tmp.v.tubb3[match(tmp.ann$bc, names(tmp.v.tubb3))]
    ntotalneu <- tmp.ann %>% filter(class =="NEURON") %>% summarise(ntotalneu=n())
    ndrop = ntotalneu$ntotalneu - nbn
    todropneu <- tmp.ann %>% filter(class =="NEURON") %>%
      arrange(tubb3) %>% slice_min(tubb3, n = ndrop)
    
    frontal.mx <- frontal.mx[ , ! colnames(frontal.mx) %in% todropneu$bc ]
    frontal.ann <- frontal.ann[! frontal.ann$bc %in% todropneu$bc,  ]
    
    print(all(colnames(frontal.mx) == frontal.ann$bc))  # ok equality matrix et annot
    if (ncol(frontal.mx)==nrow(frontal.ann)){
      saveRDS(list("matrix"=frontal.mx, "annots"=frontal.ann), 
              OUTFILEV2.frontal)
      
      # D. same reduced in python friendly format s:
      scipy$sparse$save_npz("FRONTALV2_nfull_npz.npz", frontal.mx)
      write.table(frontal.ann, "FRONTALV2_nfullANN.csv",col.names=TRUE,
                  sep=",", row.names=FALSE)
      
    }else{ print("discordant ncol .mx and nrow .ann, saveRDS V2 failed") }
  } # end if neurons > 30 000  
 
} # end if not exists rds file

