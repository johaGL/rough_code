library(Cardinal)
# mywdir <- "spatial_thesis/data_and_analyses/spatial_metabolomics/"
# dir_imzml <- "datametabo_katia/test3-gbm-lac-dan-1-1/test3-gbm-lac-dan-1-1"
# data <- readImzML(paste0(mywdir,dir_imzml))
# 
# image_data$Breaks=1
# count=1
# for(i in 2:nrow(image_data)){
#   
#   if(c(image_data$y[i]-image_data$y[i-1])>=2 | c(image_data$y[i]-image_data$y[i-1])<=c(-2) ){
#     count=count+1
#     image_data$Breaks[i]=count
#     message(count)
#   }else{image_data$Breaks[i]=count}
#   
# }
# 
# mz <- data@featureData@mz
# spectra <- as.matrix(spectraData(data)[[1]])
# rownames(spectra) <- mz
# annotation <- data.frame(mz=as.numeric(rownames(spectra)), anno=0)
# dim(annotation)
mywdir <- "~/spatial_thesis/data_and_analyses/spatial_metabolomics/"
file_mass <- "datametabo_katia/test3-gbm-lac-dan-1-1/test3-gbm-lac-dan-1-1.imzML"
m <- readMSIData(paste0(mywdir,file_mass))
m
# “MassDataFrame” object: elements of 'mz' must be unique : repeated mz ? check:
fd = featureData(m)
head(fd@mz)

n_occur <- data.frame(table(fd@mz))
reps_elems <- n_occur[n_occur$Freq > 1]

# #pca <- Cardinal::PCA(m, ncom=3) # ! danger, deprecated and kills laptop
m_mean <- summarizeFeatures(m, "mean")

# m_ref <- m_mean %>% 
#   peakPick(SNR=3) %>% 
#   peakAlign(ref="mean",
#             tolerance=0.5,
#             units="mz") %>%
#   peakFilter() %>%
#   process()

m_ref <- peakPick(m_mean, SNR=3, method= "simple")
m_ref <- peakAlign(m_ref,  tolerance=0.5, units="mz" )
m_ref <- m_ref %>% peakFilter()
m_ref <- m_ref %>% process()

m_peaks <- m %>% 
  normalize(method='tic') %>%
  peakBin(ref=mz(m_ref),
          tolerance=0.5,
          units="mz") %>%
  process()

m_peaks

set.seed(1)

m_clus <- spatialShrunkenCentroids(m_peaks,
                                   method="adaptive",
                                   r=2,
                                   s=c(20,30,40, 50), k=10)

summary(m_clus)
saveRDS(m_clus, paste0(mywdir, "results/katia_explore_clus_object.rds"))

image(m_clus, xlim=c(400, 1300), ylim=c(100, 300))

# depict the top features

topfeats2 <- topFeatures(m_clus, model=list(s=20), class==2)

image(m, mz= topfeats2$mz[5],
      xlim=c(400, 1300), ylim=c(100, 300) )

topfeats3 <- topFeatures(m_clus,model=list(s=20), 
                         class==3)

image(m, mz= topfeats3$mz[1],
      xlim=c(400, 1300), ylim=c(100, 300) )


image(m, mz= topfeats3$mz[5],
      xlim=c(400, 1300), ylim=c(100, 300) )


topfeats4 <- topFeatures(m_clus,model=list(s=20), 
                         class==4)

image(m, mz= topfeats4$mz[5],
      xlim=c(400, 1300), ylim=c(100, 300) )


topfeats5 <- topFeatures(m_clus,model=list(s=20), 
                         class==5)

image(m, mz= topfeats5$mz[5],
      xlim=c(400, 1300), ylim=c(100, 300) )
