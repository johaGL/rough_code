
library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(viridis)

mywdir = "~/"
setwd(mywdir)
wdadir = "spatialtranscripto_germanteam/10XVisium_2/10XVisium_2/"
ODIR="~/spatial_hei_results/"
odirqc = paste0(ODIR,"QC/")
odirclus = paste0(ODIR,"clus/")
dir.create(odirqc)
dir.create(odirclus)
tisswords_l = list("C" = "Cortex", "T" = "Tumor", "TC" = "Tumor_core", "TI" = "Tumor_infiltration")

#####################
# selecting patients and markers
####################
patients_all = list("UKF242" = c("C","T"), "UKF265"=c( "C", "T"))

patients_all = list("UKF242" = c("C","T"), "UKF248"=c("C", "T"),
                    "UKF256"= c("C", "TC", "TI") ,
                    "UKF259"= c("C","T"),
                    "UKF265"=c("C","T"), "UKF313"=c("C","T"),  "UKF334"=c("C","T"),
                    # group 2 :
                    "UKF241" = c("C"),
                    "UKF243" = c("T"), "UKF251" = c("T"),
                    "UKF255" = c("T"), "UKF260" = c("T"), "UKF262" = c("T"),
                    "UKF266" = c("T"), "UKF269" = c("T"), "UKF275" = c("T"),
                    "UKF296" = c("T"), "UKF304" = c("T"))

markers_l = c("CHI3L1", "CA9", "PC")  # "GFAP", "TUBB3", "VCAM1", 

ploM_ = list()
mark_minmax = list()
#initialize mark_minmax and ploM
for (m in markers_l){
  mark_minmax[[m]]  = c(0,0)
  ploM_[[m]] = NULL
}



# histology : 
histolog_pics_ = list()
for (px in names(patients_all)){
  for (tiss in patients_all[[px]]){
    px_tiss = paste0(px,"_", tiss)
    pxdir = paste0("#",px_tiss,"_ST/")
    image_file = paste0(wdadir, pxdir, "outs/spatial/tissue_hires_image.png")
    print(image_file)
    oriimage = png::readPNG(image_file, native=TRUE)
    tissword = tisswords_l[[tiss]]
    mons = ggplot() + annotate("text", x=0, y=1, label=paste(px_tiss, '\n', tissword), angle = 90) +
      theme_void()
    mins = mons + oriimage + plot_layout(widths=c(1,10))
    histolog_pics_[[px_tiss]] = mins
  }
}


histolog = patchwork::wrap_plots(histolog_pics_, ncol=1)
# pdf("foooho.pdf")
# print(histolog)
# dev.off()

################## markers plots

update_minmax <- function(a, b, min_max_pair){
  if (a < min_max_pair[1]){
    min_max_pair[1] = a
  }
  if (b > min_max_pair[2]){
    min_max_pair[2] = b
  }
  return(min_max_pair)
}

thismarker_update <- function(oo, mark_minmax, marker){
  tryCatch({
    a = min(oo@assays$SCT[marker,])
    b = max(oo@assays$SCT[marker,])
    mark_minmax[[marker]] = update_minmax(a, b, mark_minmax[[marker]])
    return(mark_minmax)
  }, error=function(cond){
    message(cond)
    return(mark_minmax)} )
}

for (px in names(patients_all)){
  for (tiss in patients_all[[px]]){
    px_tiss = paste0(px,"_", tiss)
    image_file = paste0(wdadir, pxdir, "outs/spatial/tissue_hires_image.png")
    oriimage = png::readPNG(image_file, native=TRUE)
    print("as clustered rds exists, open it")
    if (file.exists(paste0(ODIR,px_tiss,"_clus"))){
      oo <- readRDS(paste0(ODIR,px_tiss,"_clus"))
      
      # across seurat objects, update the min and max of the selected markers
      for (marker in markers_l){
        mark_minmax = thismarker_update(oo, mark_minmax, marker)
        
        tmppatch <- SpatialFeaturePlot(oo, 
                                       features =marker,
                                       pt.size.factor = 1.3, # spot size
                                       alpha = c(0.2, 1), # spot alpha, min and max transparency
                                       stroke = 0, image.alpha=0.6)  + 
          labs(caption=paste(px_tiss, marker)) + theme(
            plot.caption = element_text(size=9)) 
          
        
        ploM_[[marker]][[px_tiss]] <- tmppatch
      } # end for marker
      
    }else{
      badplot <- ggplot() + annotate("text", x=1, y=1,
                                     size = 3, label="quality issues") +
        theme_void()
      for (marker in markers_l){
        ploM_[[marker]][[px_tiss]] <- badplot
      }
      
    } # end if
    
  } # end for tiss
}


finalouf = list()
finalouf[["histo_col"]] = histolog
for (marker in markers_l){
  finalouf[[marker]] = patchwork::wrap_plots(ploM_[[marker]], ncol=1) + plot_layout(guides="collect") &
    scale_fill_distiller(palette = "Spectral",#"RdYlBu",
                         limits = mark_minmax[[marker]]) &
    theme(legend.position="none")
  #theme(plot.margin= margin( 0,0, 0, 0, "cm"))
}

wrapped_end = patchwork::wrap_plots(finalouf, ncol = 1+length(markers_l)) # 1 for the histology

bigo = wrapped_end + plot_layout(guides="collect")  & 
  theme(legend.position="top")
ggsave(file="Nantes_samples_scale_by_column.svg", plot=bigo, width=30, height=100, limitsize = FALSE )
ggsave(file="Nantes_samples_scale_by_column.png", plot=bigo, width=30, height=100, limitsize = FALSE)


# new !! unique scale-------------------------------------

# use the same markers

maximal_offall = 0
for (k in markers_l){
  if (mark_minmax[[k]][2] > maximal_offall){
    maximal_offall = mark_minmax[[k]][2]
  }
}

finalouf = list()
finalouf[["histo_col"]] = histolog
for (marker in markers_l){
  finalouf[[marker]] = patchwork::wrap_plots(ploM_[[marker]], ncol=1) + plot_layout(guides="collect") &
    guides(colour = guide_legend(title="foo")) &
    scale_fill_distiller(palette = "Spectral", limits = mark_minmax[[marker]]) &
    theme(legend.position="none")
  #theme(plot.margin= margin( 0,0, 0, 0, "cm"))
}
wrapped_end = patchwork::wrap_plots(finalouf, ncol = 1+length(markers_l)) # 1 for the histology

bigU = wrapped_end + plot_layout(guides="collect")  & 
  scale_fill_distiller(palette = "Spectral",  limits = c(0, maximal_offall))  &
  theme(legend.position="top")
   # scale_fill_distiller(palette = "PuOr", trans = "log", limits = c(0, maximal_offall))  &
 
ggsave(file="Nantes_samples_unique_scaleX.svg", plot=bigU, width=20, height=100, limitsize = FALSE )
ggsave(file="Nantes_samples_unique_scaleX.png", plot=bigU, width=20, height=100, limitsize = FALSE )

# END


