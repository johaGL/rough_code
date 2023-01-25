
library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)

mywdir = "~/"
setwd(mywdir)
wdadir = "visium_coauthorsBMC/10XVisium_2/10XVisium_2/"
ODIR="~/spatial_hei/"
odirqc = paste0(ODIR,"QC/")
odirclus = paste0(ODIR,"clus/")
dir.create(odirqc)
dir.create(odirclus)
tisswords_ = list("C" = "Cortex", "T" = "Tumor", "TC" = "Tumor_core", "TI" = "Tumor_infiltration")

#####################
# print the images, minimal test:
####################
px = "UKF241"
tiss = "C"
px_tiss = paste0(px,"_", tiss)

pxdir = paste0("#",px_tiss,"_ST/")
print(pxdir)
image_file = paste0(wdadir, pxdir, "outs/spatial/tissue_hires_image.png")
oriimage = png::readPNG(image_file, native=TRUE)
fakedata = data.frame("a1"=c(6,7),"a2"=c(100,200))
foo = ggplot(fakedata) + geom_point(aes(x=a1,y=a2))
foo + foo
hey = foo + oriimage
hey
#ggsave(file = "test.svg", plot=hey)


# patients_all = list("UKF242" = c("C","T"), "UKF248"=c("C", "T"),
#                            "UKF256"= c("C", "TC", "TI") ,
#                          "UKF259"= c("C","T"),
#                          "UKF265"=c("C","T"), "UKF313"=c("C","T"),  "UKF334"=c("C","T")
#                         )    

####################

# title??
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

markers_ = c("LDHA", "LDHB", "SLC16A1", "SLC16A3", "GFAP", "TUBB3", "VCAM1")  # "GFAP", "TUBB3", "VCAM1", 

ploM_ = list()
mark_minmax = list()
#initialize mark_minmax and ploM
for (m in markers_){
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
    tissword = tisswords_[[tiss]]
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
    if (file.exists(paste0(ODIR,px_tiss,"_clus"))){
      oo <- readRDS(paste0(ODIR,px_tiss,"_clus"))
      
      # across seurat objects, update the min and max of the selected markers
      for (marker in markers_){
        mark_minmax = thismarker_update(oo, mark_minmax, marker)
        
        tmppatch <- SpatialFeaturePlot(oo, 
                                       features =marker,
                                       stroke = 0.1, image.alpha=0)  + 
          labs(caption=paste(px_tiss, marker)) + theme(
            plot.caption = element_text(size=9)) 
          
        
        ploM_[[marker]][[px_tiss]] <- tmppatch
      } # end for marker
      
    }else{
      badplot <- ggplot() + annotate("text", x=1, y=1,
                                     size = 3, label="quality issues") +
        theme_void()
      for (marker in markers_){
        ploM_[[marker]][[px_tiss]] <- badplot
      }
      
    } # end if
    
  } # end for tiss
}


finalouf = list()
finalouf[["histo_col"]] = histolog
for (marker in markers_){
  finalouf[[marker]] = patchwork::wrap_plots(ploM_[[marker]], ncol=1) + plot_layout(guides="collect") &
    scale_fill_distiller(palette = "RdYlBu", limits = mark_minmax[[marker]]) &
    theme(legend.position="none")
  #theme(plot.margin= margin( 0,0, 0, 0, "cm"))
}

wrapped_end = patchwork::wrap_plots(finalouf, ncol = 1+length(markers_)) # 1 for the histology

bigo = wrapped_end + plot_layout(guides="collect")  & 
  theme(legend.position="top")
ggsave(file="samples_scale_by_column.svg", plot=bigo, width=30, height=100, limitsize = FALSE )
ggsave(file="samples_scale_by_column.png", plot=bigo, width=30, height=100, limitsize = FALSE)


# new !! unique scale-------------------------------------

# use only a part of the markers: 

lessmarkers = c("LDHA", "LDHB", "SLC16A1", "SLC16A3")

maximal_offall = 0
for (k in lessmarkers){
  if (mark_minmax[[k]][2] > maximal_offall){
    maximal_offall = mark_minmax[[k]][2]
  }
}

finalouf = list()
finalouf[["histo_col"]] = histolog
for (marker in lessmarkers){
  finalouf[[marker]] = patchwork::wrap_plots(ploM_[[marker]], ncol=1) + plot_layout(guides="collect") &
    guides(colour = guide_legend(title="foo")) &
    scale_fill_distiller(palette = "RdYlBu", limits = mark_minmax[[marker]]) &
    theme(legend.position="none")
  #theme(plot.margin= margin( 0,0, 0, 0, "cm"))
}
wrapped_end = patchwork::wrap_plots(finalouf, ncol = 1+length(lessmarkers)) # 1 for the histology

bigU = wrapped_end + plot_layout(guides="collect")  & 
  scale_fill_distiller(palette = "RdYlBu",  limits = c(0, maximal_offall))  &
  theme(legend.position="top")
   # scale_fill_distiller(palette = "PuOr", trans = "log", limits = c(0, maximal_offall))  &
 
ggsave(file="samples_unique_scaleX.png", plot=bigU, width=20, height=100, limitsize = FALSE )

# 
# ### new : log scale:
# bigU = wrapped_end + plot_layout(guides="collect")  & 
#   scale_fill_distiller(palette = "RdYlBu", trans="log", limits = c(0, maximal_offall))  &
#   theme(legend.position="top")
# # scale_fill_distiller(palette = "PuOr", trans = "log", limits = c(0, maximal_offall))  &
# 
# ggsave(file="bigoNEW_log.pdf", plot=bigU, width=30, height=100, limitsize = FALSE )


# vcam = patchwork::wrap_plots(ploM_[["VCAM1"]], ncol=1) + plot_layout(guides="collect") &
#   scale_fill_distiller(palette = "PRGn", limits = mark_minmax[["VCAM1"]]) &
#   theme(legend.position="top")
# 
# gfap = patchwork::wrap_plots(ploM_[["GFAP"]], ncol=1) + plot_layout(guides="collect") &
#   scale_fill_distiller(palette = "PRGn", limits = mark_minmax[["GFAP"]]) &
#   theme(legend.position="top")
# 
# bigu = histolog | vcam | gfap | plot_layout(guides="collect")  & 
#   theme(legend.position="top")
# 
# ggsave(file="bigo.svg", plot=bigu, width=20, height=30 )



# END

## waste : 

# 
# for (px in names(patients_all)){
#   for (tiss in patients_all[[px]]){
#     px_tiss = paste0(px,"_", tiss)
#     image_file = paste0(wdadir, pxdir, "outs/spatial/tissue_hires_image.png")
#     oriimage = png::readPNG(image_file, native=TRUE)
#     tryCatch({
#       oo <- readRDS(paste0(ODIR,px_tiss,"_clus"))
#       
#       a = min(oo@assays$SCT[markers_[1],])
#       b = max(oo@assays$SCT[markers_[1],])
#       mark_minmax[[markers_[1]]] = update_minmax(a, b, mark_minmax[[markers_[1]]])
#       
#       ploM_[[px_tiss]] <- SpatialFeaturePlot(oo, features = markers_[1]) # +
#               #     scale_fill_continuous(limits=c(-2,2) )
#       print(oo)
#     }, error=function(cond){
#       badplot <- ggplot() + annotate("text", x=1, y=1,
#                                      size = 7, label="bad quality, not plottable") +
#         theme_void()
#       ploM_[[px_tiss]] <- badplot
#       print("no clus file")
#       message(cond)
#       next
#     })
#   }
# }
