
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
ggsave(file = "test.svg", plot=hey)


# patients_all = list("UKF242" = c("C","T"), "UKF248"=c("C", "T"),
#                            "UKF256"= c("C", "TC", "TI") ,
#                          "UKF259"= c("C","T"),
#                          "UKF265"=c("C","T"), "UKF313"=c("C","T"),  "UKF334"=c("C","T")
#                         )    

####################

# title??
####################
patients_all = list("UKF242" = c("C","T"), "UKF265"=c("T", "C"))

fooli = list()
markers_ = c("VCAM1")
mark_minmax = list("VCAM1" = c(0,0))
#initialize mark_minmax
for (m in markers_){
  
}


update_minmax <- function(a, b, min_max_pair){
  if (a < min_max_pair[1]){
    min_max_pair[1] = a
  }
  if (b > min_max_pair[2]){
    min_max_pair[2] = b
  }
  return(min_max_pair)
}


# histology : 
fiili = list()
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
    fiili[[px_tiss]] = mins
  }
}


histolog = patchwork::wrap_plots(fiili, ncol=1)
pdf("foooho.pdf")
print(histolog)
dev.off()
# markers

for (px in names(patients_all)){
  for (tiss in patients_all[[px]]){
    px_tiss = paste0(px,"_", tiss)
    image_file = paste0(wdadir, pxdir, "outs/spatial/tissue_hires_image.png")
    oriimage = png::readPNG(image_file, native=TRUE)
    if (file.exists(paste0(ODIR,px_tiss,"_clus"))){
      oo <- readRDS(paste0(ODIR,px_tiss,"_clus"))
      
      # across seurat objects, update the min and max of the selected markers
      a = min(oo@assays$SCT[markers_[1],])
      b = max(oo@assays$SCT[markers_[1],])
      mark_minmax[[markers_[1]]] = update_minmax(a, b, mark_minmax[[markers_[1]]])
      
      fooli[[px_tiss]] <- SpatialFeaturePlot(oo, features = markers_[1]) # + 
        # labs(caption=px_tiss) + theme(
        #   plot.caption = element_text(size=5)
        # )
      #     scale_fill_continuous(limits=c(-2,2) )
      print(oo)

    }else{
      badplot <- ggplot() + annotate("text", x=1, y=1,
                                     size = 3, label="quality issues,\n not plottable") +
        theme_void()
      fooli[[px_tiss]] <- badplot
    }

  }
}



#fooli$UKF242_C + fooli$UKF242_T + fooli$UKF265_C + plot_layout(guides="collect") &
#  scale_fill_continuous(type = "viridis", limits = mark_minmax[[markers_[1]]] ) 

#fooli$UKF242_C + fooli$UKF242_T + fooli$UKF265_C + plot_layout(guides="collect") &
#  scale_fill_distiller(palette = "Spectral", limits = mark_minmax[[markers_[1]]])



onemark = patchwork::wrap_plots(fooli, ncol=1) + plot_layout(guides="collect") &
  scale_fill_distiller(palette = "PRGn", limits = mark_minmax[[markers_[1]]])


bigo = histolog | onemark
ggsave(file="bigo.svg", plot=bigo )

pdf("more_crazyness.pdf")
bigo
dev.off()

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
#       fooli[[px_tiss]] <- SpatialFeaturePlot(oo, features = markers_[1]) # +
#               #     scale_fill_continuous(limits=c(-2,2) )
#       print(oo)
#     }, error=function(cond){
#       badplot <- ggplot() + annotate("text", x=1, y=1,
#                                      size = 7, label="bad quality, not plottable") +
#         theme_void()
#       fooli[[px_tiss]] <- badplot
#       print("no clus file")
#       message(cond)
#       next
#     })
#   }
# }
