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

