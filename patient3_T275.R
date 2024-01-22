library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(viridis)
ODIR = "~/spatial_thesis/data_and_analyses/spatial_transcriptomics/spatialtranscripto_germanteam/spatial_hei_results/"
px_tiss = "UKF275_T"
alternative_name = "patient3"
obj <- readRDS(paste0(ODIR,px_tiss,"_clus")) # TODO go back to spatial_transc_v2.R, add extension .rds

# plot the clusters in the spatial image
# nb_clusters = 9 # TODO: automatize this !  
# named_vector_clusters <- c("blue", "lightblue", "green", 
#       "gold", "salmon","gray","black", "pink","cyan", "red")
# 
# names(named_vector_clusters) <- seq(0:nb_clusters)


p0 <- SpatialPlot(obj,  
                     label.size = 3, 
                     alpha=0, 
                     label= FALSE,
                     label.box=FALSE) +
  scale_fill_brewer(palette="Paired")  # alpha zero yields H&E alone

# the scatterplot
p1 <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size=0.5) +
  #scale_color_manual(values=named_vector_clusters)
  scale_colour_brewer(palette="Paired") +
  theme_test() + labs(title=paste0("Seurat clusters, ", px_tiss))

p2 <- SpatialDimPlot(obj, label = TRUE, 
                     label.size = 3, 
                     alpha=0.7) +
  scale_fill_brewer(palette="Paired")

p3 <- SpatialDimPlot(obj, label = TRUE, 
                     label.size = 3, 
                     alpha=0) +
  scale_fill_brewer(palette="Paired")  # alpha zero yields H&E alone
  
svg(paste0(ODIR,"other_images/1-histol_HE_",px_tiss,".svg"), height=10, width=12)
p0
dev.off()

svg(paste0(ODIR,"other_images/2-unsup_clustering_spots_",px_tiss,".svg"), height=7, width=8)
p1
dev.off()

svg(paste0(ODIR,"other_images/3-viz-v1_spots_per_cluster_",px_tiss,".svg"), height=10, width=12)
p2
dev.off()

svg(paste0(ODIR,"other_images/4-histol_plus_clusters-",px_tiss,".svg"), height=10, width=12)
p3
dev.off()

# markers_cluster2_df <- FindMarkers(obj, ident.1 = 2, min.pct = 0.25)

all_markers <- FindAllMarkers(obj, 
                              only.pos=TRUE, # set to False if thomas wants to see down markers as well
                              min.pct = 0.25,
                              logfc.threshold=0.1)

top10 <- all_markers %>%
         group_by(cluster) %>%
         top_n(n = 10, wt = avg_log2FC)

pdf(paste0(ODIR,"other_images/5-heatmap_",px_tiss,".pdf"), height=17, width=12)
DoHeatmap(obj, features = top10$gene) +
  scale_colour_brewer(palette="Paired") +
  scale_fill_viridis()
dev.off()

output_df <- all_markers
rownames(output_df) <- seq(1:dim(output_df)[1])
output_df <- tibble::rownames_to_column(output_df, "id")
write.table(as.data.frame(output_df), 
            file=paste0(ODIR,
                        "other_images/all-markers_per_cluster_patient3_", px_tiss,".tsv"), 
            sep="\t", row.names=FALSE,
            col.names=TRUE )

# filter that table
df <- read.table(file=paste0(ODIR,
                  "other_images/all-markers_per_cluster_patient3_", px_tiss, ".tsv"), 
sep="\t", header=TRUE )
top30 <- df %>% group_by(cluster) %>% top_n(n=30, wt= avg_log2FC)
for (i in 0:max(top30$cluster)){
  tmp_df <- top30 %>% dplyr::filter(cluster == i)
  tmp_df$id <- tmp_df$gene
  colnames(tmp_df)[1] <- "gene_symbol"
  write.table(as.data.frame(tmp_df), 
              file=paste0(ODIR,
                          "other_images/markers_cluster-",i,
                          "_patient3_", px_tiss,
                          ".tsv"), 
              sep="\t", row.names=FALSE,
              col.names=TRUE )
}
