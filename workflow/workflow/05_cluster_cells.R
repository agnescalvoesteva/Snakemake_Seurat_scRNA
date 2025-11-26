#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)

config <- snakemake@config
input_file <- snakemake@input$pca
output_obj <- snakemake@output$clustered
output_umap <- snakemake@output$umap   # Now Snakemake controls the path

dir.create(dirname(output_obj), recursive = TRUE, showWarnings = FALSE)

# Load object
pbmc <- readRDS(input_file)

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:config$neighbors$dims)
pbmc <- FindClusters(pbmc, resolution = config$clustering$resolution)

# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:config$umap$dims)

# Save UMAP
umap_plot <- DimPlot(pbmc, reduction = "umap", 
                     label = config$umap$label_clusters,
                     pt.size = config$umap$point_size)
ggsave(output_umap, umap_plot, width = 8, height = 6)

# Save object with clusters and UMAP
saveRDS(pbmc, file = output_obj)

cat("Clustering completed. Number of clusters:", length(unique(Idents(pbmc))), "\n")
