#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)

config <- snakemake@config
input_file <- snakemake@input$norm_hvg     # matches Snakefile
output_file <- snakemake@output$pca        # matches Snakefile

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

pbmc <- readRDS(input_file)

# Scale data
if (config$scale_data$features == "all") {
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
}

# PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = config$pca$dims)

# Save object with PCA
saveRDS(pbmc, file = output_file)

cat("Scaling and PCA completed\n")
