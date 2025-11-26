#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(patchwork)

config <- snakemake@config

# 🔹 Use correct keys
input_file <- snakemake@input$filtered
output_obj  <- snakemake@output$norm_hvg

# Output directory
dir.create(dirname(output_obj), recursive = TRUE, showWarnings = FALSE)

# Load filtered object
pbmc <- readRDS(input_file)

# Normalization
pbmc <- NormalizeData(pbmc,
                     normalization.method = config$normalization$method,
                     scale.factor = config$normalization$scale_factor)

# Find variable features
pbmc <- FindVariableFeatures(pbmc,
                            selection.method = config$variable_features$method,
                            nfeatures = config$variable_features$nfeatures)

# Save object
saveRDS(pbmc, file = output_obj)

cat("Normalization and HVG completed\n")
