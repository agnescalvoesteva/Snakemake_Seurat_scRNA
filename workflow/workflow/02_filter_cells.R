#!/usr/bin/env Rscript

# Load required libraries
library(Seurat)
library(ggplot2)

# Load configuration and paths from Snakemake
config <- snakemake@config
input_file <- snakemake@input$rds
output_obj <- snakemake@output$filtered

# Create output directory if it does not exist
dir.create(dirname(output_obj), recursive = TRUE, showWarnings = FALSE)

# Load previously generated Seurat object
pbmc <- readRDS(input_file)

# Filter cells based on defined criteria (number of features and mitochondrial percentage)
pbmc_filtered <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Generate statistics before and after filtering
filter_stats <- data.frame(
  stage = c("Pre-filter", "Post-filter"),
  cells = c(ncol(pbmc), ncol(pbmc_filtered)),
  features = c(mean(pbmc$nFeature_RNA), mean(pbmc_filtered$nFeature_RNA))
)

# Print statistics to console
print(filter_stats)

# Save the filtered Seurat object
saveRDS(pbmc_filtered, file = output_obj)

# Final confirmation message
cat("Filtering completed. Cells:", ncol(pbmc_filtered), "\n")
