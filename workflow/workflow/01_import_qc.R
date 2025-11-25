#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(patchwork)

# Read config.yaml
config <- snakemake@config
sample <- snakemake@wildcards$sample
input_file <- snakemake@input$h5
output_obj <- snakemake@output$rds
output_plot <- snakemake@output$qc_plot

# Make a carpet if id not exist one before that. 
dir.create(dirname(output_plot), recursive = TRUE, showWarnings = FALSE)

cat("Procesando muestra:", sample, "\n")
cat("Archivo de entrada:", input_file, "\n")

# Try catch the errors
tryCatch({
  # Download the data of input_file (h5_file)
  pbmc.data <- Read10X_h5(input_file)

  # Create a seurat object 
  pbmc <- CreateSeuratObject(
    counts = pbmc.data,
    project = config$project$name,
    min.cells = config$qc$min_cells,
    min.features = config$qc$min_features
  )

  # calculate the percentatge 
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = config$qc$mt_pattern)

  # Features of the graphic
  if (is.null(config$plots$violin_features) || length(config$plots$violin_features) == 0) {
    stop("config$plots$violin_features no está definido o está vacío")
  }

  # Visualization of QC
  qc_plot <- VlnPlot(pbmc, features = config$plots$violin_features, ncol = 3) +
    ggtitle(paste("QC Metrics -", sample))

  # Save the graphic
  ggsave(filename = output_plot, plot = qc_plot, width = 12, height = 8)

  # Save the RDS
  saveRDS(pbmc, file = output_obj)

  cat("QC completado para:", sample, "\n")

}, error = function(e) {
  cat("ERROR en la muestra", sample, ":\n")
  cat(e$message, "\n")
  stop(e)  
})
