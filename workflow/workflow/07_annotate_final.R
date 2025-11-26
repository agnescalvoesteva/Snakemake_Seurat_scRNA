#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

config <- snakemake@config

# Inputs
input_file <- snakemake@input$clustered
markers_file <- snakemake@input$markers

# Outputs
output_rds <- snakemake@output$final_rds
output_umap <- snakemake@output$umap_img
output_heatmap <- snakemake@output$heatmap

output_dir <- dirname(output_rds)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load object
pbmc <- readRDS(input_file)

# 1. Annotate clusters
new.cluster.ids <- config$rename_clusters
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# 2. QC Violin plots
vln_file <- file.path(output_dir, "VlnPlot_QC.png")
png(vln_file, width = 2400, height = 1800, res = 200)
print(VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()

# 3. FeatureScatter plots
fs_file <- file.path(output_dir, "FeatureScatter.png")
png(fs_file, width = 2400, height = 1800, res = 200)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)
dev.off()

# 4. Variable feature plots
vf_file <- file.path(output_dir, "VariableFeatures.png")
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png(vf_file, width = 2400, height = 1800, res = 200)
print(plot1 + plot2)
dev.off()

# 5. PCA plots
pca_file <- file.path(output_dir, "PCA.png")
png(pca_file, width = 2400, height = 1800, res = 200)
print(DimPlot(pbmc, reduction = "pca") + NoLegend())
dev.off()

# 6. DimHeatmaps
heatmap_file <- file.path(output_dir, "DimHeatmap.png")
png(heatmap_file, width = 2400, height = 1800, res = 200)
print(DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE))
dev.off()

# 7. ElbowPlot
elbow_file <- file.path(output_dir, "ElbowPlot.png")
png(elbow_file, width = 2400, height = 1800, res = 200)
print(ElbowPlot(pbmc))
dev.off()

# 8. Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)
umap_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = config$umap$point_size) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 12))
ggsave(output_umap, plot = umap_plot, width = 12, height = 8, dpi = 300)

# 9. FeaturePlot / Heatmap
features <- config$plots$feature_plot
feature_plot <- FeaturePlot(pbmc, features = features)
png(output_heatmap, width = 2400, height = 2000, res = 200)
print(feature_plot)
dev.off()

# Save final object
saveRDS(pbmc, file = output_rds)

cat("Annotation completed and all plots saved\n")
