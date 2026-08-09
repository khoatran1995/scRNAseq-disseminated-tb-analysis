# ==============================================================================
# Phase 2: Normalization, PCA, Clustering & UMAP Dimensionality Reduction
# Dataset: Disseminated Tuberculosis (DTB) 10X scRNA-seq
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)

# Create results directory if it doesn't exist
if (!dir.exists("results")) dir.create("results")

# ---- 1. Load Clean QC-Filtered Data ----
seurat_list <- readRDS("seurat_list_qc_filtered.rds")

# ---- 2. Merge All Samples into One Seurat Object ----
dtb_merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))

# ---- 3. Normalize Data & Identify Top Variable Genes ----
# Log-normalize counts to account for sequencing depth differences
dtb_merged <- NormalizeData(dtb_merged, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify top 2,000 highly variable genes (HVGs) driving cell heterogeneity
dtb_merged <- FindVariableFeatures(dtb_merged, selection.method = "vst", nfeatures = 2000)

# ---- 4. Scale Data & Run PCA ----
# Shift gene expression so mean = 0 and variance = 1
dtb_merged <- ScaleData(dtb_merged)

# Perform Principal Component Analysis (PCA) on top 2,000 HVGs
dtb_merged <- RunPCA(dtb_merged, npcs = 30, verbose = FALSE)

# ---- 5. Cell Clustering & UMAP Dimensionality Reduction ----
# Construct K-nearest neighbor (KNN) graph using top 20 PCs
dtb_merged <- FindNeighbors(dtb_merged, dims = 1:20)

# Identify cell clusters (resolution = 0.5)
dtb_merged <- FindClusters(dtb_merged, resolution = 0.5)

# Generate 2D UMAP embedding
dtb_merged <- RunUMAP(dtb_merged, dims = 1:20)

# ---- 6. Visualize UMAP Plots ----
# Plot A: Colored by Cell Clusters
p1 <- DimPlot(dtb_merged, reduction = "umap", label = TRUE, pt.size = 0.1) + 
      ggtitle("Single-Cell Clusters (Res 0.5)")

# Plot B: Colored by Sample ID (Checks sample distribution across clusters)
p2 <- DimPlot(dtb_merged, reduction = "umap", group.by = "sample", pt.size = 0.1) + 
      ggtitle("Cells by Sample ID")

# Display side-by-side
p_combined <- p1 + p2
print(p_combined)

# Export UMAP figure for GitHub documentation
ggsave("results/umap_clusters_and_samples.png", plot = p_combined, width = 12, height = 5, dpi = 300)

# ---- 7. Save Processed Object ----
saveRDS(dtb_merged, file = "dtb_merged_processed.rds")
cat("✅ Phase 2 complete! Processed object saved to dtb_merged_processed.rds\n")