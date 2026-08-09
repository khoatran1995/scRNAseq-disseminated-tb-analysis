# ==============================================================================
# Phase 1: Data Loading & Quality Control (QC) Filtering
# Dataset: Disseminated Tuberculosis (GSE287288) - 7 Samples (DTB_01 to DTB_07)
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)

# Create results directory if it doesn't exist
if (!dir.exists("results")) dir.create("results")

# ---- 1. Load 7 DTB Samples ----
base_dir <- "/GSE287288_dataset"
samples  <- sprintf("DTB_%02d", 1:7)

seurat_list <- lapply(samples, function(s) {
  counts <- Read10X(data.dir = file.path(base_dir, s))    #load the matrix of counts genes x cells
  obj <- CreateSeuratObject(counts = counts, project = s, min.cells = 3, min.features = 200)
  obj$sample <- s
  obj
})
names(seurat_list) <- samples

# ---- 2. Compute Mitochondrial Percentage (percent.mt) ----
seurat_list <- lapply(seurat_list, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
})

# ---- 3. Merge Objects to Inspect QC Metrics Before Filtering ----
merged_qc <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = samples)

p_qc <- VlnPlot(merged_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = "sample", ncol = 3, pt.size = 0)
print(p_qc)

# Save QC violin plot for GitHub documentation
ggsave("results/qc_violin_plots.png", plot = p_qc, width = 12, height = 4.5, dpi = 300)

# ---- 4. Set Filtering Thresholds ----
qc_min_features <- 200
qc_max_features <- 6000
qc_max_mt       <- 15

# ---- 5. Apply Cell Quality Filtering ----
seurat_list_filtered <- lapply(seurat_list, function(obj) {
  subset(obj, subset = nFeature_RNA > qc_min_features &
                        nFeature_RNA < qc_max_features &
                        percent.mt < qc_max_mt)
})

# ---- 6. Check Cell Counts Before vs. After QC ----
cat("\n=== CELL COUNT RETENTION SUMMARY ===\n")
for (s in samples) {
  before <- ncol(seurat_list[[s]])
  after  <- ncol(seurat_list_filtered[[s]])
  retention <- round((after / before) * 100, 2)
  cat(sprintf("%s : %5d -> %5d cells (%6.2f%% retained)\n", s, before, after, retention))
}

# ---- 7. Save Cleaned Data ----
saveRDS(seurat_list_filtered, file = "seurat_list_qc_filtered.rds")
cat("✅ Phase 1 complete! Cleaned list saved to seurat_list_qc_filtered.rds\n")