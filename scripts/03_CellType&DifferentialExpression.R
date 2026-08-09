# ==============================================================================
# Cell Type Annotation & Differential Expression Analysis
# Dataset: Disseminated Tuberculosis (DTB) 10X scRNA-seq
# ==============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)

# Create results directory if it doesn't exist
if (!dir.exists("results")) dir.create("results")

# ---- 1. Load Processed Seurat Object ----
dtb_merged <- readRDS("dtb_merged_processed.rds")

# ---- 2. Visualize Canonical Immune Lineage Markers ----
marker_genes <- c("CD3D", "CD4", "CD8A", "CD14", "LYZ", "FCGR3A", "MS4A1", "NKG7", "PPBP")

# 2A. FeaturePlot: Overlay marker gene expression onto UMAP
p_features <- FeaturePlot(dtb_merged, features = marker_genes, ncol = 3, pt.size = 0.1)
ggsave("results/canonical_markers_featureplot.png", plot = p_features, width = 12, height = 10, dpi = 300)

# 2B. DotPlot: Marker expression intensity & detection rate per cluster
p_dotplot <- DotPlot(dtb_merged, features = marker_genes) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Canonical Immune Markers per Cluster")
ggsave("results/canonical_markers_dotplot.png", plot = p_dotplot, width = 9, height = 6, dpi = 300)

# ---- 3. Identify Top Cluster Marker Genes ----
# Join layers (Seurat v5 compatibility) before running differential expression
dtb_merged <- JoinLayers(dtb_merged)

cluster_markers <- FindAllMarkers(dtb_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top2_markers <- cluster_markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

print(top2_markers, n = 50)
write.csv(cluster_markers, file = "results/all_cluster_markers.csv", row.names = FALSE)

# ---- 4. Assign Cell Type Identities to All 24 Clusters ----
# Complete mapping for clusters 0 through 23, this is just a quick assignment based on canonical markers and known immune cell types. Further validation may be needed.
new_cluster_ids <- c(
  "0"  = "CD8+ T cells",
  "1"  = "CD14+ Monocytes",
  "2"  = "CD4+ T cells",
  "3"  = "CD4+ T cells",
  "4"  = "B cells",
  "5"  = "CD4+ T cells",
  "6"  = "NK cells",
  "7"  = "CD8+ T cells",
  "8"  = "CD14+ Monocytes",
  "9"  = "CD14+ Monocytes",
  "10" = "CD14+ Monocytes",
  "11" = "FCGR3A+ Monocytes",
  "12" = "Unassigned/Other",
  "13" = "NK cells",
  "14" = "CD8+ T cells",
  "15" = "NK cells",
  "16" = "Unassigned/Other",
  "17" = "CD14+ Monocytes",
  "18" = "Platelets",
  "19" = "NK cells",
  "20" = "CD14+ Monocytes",
  "21" = "CD4+ T cells",
  "22" = "Unassigned/Other",
  "23" = "Unassigned/Other"
)

# Reset identity to numeric clusters and apply cell-type mapping
Idents(dtb_merged) <- "seurat_clusters"
dtb_merged <- RenameIdents(dtb_merged, new_cluster_ids)
dtb_merged$cell_type <- Idents(dtb_merged)

# ---- 5. Final Annotated UMAP Plot ----
p_annotated <- DimPlot(dtb_merged, reduction = "umap", group.by = "cell_type", label = TRUE, pt.size = 0.2, repel = TRUE) + 
               ggtitle("Disseminated TB Immune Cell Landscape") +
               theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

print(p_annotated)
ggsave("results/annotated_immune_umap.png", plot = p_annotated, width = 10, height = 7, dpi = 300)

# Save annotated Seurat object
saveRDS(dtb_merged, file = "dtb_annotated_processed.rds")

# ---- 6. Differential Expression Analysis (CD14+ Monocytes vs CD4+ T cells) ----
monocyte_vs_tcell_de <- FindMarkers(
  dtb_merged, 
  ident.1 = "CD14+ Monocytes", 
  ident.2 = "CD4+ T cells",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Display top 10 differentially expressed genes
head(monocyte_vs_tcell_de, 10)

# Save DE results to CSV
write.csv(monocyte_vs_tcell_de, file = "results/DE_CD14Monocytes_vs_CD4Tcells.csv")
cat("✅ Cell type annotation & Differential Expression analysis complete!\n")