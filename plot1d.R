# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)

# Define the genes that best identify each cell type
marker_genes <- c(
  "Ltf",  # Neutrophil
  "Ccr2",  # Monocyte
  "Ms4a7",  # Macrophage
  "Gpr84",  # Microglia
  "Mki67",  # Div-Myeloid
  "Col1a1",  # Fibroblast
  "Pecam1",  # Endothelial
  "Higd1b",  # Pericyte
  "Tnr",  # OPC
  "Ermn",  # Oligodendrocyte
  "Slc6a11",  # Astrocyte
  "Pifo"  # Ependymal
)
# Normalize data if not done already
if (!"SCT" %in% names(seurat_obj@assays)) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
}

# Run PCA
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)

# Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# Verify UMAP is available
print(seurat_obj@reductions)


# Generate UMAP feature plots
umap_plots <- lapply(marker_genes, function(gene) {
  FeaturePlot(
    seurat_obj, 
    features = gene, 
    cols = c("lightgrey", "red"), 
    min.cutoff = "q10", 
    max.cutoff = "q90"
  ) + ggtitle(gene)
})

# Arrange plots in a grid
final_plot <- wrap_plots(umap_plots, ncol = 4)

# Display the final UMAP plot grid
print(final_plot)
