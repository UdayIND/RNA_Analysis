# Load necessary libraries
library(Seurat)
library(ggplot2)
library(patchwork)

# Define **30 marker genes** from the paper
marker_genes <- c(
  # Microglia markers
  "P2ry12", "Tmem119", "Siglech",  
  # Peripheral Myeloid Cells
  "Ccr2", "Ly6c2",                 
  # Homeostatic Microglia (H-Microglia)
  "P2ry12", "Siglech", "Tmem119",  
  # Inflammatory Microglia
  "Igf1",  
  # Migrating Microglia
  "Msr1", "Igf1",                   
  # Dividing Myeloid Cells
  "Mki67", "Top2a",                 
  # Monocytes
  "Ccr2", "Ly6c2",                 
  # Macrophages
  "Cd63",  
  # Additional myeloid and inflammation-related genes from the study
  "Apoe", "Cd74", "Csf1r", "Cst3", "H2-Ab1", "H2-Aa", "H2-DMb1", "Ifitm3", "Trem2",
  "Fcgr1", "Itgax", "C1qa", "C1qb", "C1qc", "Spp1", "Cd163", "Ly6a", "Ly86", "Cd14"
)

# Ensure all genes are present in the dataset
available_genes <- marker_genes[marker_genes %in% rownames(myeloid_seurat_obj)]

# Generate UMAP Feature Plots for all available genes
umap_plots <- lapply(available_genes, function(gene) {
  FeaturePlot(
    myeloid_seurat_obj, 
    features = gene, 
    cols = c("lightgrey", "red"),  # Color: Grey (low expression) → Red (high expression)
    min.cutoff = "q10", max.cutoff = "q90",
    reduction = "umap"
  ) + ggtitle(gene)
})

# Arrange plots in a grid layout (adjust columns for better display)
final_plot <- wrap_plots(umap_plots, ncol = 5)

# Display the UMAP plots
print(final_plot)
