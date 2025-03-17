# Define myeloid cell types based on metadata
myeloid_cells <- c("Neutrophil", "Monocyte", "Macrophage", "Microglia", 
                   "Div-Myeloid", "Dendritic")

# Subset Seurat object to retain only myeloid cells
myeloid_seurat_obj <- subset(seurat_obj, subset = celltype %in% myeloid_cells)

# Normalize and scale the data
myeloid_seurat_obj <- NormalizeData(myeloid_seurat_obj)
myeloid_seurat_obj <- FindVariableFeatures(myeloid_seurat_obj)
myeloid_seurat_obj <- ScaleData(myeloid_seurat_obj)

# Run PCA
myeloid_seurat_obj <- RunPCA(myeloid_seurat_obj, npcs = 30)

# Run UMAP
myeloid_seurat_obj <- RunUMAP(myeloid_seurat_obj, dims = 1:30)

# Run clustering
myeloid_seurat_obj <- FindNeighbors(myeloid_seurat_obj, dims = 1:30)
myeloid_seurat_obj <- FindClusters(myeloid_seurat_obj, resolution = 0.5)


# Assign colors to myeloid subtypes
myeloid_colors <- c(
  "Neutrophil" = "#E41A1C",
  "Monocyte" = "#377EB8",
  "Macrophage" = "#4DAF4A",
  "Microglia" = "#984EA3",
  "Div-Myeloid" = "#FF7F00",
  "Dendritic" = "#FFFF33"
)

# Plot UMAP of myeloid cells
DimPlot(myeloid_seurat_obj, 
        reduction = "umap", 
        group.by = "celltype", 
        cols = myeloid_colors,
        label = TRUE, 
        repel = TRUE) +
  ggtitle("UMAP of Myeloid Cells After SCI") +
  theme_minimal()
