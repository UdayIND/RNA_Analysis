# Load necessary libraries
install.packages("Seurat")   # If not already installed
install.packages("ggplot2")
install.packages("dplyr")
install.packages("DoubletFInder")
# install.packages("SingleCellExperiment")
install.packages("BiocManager")
BiocManager::install("GEOquery")   # For fetching GEO datasets
BiocManager::install("DoubletFinder")   
library(Seurat)
library(ggplot2)
library(dplyr)
library(GEOquery)
library(Matrix)


# Define file paths (update these based on your directory)
mtx_file <- "/N/u/saprem/Quartz/Downloads/GSE162610_sci_mat.mtx"
barcodes_file <- "/N/u/saprem/Quartz/Downloads/GSE162610_barcodes.tsv"
genes_file <- "/N/u/saprem/Quartz/Downloads/GSE162610_genes.tsv"
cell_metadata_file <- "/N/u/saprem/Quartz/Downloads/GSE162610_barcode_metadata.tsv"
gene_metadata_file <- "/N/u/saprem/Quartz/Downloads/GSE162610_gene_metadata.tsv"

# Read expression matrix (sparse format)
expr_matrix <- readMM(mtx_file)

# Read barcodes (cell IDs)
barcodes <- read.table(barcodes_file, header = FALSE, stringsAsFactors = FALSE)

# Read gene names
genes <- read.table(genes_file, header = FALSE, stringsAsFactors = FALSE)

# Read cell metadata
cell_metadata <- read.table(cell_metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Read gene metadata
gene_metadata <- read.table(gene_metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Assign row and column names to the matrix
rownames(expr_matrix) <- genes$V2  # Assuming gene symbols are in the second column
colnames(expr_matrix) <- barcodes$V1  # Assigning cell names

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = expr_matrix, project = "SCI_scRNA", min.cells = 3, min.features = 200)

# Add cell metadata to Seurat object
seurat_obj <- AddMetaData(seurat_obj, metadata = cell_metadata)

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj)

# Identify highly variable genes
seurat_obj <- FindVariableFeatures(seurat_obj)

# Scale the data
seurat_obj <- ScaleData(seurat_obj)

# Calculate mitochondrial gene percentage (for quality filtering)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Filter out low-quality cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)

# Check QC metrics
library(ggplot2)
VlnPlot(seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0.2,  # Smaller point size for better visibility
        cols = c("#1f77b4", "#ff7f0e", "#2ca02c")  # Custom color scheme
) +
  theme_minimal() +  # Clean background
  theme(
    axis.text.x = element_text(size = 14, face = "bold", angle = 0),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  labs(
    title = "Quality Control Metrics",
    x = "Cell Identity",
    y = "Expression Level"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))  # Adjust y-axis scale

#Normalize Data (if not already normalized)
seurat_obj <- NormalizeData(seurat_obj)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
seurat_obj <- ScaleData(seurat_obj)  # Re-run scaling

# Run PCA for dimensionality reduction
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)


seurat_obj@meta.data$cell_type <- seurat_obj@meta.data$celltype  # Change to uniform name
Idents(seurat_obj) <- seurat_obj$cell_type

# Define custom color palette matching the figure
cell_type_colors <- c(
  "Neutrophil" = "#8B0000", "Monocyte" = "#A52A2A", "Macrophage" = "#DC143C",
  "Dendritic" = "#FF4500", "Microglia" = "#FFD700", "Div-Myeloid" = "#8B4513",
  "Fibroblast" = "#4682B4", "Endothelial" = "#008000", "Pericyte" = "#556B2F",
  "OPC" = "#00CED1", "Oligodendrocyte" = "#1E90FF", "Astrocyte" = "#4B0082",
  "Ependymal" = "#9400D3", "Lymphocyte" = "#FF00FF", "Neuron" = "#FFDAB9"
)
# Ensure PCA has been run before UMAP
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)

# Compute UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# Verify that UMAP is now in the Seurat object
seurat_obj@reductions  # Check if "umap" appears

# Generate UMAP plot
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = "cell_type") +
  scale_color_manual(values = cell_type_colors) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  labs(
    title = "UMAP Plot of Spinal Cord Cells",
    x = "UMAP 1",
    y = "UMAP 2"
  )
# Count number of cells per cell type
cell_counts <- table(seurat_obj$cell_type)

# Convert to a dataframe for display
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("Cell Type", "Count")

# Display table in console
print(cell_counts_df)

######## Figure 1b

# Load required libraries
library(ggplot2)
library(dplyr)
library(Seurat)

# Extract UMAP coordinates
umap_data <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
colnames(umap_data) <- c("UMAP_1", "UMAP_2")  # Standardize names

# Assign time point metadata
umap_data$time_point <- seurat_obj@meta.data$time_point

# Remove NA values in time point column
umap_data <- na.omit(umap_data)

# Order time points to match reference plot layout
umap_data$time_point <- factor(umap_data$time_point, levels = c("Uninjured", "1dpi", "3dpi", "7dpi"))

# Create UMAP density plot with correct colors and layout
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_bin2d(bins = 100) +  # Create 2D histogram bins for density
  scale_fill_gradient(low = "yellow", high = "red", trans = "log2") +  # Match paper's density color scale
  facet_wrap(~time_point, ncol = 2) +  # Ensure correct 2x2 layout
  theme_minimal() +  # Clean background
  theme(
    panel.grid = element_blank(),  # Remove unnecessary grid lines
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "UMAP of All Cells Split by Time Point (Log2 Density)",
    x = "UMAP 1",
    y = "UMAP 2",
    fill = "log2(Density)"  # Correct legend label
  )
