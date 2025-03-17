# Load necessary libraries
library(Seurat)
library(ggplot2)
library(MASS)
library(patchwork)

# Step 1: Ensure UMAP exists, otherwise compute it
if (!"umap" %in% names(seurat_obj@reductions)) {
  print("UMAP reduction missing, running UMAP...")
  seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)  # Ensure PCA is computed
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
}

# Step 2: Subset Myeloid Cells (Neutrophils, Monocytes, Macrophages, etc.)
myeloid_cells <- c("Neutrophil", "Monocyte", "Macrophage", "Microglia", "Div-Myeloid", "Dendritic")

# Subset Seurat object for myeloid cells only
myeloid_seurat_obj <- subset(seurat_obj, subset = celltype %in% myeloid_cells)

# Ensure time is a factor with correct levels
myeloid_seurat_obj$time <- factor(myeloid_seurat_obj$time, levels = c("Uninjured", "1dpi", "3dpi", "7dpi"))

# Step 3: Extract UMAP Coordinates
umap_embeddings <- Embeddings(myeloid_seurat_obj, "umap")

# Ensure correct column names
umap_df <- data.frame(
  UMAP_1 = umap_embeddings[,1], 
  UMAP_2 = umap_embeddings[,2],
  time = myeloid_seurat_obj$time
)

# Step 4: Verify the Extracted Data
print("Checking extracted UMAP dataframe:")
print(head(umap_df))
print(table(umap_df$time))  # Ensure all time points are represented

# Step 5: Compute Density Function
compute_density <- function(df) {
  # Remove any NA or infinite values
  df <- df[is.finite(df$UMAP_1) & is.finite(df$UMAP_2), ]
  
  # Ensure at least 2 points exist before KDE computation
  if (nrow(df) < 2) {
    return(data.frame(UMAP_1 = numeric(0), UMAP_2 = numeric(0), density = numeric(0), time = unique(df$time)))  # Structured empty dataframe
  }
  
  # Compute KDE density
  dens <- kde2d(df$UMAP_1, df$UMAP_2, n = 100)
  dens_data <- data.frame(
    expand.grid(UMAP_1 = dens$x, UMAP_2 = dens$y),
    density = as.vector(dens$z)
  )
  return(dens_data)
}

# Step 6: Compute Density for Each Time Point
density_list <- lapply(levels(umap_df$time), function(tp) {
  tp_data <- subset(umap_df, time == tp)
  print(paste("Processing time point:", tp, "with", nrow(tp_data), "cells"))
  
  if (nrow(tp_data) < 2) {
    print(paste("Skipping", tp, "due to insufficient data"))
    return(NULL)
  }
  
  dens_df <- tryCatch({
    compute_density(tp_data)
  }, error = function(e) {
    print(paste("Error processing", tp, ":", e$message))
    return(NULL)
  })
  
  if (!is.null(dens_df) && nrow(dens_df) > 0) {
    dens_df$time <- tp
  }
  
  return(dens_df)
})

# Step 7: Remove NULL values before merging
density_list <- density_list[!sapply(density_list, is.null)]

# Ensure density_list is not empty before merging
if (length(density_list) == 0) {
  stop("ERROR: No density data generated. Check if UMAP embeddings exist.")
} else {
  density_df <- do.call(rbind, density_list)
}

# Step 8: Generate the UMAP Density Plot
# Open Cairo graphics device
Cairo::CairoPDF("umap_density_plot.pdf", width = 10, height = 7)
# Open Cairo graphics device (saving as PDF)
# Load required libraries
library(Cairo)
library(ggplot2)

# Save as high-resolution PNG to avoid transparency issues
CairoTIFF("umap_density_plot.tiff", width = 10, height = 7, units = "in", res = 300)# Generate the UMAP density plot
plot_figure <- ggplot() +
  geom_point(data = umap_df, aes(x = UMAP_1, y = UMAP_2), color = "gray80", alpha = 0.3) +  # Background points
  geom_tile(data = density_df, aes(x = UMAP_1, y = UMAP_2, fill = density), alpha = 0.8) +  # Density overlay
  scale_fill_gradient(low = "gray", high = "darkred", na.value = "gray80") +  # Handle NA densities
  facet_wrap(~ time, ncol = 2) +
  theme_minimal() +
  labs(title = "UMAP of Myeloid Cells Split by Time Point",
       fill = "Cell Density") +
  theme(strip.text = element_text(size = 14, face = "bold"))

# Print the plot
print(plot_figure)

# Close the graphics device to save the PNG file
dev.off()

# Also print the plot in the R session
plot_figure
