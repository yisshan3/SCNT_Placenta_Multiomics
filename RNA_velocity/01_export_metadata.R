# ==============================================================================
# Description: Extract and save cell embeddings, barcodes, and cell type metadata from the Seurat object for downstream RNA velocity analysis (scVelo).
# ==============================================================================


# Setup & Load Packages --------------------------------------------------------
set.seed(133)
setwd("/home4/ssyi/Mouse_Placenta")
rm(list = ls())

library(Seurat)
library(Signac)
library(tidyverse)


# Define Export Function -------------------------------------------------------
export_scvelo_meta <- function(seurat_obj, group_name) {
  # UMAP coordinates
  write.csv(Embeddings(seurat_obj, reduction = "wnn.umap"), file = file.path(out_dir, paste0(group_name, "_cell_embeddings.csv")))
  
  # cell barcodes
  write.csv(Cells(seurat_obj), file = file.path(out_dir, paste0(group_name, "_cellID_obs.csv")), row.names = FALSE)
  
  # cell type annotations
  write.csv(seurat_obj@meta.data[, 'Tropho_celltype', drop = FALSE], file = file.path(out_dir, paste0(group_name, "_cell_celltype.csv")))
}


# Load Data & Execute Export ---------------------------------------------------
out_dir <- "codes/4.RNA_velocity"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Export WT group data
WT_tropho_obj <- readRDS("codes/3.Subtype_anno/WT_tropho_obj.rds.gz")
export_scvelo_meta(WT_tropho_obj, "WT")

# Export SCNT group data
SCNT_tropho_obj <- readRDS("codes/3.Subtype_anno/SCNT_tropho_obj.rds.gz")
export_scvelo_meta(SCNT_tropho_obj, "SCNT")

