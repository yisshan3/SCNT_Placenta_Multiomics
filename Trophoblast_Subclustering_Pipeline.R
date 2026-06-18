# ==============================================================================
# Description: Sub-clustering of Trophoblast cells, Harmony & WNN integration, high-resolution annotation, DEGs (WT vs SCNT), and stat analysis.
# ==============================================================================


# Setup & Load Packages --------------------------------------------------------
set.seed(133)
setwd("/home4/ssyi/Mouse_Placenta")
rm(list = ls())

library(Seurat)
library(Signac)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(harmony)
library(SingleR)
library(celldex)
library(SummarizedExperiment)
library(scater)
library(clustree)
library(cowplot)
library(ggrastr)
library(data.table)


# Load Data & Metadata Processing ----------------------------------------------
integrated_obj <- readRDS("codes/2.Celltype_anno/integrated_obj.rds.gz")
meta <- integrated_obj@meta.data
tropho_cells <- rownames(meta[meta$Major_celltype == "Trophoblast", ])
tropho_obj <- subset(integrated_obj, cells = tropho_cells)

tropho_obj@meta.data$treatment <- gsub("\\d+", "", tropho_obj@meta.data$orig.ident)
tropho_obj@meta.data$treatment <- factor(tropho_obj@meta.data$treatment, levels = c("WT", "SCNT"))
tropho_obj@meta.data$time <- gsub("\\D+", "", tropho_obj@meta.data$orig.ident)
tropho_obj@meta.data$time <- factor(tropho_obj@meta.data$time, levels = c("125", "145", "185"))
tropho_obj@meta.data$Index <- rownames(tropho_obj@meta.data)


# Upstream Processing ----------------------------------------------------------
# RNA processing
DefaultAssay(tropho_obj) <- "RNA"
tropho_obj <- NormalizeData(tropho_obj) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE)     

# ATAC processing
DefaultAssay(tropho_obj) <- "ATAC"
tropho_obj <- RunTFIDF(tropho_obj)   
tropho_obj <- FindTopFeatures(tropho_obj, min.cutoff = 'q0')
tropho_obj <- RunSVD(tropho_obj)    

# Harmony RNA
DefaultAssay(tropho_obj) <- "RNA"
tropho_obj <- RunHarmony(
  object = tropho_obj,
  group.by.vars = 'orig.ident', reduction = 'pca', assay.use = 'RNA',
  project.dim = FALSE, reduction.save = "harmony_rna", plot_convergence = TRUE
)

# Harmony ATAC
DefaultAssay(tropho_obj) <- "ATAC"
tropho_obj <- RunHarmony(
  object = tropho_obj,
  group.by.vars = 'orig.ident', reduction = 'lsi', assay.use = 'ATAC',
  project.dim = FALSE, reduction.save = "harmony_atac", plot_convergence = TRUE
)

# WNN Integration
ElbowPlot(tropho_obj, ndims = 50, reduction = "harmony_rna")
ElbowPlot(tropho_obj, ndims = 50, reduction = "harmony_atac")

tropho_obj <- FindMultiModalNeighbors(
  object = tropho_obj,
  reduction.list = list("harmony_rna", "harmony_atac"),
  dims.list = list(1:50, 2:50), verbose = TRUE
)     

DefaultAssay(tropho_obj) <- "RNA"
tropho_obj <- RunUMAP(tropho_obj, reduction = 'harmony_rna', dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

DefaultAssay(tropho_obj) <- "ATAC"
tropho_obj <- RunUMAP(tropho_obj, reduction = 'harmony_atac', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

tropho_obj <- RunUMAP(tropho_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", verbose = TRUE, assay = "RNA")

tropho_obj <- FindClusters(tropho_obj, graph.name = 'wsnn', algorithm = 3, resolution = c(0.3, 0.5, 0.6, 0.8), verbose = FALSE)


# Subtype Annotation -----------------------------------------------------------
# 01.SingleR Reference-based Annotation
load("Public_data/elife/AllStages_TrophoblastNuclei_obj.Rdata")   
elife_tropho_count <- mouse.troph.combined@assays$RNA@counts
elife_tropho_cluster <- data.frame(Index = rownames(mouse.troph.combined@meta.data), cluster = mouse.troph.combined@active.ident)     
elife_tropho_cluster$Index <- NULL
elife_tropho_cluster <- elife_tropho_cluster[colnames(elife_tropho_count), , drop = FALSE]

elife_tropho_SE <- SummarizedExperiment(assays = list(counts = elife_tropho_count), colData = elife_tropho_cluster)   
elife_tropho_SE <- scater::logNormCounts(elife_tropho_SE) 

my_tropho_count <- tropho_obj@assays$RNA@counts
my_tropho_SE <- SummarizedExperiment(assays = list(counts = my_tropho_count))    
my_tropho_SE <- scater::logNormCounts(my_tropho_SE)

# Intersect common genes
common_gene <- intersect(rownames(my_tropho_SE), rownames(elife_tropho_SE))  
elife_tropho_SE <- elife_tropho_SE[common_gene, ]
my_tropho_SE <- my_tropho_SE[common_gene, ]

tropho_singleR_res <- SingleR(test = my_tropho_SE, ref = elife_tropho_SE, labels = elife_tropho_SE$X)

anno_df <- data.frame(Index = rownames(tropho_singleR_res), ref_label_from_elifeTropho = tropho_singleR_res$labels)
tropho_obj@meta.data <- tropho_obj@meta.data %>% inner_join(anno_df, by = "Index")      
rownames(tropho_obj@meta.data) <- tropho_obj@meta.data$Index

DimPlot(tropho_obj, reduction = "wnn.umap", group.by = 'ref_label_from_elifeTropho', label = TRUE, pt.size = 0.5) + NoLegend()

# Marker expression based on singleR annotation
cell_levels <- c("LaTP", "SynTII Precursor", "SynTII", "LaTP 2", "SynTI Precursor", "SynTI", "S-TGC Precursor", "S-TGC", "JZP 1", "JZP 2", "Glycogen Cells", "SpT Precursor", "SpT")
tropho_obj$ref_label_from_elifeTropho <- factor(tropho_obj$ref_label_from_elifeTropho, levels = cell_levels)
Idents(tropho_obj) <- factor(tropho_obj$ref_label_from_elifeTropho, levels = cell_levels)
DotPlot(tropho_obj, features = c("Met", "Epcam", "Ror2", "Lgr5", "Tcf7l1", "Gcm1", "Synb", "Egfr", "Pvt1", "Epha4", "Tgfa", "Eps8", "Tfrc", "Glis1", "Stra6", "Slc16a1", "Hand1", "Pparg", "Ctsq", "Ctsj", "Lepr", "Nos1ap", "Cdh4", "Prune2", "Ncam1", "Pcdh12", "Prl7b1", "Igfbp7", "Plac8", "Pla2g4d", "Mitf", "Flt1", "Prl8a9", "Slco2a1"), 
        group.by = "ref_label_from_elifeTropho") +
  scale_colour_gradient2(low = "gray", mid = "white", high = "firebrick3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 02.Manually Define Subtypes
DefaultAssay(tropho_obj) <- "RNA"
DimPlot(tropho_obj, reduction = "wnn.umap", group.by = 'wsnn_res.0.6', label = TRUE, pt.size = 0.5) + NoLegend()

DotPlot(tropho_obj, features = c("Met", "Epcam", "Egfr", "Ror2", "Lgr5", "Pvt1", "Tcf7l1", "Gcm1", "Synb", "Epha4", "Tgfa", "Eps8", "Tfrc", "Glis1", "Stra6", "Slc16a1", "Hand1", "Lifr", "Epas1", "Timp1", "Pparg", "Ctsq", "Ctsj", "Lepr", "Nos1ap", "Cdh4", "Prune2", "Ncam1", "Pcdh12", "Prl7b1", "Igfbp7", "Pla2g4d", "Plac8", "Mitf", "Flt1", "Prl8a9", "Slco2a1", "Ascl2", "Camk1d", "Cmss1", "Prl4a1", "Prl7d1", "Sfmbt2", "Prl3b1", "Prl8a8", "Gata4", "Pdpn", "Acta2", "Pdgfrb", "Col1a1"), 
        group.by = "wsnn_res.0.6") +
  scale_colour_gradient2(low = "gray", mid = "white", high = "firebrick3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Assign general cell types
tropho_obj_meta <- tropho_obj@meta.data %>% 
  mutate(General_celltype = case_when(
    wsnn_res.0.6 %in% c(0) ~ 'JZP1/2',
    wsnn_res.0.6 %in% c(1, 4, 13, 14) ~ 'GlyT',
    wsnn_res.0.6 %in% c(10) ~ 'Unknown1',
    wsnn_res.0.6 %in% c(15) ~ 'Unknown2',
    wsnn_res.0.6 %in% c(8) ~ 'SpT/Precursor',
    wsnn_res.0.6 %in% c(6) ~ 'SpT',
    wsnn_res.0.6 %in% c(11) ~ 'LaTP1/LaTP2/SynTII Precursor',
    wsnn_res.0.6 %in% c(9) ~ 'SynTII',
    wsnn_res.0.6 %in% c(12) ~ 'SynTI/Precursor',
    wsnn_res.0.6 %in% c(3, 5) ~ 'SynTI',
    wsnn_res.0.6 %in% c(2) ~ 'S-TGC/Precursor',
    wsnn_res.0.6 %in% c(7) ~ 'S-TGC'
  ))


# 03.Subset and Refine Ambiguous Clusters
# Subset: LaTP_SynT2pre
LaTP_SynT2pre <- subset(tropho_obj, wsnn_res.0.6 %in% c(11))
LaTP_SynT2pre <- FindMultiModalNeighbors(LaTP_SynT2pre, reduction.list = list("harmony_rna", "harmony_atac"), dims.list = list(1:50, 2:50), verbose = FALSE)
LaTP_SynT2pre <- FindClusters(LaTP_SynT2pre, graph.name = 'wsnn', algorithm = 3, resolution = c(0.3,0.4,0.5,0.6,0.8,1,1.2,1.5), verbose = FALSE)
DimPlot(LaTP_SynT2pre, reduction = "wnn.umap", label = TRUE, label.size = 3, group.by = "wsnn_res.1.5")
DotPlot(LaTP_SynT2pre, features = c("Met", "Ror2", "Lgr5",  "Epcam", "Pvt1", "Tcf7l1", "Egfr", "Gcm1", "Synb", "Vegfa", "Gcgr"), group.by = "wsnn_res.1.5") +
  scale_colour_gradient2(low = "gray", mid = "white", high = "firebrick3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
LaTP_SynT2pre_meta <- LaTP_SynT2pre@meta.data %>% 
  mutate(Tropho_celltype = case_when(
    wsnn_res.1.5 %in% c(0, 9) ~ 'LaTP',
    wsnn_res.1.5 == 8 ~ 'SynTII Precursor',
    wsnn_res.1.5 %in% c(1, 2, 3, 4, 5, 6, 7, 10) ~ 'LaTP2'
  ))

# Subset: JZP
JZP <- subset(tropho_obj, wsnn_res.0.6 %in% c(0))
JZP <- FindMultiModalNeighbors(JZP, reduction.list = list("harmony_rna", "harmony_atac"), dims.list = list(1:50, 2:50), verbose = FALSE)
JZP <- FindClusters(JZP, graph.name = 'wsnn', algorithm = 3, resolution = c(0.3,0.4,0.5,0.6,0.8,1,1.2,1.5), verbose = FALSE)
DimPlot(JZP, reduction = "wnn.umap", label = TRUE, label.size = 3, group.by = "wsnn_res.0.3")
DotPlot(JZP, features = c("Cdh4", "Pvt1", "Prune2", "Ncam1", "Pcdh12", "Prl7b1", "Igfbp7", "Pla2g4d", "Plac8"), group.by = "wsnn_res.0.3") +
  scale_colour_gradient2(low = "gray", mid = "white", high = "firebrick3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
JZP_meta <- JZP@meta.data %>% 
  mutate(Tropho_celltype = case_when(
    wsnn_res.0.3 %in% c(1, 2) ~ 'JZP1',
    wsnn_res.0.3 %in% c(0, 3) ~ 'JZP2'
  ))

# Subset: SpT_pre
SpT_pre <- subset(tropho_obj, wsnn_res.0.6 %in% c(8))
SpT_pre <- FindMultiModalNeighbors(SpT_pre, reduction.list = list("harmony_rna", "harmony_atac"), dims.list = list(1:50, 2:50), verbose = FALSE)
SpT_pre <- FindClusters(SpT_pre, graph.name = 'wsnn', algorithm = 3, resolution = c(0.3,0.4,0.5,0.6,0.8,1,1.2,1.5), verbose = FALSE)
DimPlot(SpT_pre, reduction = "wnn.umap", label = TRUE, label.size = 3, group.by = "wsnn_res.0.5")
DotPlot(SpT_pre, features = c("Cdh4", "Prune2", "Ncam1", "Plac8", "Mitf", "Flt1", "Prl8a9", "Slco2a1", "Ascl2"), group.by = "wsnn_res.0.5") +
  scale_colour_gradient2(low = "gray", mid = "white", high = "firebrick3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
SpT_pre_meta <- SpT_pre@meta.data %>% 
  mutate(Tropho_celltype = case_when(
    wsnn_res.0.5 %in% c(0, 4) ~ 'SpT Precursor',
    wsnn_res.0.5 %in% c(1, 2, 3) ~ 'SpT'
  ))

# Subset: TGC_pre
TGC_pre <- subset(tropho_obj, wsnn_res.0.6 %in% c(2))
TGC_pre <- FindMultiModalNeighbors(TGC_pre, reduction.list = list("harmony_rna", "harmony_atac"), dims.list = list(1:50, 2:50), verbose = FALSE)
TGC_pre <- FindClusters(TGC_pre, graph.name = 'wsnn', algorithm = 3, resolution = c(0.3,0.4,0.5,0.6,0.8,1,1.2,1.5), verbose = FALSE)
DimPlot(TGC_pre, reduction = "wnn.umap", label = TRUE, label.size = 3, group.by = "wsnn_res.0.5")
DotPlot(TGC_pre, features = c("Hand1", "Lifr", "Epas1", "Timp1", "Prl3d1", "Prl3b1", "Pparg", "Ctsq", "Ctsj", "Lepr", "Nos1ap"), group.by = "wsnn_res.0.5") +
  scale_colour_gradient2(low = "gray", mid = "white", high = "firebrick3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
TGC_pre_meta <- TGC_pre@meta.data %>% 
  mutate(Tropho_celltype = case_when(
    wsnn_res.0.5 %in% c(0, 3) ~ 'S-TGC Precursor',
    wsnn_res.0.5 %in% c(1, 2) ~ 'S-TGC'
  ))

# Subset: SynT1_pre
SynT1_pre <- subset(tropho_obj, wsnn_res.0.6 %in% c(12))
SynT1_pre <- FindMultiModalNeighbors(SynT1_pre, reduction.list = list("harmony_rna", "harmony_atac"), dims.list = list(1:50, 2:50), verbose = FALSE)
SynT1_pre <- FindClusters(SynT1_pre, graph.name = 'wsnn', algorithm = 3, resolution = c(0.3,0.4,0.5,0.6,0.8,1,1.2,1.5), verbose = FALSE)
DimPlot(SynT1_pre, reduction = "wnn.umap", label = TRUE, label.size = 3, group.by = "wsnn_res.0.3")
DotPlot(SynT1_pre, features = c("Pvt1", "Tcf7l1", "Gcm1", "Synb", "Epha4", "Tgfa", "Eps8", "Tfrc", "Glis1", "Stra6", "Slc16a1"), group.by = "wsnn_res.0.3") +
  scale_colour_gradient2(low = "gray", mid = "white", high = "firebrick3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
SynT1_pre_meta <- SynT1_pre@meta.data %>% 
  mutate(Tropho_celltype = case_when(
    wsnn_res.0.3 %in% c(0) ~ 'SynTI Precursor',
    wsnn_res.0.3 %in% c(1, 2) ~ 'SynTI'
  ))


# 04.Final Integration of Annotations
merge_sub_meta <- rbind(LaTP_SynT2pre_meta[, c("Index", "Tropho_celltype")], 
                        JZP_meta[, c("Index", "Tropho_celltype")], 
                        SpT_pre_meta[, c("Index", "Tropho_celltype")], 
                        TGC_pre_meta[, c("Index", "Tropho_celltype")], 
                        SynT1_pre_meta[, c("Index", "Tropho_celltype")])
tropho_obj_meta_final <- left_join(tropho_obj_meta, merge_sub_meta, by = "Index")
tropho_obj_meta_final[is.na(tropho_obj_meta_final$Tropho_celltype), "Tropho_celltype"] <- tropho_obj_meta_final[is.na(tropho_obj_meta_final$Tropho_celltype), "General_celltype"]
rownames(tropho_obj_meta_final) <- tropho_obj_meta_final$Index
tropho_obj@meta.data <- tropho_obj_meta_final

# Remove unknown clusters
non_Unknown_cells <- tropho_obj@meta.data[!(tropho_obj@meta.data$Tropho_celltype %in% c("Unknown1", "Unknown2")), "Index"]
tropho_obj <- subset(tropho_obj, cells = non_Unknown_cells)
dir.create("codes/3.Subtype_anno/", recursive = TRUE, showWarnings = FALSE)
saveRDS(tropho_obj, 'codes/3.Subtype_anno/tropho_obj.rds.gz', compress = 'gzip')

# Sub-population
WT_tropho_obj <- subset(tropho_obj, subset = treatment == "WT")
SCNT_tropho_obj <- subset(tropho_obj, subset = treatment == "SCNT")

WT_tropho_obj <- FindMultiModalNeighbors(WT_tropho_obj, reduction.list = list("harmony_rna", "harmony_atac"), dims.list = list(1:50, 2:50), verbose = TRUE)
WT_tropho_obj <- RunUMAP(WT_tropho_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", verbose = TRUE, assay = "RNA")
SCNT_tropho_obj <- FindMultiModalNeighbors(SCNT_tropho_obj, reduction.list = list("harmony_rna", "harmony_atac"), dims.list = list(1:50, 2:50), verbose = TRUE)
SCNT_tropho_obj <- RunUMAP(SCNT_tropho_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", verbose = TRUE, assay = "RNA")

saveRDS(WT_tropho_obj, "codes/3.Subtype_anno/WT_tropho_obj.rds.gz", compress = 'gzip')
saveRDS(SCNT_tropho_obj, "codes/3.Subtype_anno/SCNT_tropho_obj.rds.gz", compress = 'gzip')


# Plotting ---------------------------------------------------------------------
# 01.DimPlot
tropho_cell_levels <- c("LaTP", "LaTP2", "SynTII Precursor", "SynTII", "SynTI Precursor", "SynTI", "S-TGC Precursor", "S-TGC", "JZP1", "JZP2", "GlyT", "SpT Precursor", "SpT")
tropho_obj$Tropho_celltype <- factor(tropho_obj$Tropho_celltype, levels = tropho_cell_levels)

Tropho_cell_colors <- c("SynTI Precursor" = "#1f78b4", "SynTI" = "#a6cee3", "S-TGC Precursor" = "#33a02c", "S-TGC" = "#b2df8a", "SpT Precursor" = "#e31a1c", "SpT" = "#fb9a99",
                        "SynTII Precursor" = "#ff7f00", "SynTII" = "#fdbf6f", "JZP1" = "#6a3d9a", "JZP2" = "#cab2d6", "LaTP" = "#ffd92f", "LaTP2" = "#ffff99", "GlyT" = "#8da0cb")
treatment_colors <- c("WT" = "#A6CEE3", "SCNT" = "#FF9FCF")

p1 <- DimPlot(tropho_obj, reduction = "wnn.umap", label = FALSE, cols = Tropho_cell_colors) + theme(aspect.ratio = 1)
p2 <- DimPlot(tropho_obj, reduction = "wnn.umap", group.by = 'Tropho_celltype', split.by = 'treatment', label = FALSE, cols = Tropho_cell_colors) + theme(aspect.ratio = 1)
p3 <- DimPlot(tropho_obj, reduction = "wnn.umap", group.by = 'treatment', label = FALSE, cols = treatment_colors) + theme(aspect.ratio = 1)
ggsave("Figures/Tropho_wnnUMAP.pdf", p1, width = 6.5, height = 5)
ggsave("Figures/Tropho_wnnUMAP_samples.pdf", p2, width = 6.5, height = 5)
ggsave("Figures/Tropho_wnnUMAP_treatment.pdf", p3, width = 10, height = 5)


# 02.DotPlot
Idents(tropho_obj) <- factor(tropho_obj$Tropho_celltype, levels = tropho_cell_levels)

p4 <- DotPlot(tropho_obj, features = c("Met", "Epcam", "Ror2", "Lgr5", "Tcf7l1", "Gcm1", "Slc16a3", "Synb", "Egfr", "Pvt1", "Syna", "Epha4", "Tgfa", "Eps8", "Tfrc", "Glis1", "Slc16a1", "Stra6", "Hand1", "Pparg", "Nos1ap", "Lepr", "Ctsj", "Ctsq", 
                                       "Cdh4", "Gjb3", "Prune2", "Ncam1", "Tpbpa", "Plac8", "Pcdh12", "Igfbp7", "Pla2g4d", "Prl7b1", "Flt1", "Mitf", "Prl3b1", "Prl8a9", "Prl8a8", "Slco2a1"), 
              group.by = "Tropho_celltype") +
  scale_colour_gradient2(low = "#1f78b4", mid = "white", high = "firebrick3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.3)
ggsave("Figures/Tropho_marker_Dotplot.pdf", p4, width = 14, height = 4.5)


# 03.FeaturePlot 
p5 <- FeaturePlot(tropho_obj, features = c("Epcam", "Lgr5", "Gcm1", "Nos1ap", "Ctsq", "Tfrc", "Tpbpa", "Pcdh12", "Flt1", "Prl8a8"), cols = c("lightgrey", "firebrick3"), reduction = "wnn.umap", ncol = 5) 
p5 <- rasterize(p5, dpi = 300)
ggsave("Figures/Tropho_marker_Featureplot.pdf", p5, width = 3.5 * 5, height = 3 * 2)

p6 <- FeaturePlot(tropho_obj, features = c("Mki67"), cols = c("lightgrey", "firebrick3"), reduction = "wnn.umap")
p6 <- rasterize(p6, dpi = 300)
ggsave("Figures/Tropho_Mki67_Featureplot.pdf", p6, width = 3.5, height = 3)

p7 <- FeaturePlot(tropho_obj, features = c("Gjb3"), cols = c("lightgrey", "firebrick3"), reduction = "wnn.umap")
p7 <- rasterize(p7, dpi = 300)
ggsave("Figures/Tropho_Gjb3_Featureplot.pdf", p7, width = 3.5, height = 3)


# Differentially Expressed Genes (DEGs) ----------------------------------------
tropho_obj$treatment_celltype <- paste(tropho_obj$treatment, tropho_obj$Tropho_celltype, sep = "_")
Idents(tropho_obj) <- tropho_obj$treatment_celltype

tropho_celltypes <- unique(tropho_obj$Tropho_celltype)

for(i in tropho_celltypes) {
  ident.1_SCNT <- paste("SCNT", i, sep = "_")
  ident.2_WT <- paste("WT", i, sep = "_")
  
  # Calculate DEGs (SCNT vs WT) per cell type
  celltype_DEGs <- FindMarkers(tropho_obj, ident.1 = ident.1_SCNT, ident.2 = ident.2_WT, 
                               only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  celltype_DEGs$gene <- rownames(celltype_DEGs)
  
  celltype_DEGs$label <- ifelse(celltype_DEGs$avg_log2FC > 0, 'SigUp', 'SigDown')
  celltype_DEGs$q_val_fdr <- p.adjust(celltype_DEGs$p_val, method = "BH") 
  
  celltype_DEGs <- celltype_DEGs %>% 
    dplyr::filter(q_val_fdr < 0.05) %>%
    dplyr::mutate(label = factor(label, levels = c("SigUp", "SigDown"))) %>%
    dplyr::arrange(label, q_val_fdr)
  
  colnames(celltype_DEGs)[colnames(celltype_DEGs) == "p_val_adj"] <- "p_val_Bonferroni"
  colnames(celltype_DEGs)[colnames(celltype_DEGs) == "q_val_fdr"] <- "p_val_BH(FDR)"
  
  celltype_DEGs <- celltype_DEGs[, c("p_val", "avg_log2FC", "pct.1", "pct.2", 
                                     "p_val_Bonferroni", "p_val_BH(FDR)", "gene", "label")]
  x <- gsub(" ", "_", i)
  output_path <- paste0('codes/3.Subtype_anno/Trophoblast/', x, '_DEGs.csv')
  write.csv(celltype_DEGs, output_path, row.names = TRUE, quote = FALSE)
}


# Statistical Proportions & Dynamics -------------------------------------------
# 01.Cell Type Proportion by Time & Group
cell_pct <- data.frame(table(tropho_obj$orig.ident, tropho_obj$Tropho_celltype))
colnames(cell_pct) <- c("Sample", "Celltype", "Num") 
cell_pct <- spread(cell_pct, "Celltype", "Num")
cell_pct[, 2:ncol(cell_pct)] <- 100 * cell_pct[, 2:ncol(cell_pct)] / rowSums(cell_pct[, 2:ncol(cell_pct)])
cell_pct <- gather(cell_pct, "Celltype", "Percentage", -"Sample")

# extraction of metadata identifiers
cell_pct$Group <- gsub("\\d+", "", cell_pct$Sample)
cell_pct$Time <- gsub("\\D+", "", cell_pct$Sample)
cell_pct$Group <- factor(cell_pct$Group, levels = c("WT", "SCNT"))
cell_pct$Celltype <- factor(cell_pct$Celltype, levels = tropho_cell_levels)

p1 <- ggbarplot(cell_pct, x = "Time", y = "Percentage", fill = "Celltype", color = "Celltype",
                palette = Tropho_cell_colors, label = FALSE, facet.by = "Group") +
  ylab("Percentage (%)") + theme(aspect.ratio = 1.9) 
ggsave("Figures/Tropho_stat_celltype.pdf", p1, width = 9, height = 6)


# 02.Precursor vs Differentiated Proportion
tropho_obj$Cell_stage <- ifelse(
  tropho_obj$Tropho_celltype %in% c("SynTII", "SynTI", "S-TGC", "GlyT", "SpT"), 
  "Differentiated", "Precursor"
)
cell_pct <- data.frame(table(tropho_obj$orig.ident, tropho_obj$Cell_stage))
colnames(cell_pct) <- c("Sample", "Cell_stage", "Num") 
cell_pct <- spread(cell_pct, "Cell_stage", "Num")
cell_pct[, 2:ncol(cell_pct)] <- 100 * cell_pct[, 2:ncol(cell_pct)] / rowSums(cell_pct[, 2:ncol(cell_pct)])
cell_pct <- gather(cell_pct, "Cell_stage", "Percentage", -"Sample")

cell_pct$Sample <- factor(cell_pct$Sample, levels = c("SCNT185", "WT185", "SCNT145", "WT145", "SCNT125", "WT125"))
cell_pct$Cell_stage <- factor(cell_pct$Cell_stage, levels = c("Precursor", "Differentiated"))

Precursor_Diff_cell_colors <- c("Precursor" = "#a1d76a", "Differentiated" = "#e9a3c9")

p3 <- ggbarplot(cell_pct, x = "Sample", y = "Percentage", fill = "Cell_stage", color = "Cell_stage",
                palette = Precursor_Diff_cell_colors, label = FALSE) +
  ylab("Percentage (%)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.7) + coord_flip()
ggsave("Figures/Tropho_Precursor_Diff_stat.pdf", p3, width = 5, height = 4)


# 03.Sample Contribution per Cell Type
cell_pct <- data.frame(table(tropho_obj$orig.ident, tropho_obj$Tropho_celltype))
colnames(cell_pct) <- c("Sample", "Celltype", "Num") 
sample_pct <- spread(cell_pct, "Sample", "Num") %>% column_to_rownames(var = "Celltype")
sample_pct["col_sum", ] <- colSums(sample_pct)
max_sum <- max(sample_pct["col_sum", ])

scale_sample_pct <- apply(sample_pct, 2, function(col) { col * max_sum / col[nrow(sample_pct)] }) %>% as.data.frame()
scale_sample_pct <- scale_sample_pct[1:(nrow(sample_pct) - 1), ] %>% rownames_to_column(var = "Celltype")
scale_sample_pct[, 2:7] <- 100 * scale_sample_pct[, 2:7] / rowSums(scale_sample_pct[, 2:7])
scale_sample_pct <- gather(scale_sample_pct, "Sample", "Percentage", -"Celltype")
scale_sample_pct$Group <- gsub("\\d+", "", scale_sample_pct$Sample)
scale_sample_pct$Time <- gsub("\\D+", "", scale_sample_pct$Sample)

scale_sample_pct$Group <- factor(scale_sample_pct$Group, levels = c("SCNT", "WT"))
scale_sample_pct$Sample <- factor(scale_sample_pct$Sample, levels = c("SCNT185", "SCNT145", "SCNT125", "WT185", "WT145", "WT125"))
scale_sample_pct$Celltype <- factor(scale_sample_pct$Celltype, levels = c("LaTP", "LaTP2", "JZP1", "JZP2", "SynTI Precursor", "S-TGC Precursor", "SynTII Precursor", "SpT Precursor", "GlyT", "S-TGC", "SynTI", "SynTII", "SpT"))

sample_colors <- c("WT125" = "#d1e5f0", "WT145" = "#67a9cf", "WT185" = "#2166ac", "SCNT125" = "#fddbc7", "SCNT145" = "#ef8a62", "SCNT185" = "#b2182b")

p1 <- ggbarplot(scale_sample_pct, x = "Celltype", y = "Percentage", fill = "Sample", color = "Sample",
                palette = sample_colors, label = FALSE) +
  ylab("Percentage (%)") + theme(aspect.ratio = 0.8) + coord_flip()
ggsave("Figures/Tropho_stat_sample.pdf", p1, width = 6, height = 7)


# 04.Trophoblast Dynamics Lineplots
md <- as.data.table(tropho_obj@meta.data)
md <- md[, .(sample_id = orig.ident, condition = treatment, celltype = Tropho_celltype)]

count_dt <- md[, .N, by = .(sample_id, condition, celltype)]
setnames(count_dt, "N", "n_cells")

grid <- CJ(sample_id = unique(count_dt$sample_id), celltype = unique(count_dt$celltype))
sample_cond <- unique(count_dt[, .(sample_id, condition)])
grid <- merge(grid, sample_cond, by = "sample_id", all.x = TRUE)

count_dt2 <- merge(grid, count_dt, by = c("sample_id", "condition", "celltype"), all.x = TRUE)
count_dt2[is.na(n_cells), n_cells := 0]
count_dt2[, total_cells := sum(n_cells), by = sample_id]
count_dt2[, prop := fifelse(total_cells > 0, n_cells / total_cells, NA_real_)]

dt <- copy(count_dt2)
dt[, time := as.integer(gsub("\\D+", "", sample_id))]
dt[, condition := factor(condition, levels = c("WT", "SCNT"))]
dt[, time := factor(time, levels = sort(unique(time)))]
dt$celltype <- factor(dt$celltype, levels = tropho_cell_levels)

p1 <- ggplot(dt, aes(x = time, y = prop, group = condition, color = condition)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~ celltype, ncol = 4, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = NULL)) +
  scale_color_manual(values = c(WT = "#A6CEE3", SCNT = "#FF9FCF")) +
  labs(x = "Time", y = "Cell proportion", color = "Condition") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    aspect.ratio = 0.8
  )
ggsave("Figures/Tropho_celltype_dev_dynamics_lineplots.pdf", p1, width = 9, height = 7)
