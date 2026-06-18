# ==============================================================================
# Description: Monocle 3 (using WNN UMAP for learn_graph), Monocle 2 (DDRTree), and trajectory-associated DEG & GO enrichment analysis.
# ==============================================================================


# Setup & Load Packages --------------------------------------------------------
set.seed(133)
setwd("/home4/ssyi/Mouse_Placenta")
rm(list = ls())

library(Seurat)
library(Signac)
library(tidyverse)
library(dplyr)
library(monocle3)
library(monocle)
library(patchwork) 
library(clusterProfiler)
library(org.Mm.eg.db)


# Load Data & Define Variables -------------------------------------------------
WT_tropho_obj <- readRDS("codes/3.Subtype_anno/WT_tropho_obj.rds.gz")
SCNT_tropho_obj <- readRDS("codes/3.Subtype_anno/SCNT_tropho_obj.rds.gz")

# Define levels and colors
tropho_cell_levels <- c("LaTP", "LaTP2", "SynTII Precursor", "SynTII", "SynTI Precursor", "SynTI", "S-TGC Precursor", "S-TGC", "JZP1", "JZP2", "GlyT", "SpT Precursor", "SpT")
Tropho_cell_colors <- c("SynTI Precursor" = "#1f78b4", "SynTI" = "#a6cee3", "S-TGC Precursor" = "#33a02c", "S-TGC" = "#b2df8a", "SpT Precursor" = "#e31a1c", "SpT" = "#fb9a99",
                        "SynTII Precursor" = "#ff7f00", "SynTII" = "#fdbf6f", "JZP1" = "#6a3d9a", "JZP2" = "#cab2d6", "LaTP" = "#ffd92f", "LaTP2" = "#ffff99", "GlyT" = "#8da0cb")


# Monocle3 Analysis ------------------------------------------------------------
# Define Monocle3 processing function
Run_Monocle3 <- function(Seurat_Obj, num_k = 20) {
  obj <- Seurat_Obj
  data <- obj@assays$RNA@counts
  cell_metadata <- obj@meta.data 
  gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
  
  cds <- new_cell_data_set(
    data,
    cell_metadata = cell_metadata,
    gene_metadata = gene_annotation
  )
  
  cds <- preprocess_cds(cds, num_dim = 50)   
  cds <- reduce_dimension(cds, preprocess_method = "PCA")  
  
  # Import WNN UMAP coordinates from Seurat
  cds.embed <- cds@int_colData$reducedDims$UMAP
  int.embed <- Embeddings(obj, reduction = "wnn.umap")   
  int.embed <- int.embed[rownames(cds.embed), ]
  cds@int_colData$reducedDims$UMAP <- int.embed
  
  cds <- cluster_cells(cds, k = num_k)
  cds <- learn_graph(cds)
  root.cell <- cds@colData[cds@colData$Tropho_celltype %in% c("LaTP", "JZP1", "SpT Precursor"), "Index"]   
  cds <- order_cells(cds, root_cells = root.cell)
  
  return(cds)
}

# Run Monocle3
WT_tropho_monocle3_cds <- Run_Monocle3(WT_tropho_obj, num_k = 20)
SCNT_tropho_monocle3_cds <- Run_Monocle3(SCNT_tropho_obj, num_k = 20)

plot_Monocle3 <- function(cds) {
  p1 <- plot_cells(cds, 
                   color_cells_by = "Tropho_celltype",
                   label_groups_by_cluster = FALSE,
                   label_leaves = FALSE,
                   label_branch_points = TRUE,
                   label_cell_groups = TRUE,
                   label_roots = FALSE,
                   group_label_size = 3) +
    scale_color_manual(values = Tropho_cell_colors)
  
  p2 <- plot_cells(cds,
                   color_cells_by = "pseudotime",
                   label_cell_groups = FALSE,
                   label_groups_by_cluster = FALSE,
                   label_leaves = FALSE,
                   label_branch_points = FALSE,
                   label_roots = FALSE)
  
  return(p1 | p2)
}

# Visualize trajectories
p1 <- plot_Monocle3(WT_tropho_monocle3_cds)
p2 <- plot_Monocle3(SCNT_tropho_monocle3_cds)
ggsave("Figures/WT_tropho_monocle3.pdf", p1, width = 4.5 * 2, height = 4)
ggsave("Figures/SCNT_tropho_monocle3.pdf", p2, width = 4.5 * 2, height = 4)


# Pseudotime Distribution Violin Plots -----------------------------------------
# WT Plot
plot_WT_dta <- tibble(
  barcode = rownames(colData(WT_tropho_monocle3_cds)),
  orig_ident = colData(WT_tropho_monocle3_cds)$orig.ident,
  Celltype = colData(WT_tropho_monocle3_cds)$Tropho_celltype,
  pseudotime = pseudotime(WT_tropho_monocle3_cds)
)

plot_WT_dta <- plot_WT_dta %>% mutate(
  region = case_when(
    Celltype %in% c("JZP1", "JZP2", "GlyT", "SpT Precursor", "SpT") ~ "JZ",
    Celltype %in% c("LaTP", "LaTP2", "SynTII Precursor", "SynTII", "SynTI Precursor", "SynTI", "S-TGC Precursor", "S-TGC") ~ "Lab"
  )
)

plot_WT_dta$Celltype <- factor(plot_WT_dta$Celltype, levels = tropho_cell_levels)
p3 <- ggplot(plot_WT_dta, aes(x = Celltype, y = pseudotime)) +
  geom_violin(aes(fill = Celltype), scale = "width", trim = TRUE) +
  geom_boxplot(aes(fill = Celltype), width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA, lwd = 0.3) +
  coord_cartesian(ylim = c(0, 20)) +
  scale_fill_manual(values = Tropho_cell_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~ region, scales = "free", space = "free") 
ggsave("Figures/WT_tropho_monocle3_pseudotime_vlnplot.pdf", p3, width = 8, height = 5)

# SCNT Plot
plot_SCNT_dta <- tibble(
  barcode = rownames(colData(SCNT_tropho_monocle3_cds)),
  orig_ident = colData(SCNT_tropho_monocle3_cds)$orig.ident,
  Celltype = colData(SCNT_tropho_monocle3_cds)$Tropho_celltype,
  pseudotime = pseudotime(SCNT_tropho_monocle3_cds)
)

plot_SCNT_dta <- plot_SCNT_dta %>% mutate(
  region = case_when(
    Celltype %in% c("JZP1", "JZP2", "GlyT", "SpT Precursor", "SpT") ~ "JZ",
    Celltype %in% c("LaTP", "LaTP2", "SynTII Precursor", "SynTII", "SynTI Precursor", "SynTI", "S-TGC Precursor", "S-TGC") ~ "Lab"
  )
)

plot_SCNT_dta$Celltype <- factor(plot_SCNT_dta$Celltype, levels = tropho_cell_levels)
p4 <- ggplot(plot_SCNT_dta, aes(x = Celltype, y = pseudotime)) +
  geom_violin(aes(fill = Celltype), scale = "width", trim = TRUE) +
  geom_boxplot(aes(fill = Celltype), width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA, lwd = 0.3) +
  coord_cartesian(ylim = c(0, 15)) +
  scale_fill_manual(values = Tropho_cell_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~ region, scales = "free", space = "free") 
ggsave("Figures/SCNT_tropho_monocle3_pseudotime_vlnplot.pdf", p4, width = 8, height = 5)


# Trajectory Gene Analysis & GO Enrichment -------------------------------------
WT_JZ_cds <- WT_tropho_monocle3_cds[, colData(WT_tropho_monocle3_cds) %>%
                                      subset(Tropho_celltype %in% c("JZP1", "JZP2", "GlyT", "SpT Precursor", "SpT")) %>%
                                      row.names]

SCNT_JZ_cds <- SCNT_tropho_monocle3_cds[, colData(SCNT_tropho_monocle3_cds) %>%
                                          subset(Tropho_celltype %in% c("JZP1", "JZP2", "GlyT", "SpT Precursor", "SpT")) %>%
                                          row.names]

WT_JZ_Track_genes <- graph_test(WT_JZ_cds, neighbor_graph = "principal_graph", cores = 10)
SCNT_JZ_Track_genes <- graph_test(SCNT_JZ_cds, neighbor_graph = "principal_graph", cores = 10)

filter_WT_JZ_Track_genes <- WT_JZ_Track_genes %>%
  as_tibble() %>%
  dplyr::filter(q_value < 0.05) %>%
  dplyr::filter(abs(morans_I) >= 0.3)

filter_SCNT_JZ_Track_genes <- SCNT_JZ_Track_genes %>%
  as_tibble() %>%
  dplyr::filter(q_value < 0.05) %>%
  dplyr::filter(abs(morans_I) >= 0.3)

WT_JZ_specific_track_genes <- setdiff(filter_WT_JZ_Track_genes$gene_short_name, filter_SCNT_JZ_Track_genes$gene_short_name)   
SCNT_JZ_specific_track_genes <- setdiff(filter_SCNT_JZ_Track_genes$gene_short_name, filter_WT_JZ_Track_genes$gene_short_name)  
write_tsv(as.data.frame(SCNT_JZ_specific_track_genes), "codes/5.Monocle/SCNT_JZ_diff_Track_133genes.txt", col_names = FALSE)

# GO Enrichment
enrich_go <- enrichGO(
  gene = SCNT_JZ_specific_track_genes,
  OrgDb = org.Mm.eg.db,
  keyType = 'SYMBOL',
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
go_data <- as_tibble(enrich_go@result)
go_data$GeneRatio <- as.numeric(str_split_i(go_data$GeneRatio, '/', 1)) / as.numeric(str_split_i(go_data$GeneRatio, '/', 2))
go_data$BgRatio <- as.numeric(str_split_i(go_data$BgRatio, '/', 1)) / as.numeric(str_split_i(go_data$BgRatio, '/', 2))
go_data$fold_enrichment <- go_data$GeneRatio / go_data$BgRatio
go_data$log10_p <- -log10((go_data$pvalue))
go_data_top10 <- go_data %>% arrange(desc(log10_p)) %>% head(10)
go_data_plot <- rbind(
  go_data_top10[!go_data_top10$Description %in% c('nuclear chromosome segregation', 'chromosome segregation', 'positive regulation of cell cycle process', 'regulation of mitotic cell cycle phase transition'), ],
  go_data[go_data$Description %in% c('mitotic G2 DNA damage checkpoint signaling', 'double-strand break repair', 'G2/M transition of mitotic cell cycle', 'double-strand break repair via homologous recombination'), ]
)

p5 <- ggplot(data = arrange(go_data_plot, desc(log10_p)),
             aes(x = reorder(Description, log10_p), y = log10_p, fill = fold_enrichment)) +
  geom_bar(stat = 'identity', width = 0.9) +
  theme_classic() +
  geom_text(aes(label = Description), y = 0, hjust = 0, color = "black", size = 5) +
  scale_fill_gradient(name = "Fold enrichment", low = "#fee5d9", high = "#cb181d") +
  theme(axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15),
        legend.text = element_text(colour = "black", size = 13),
        legend.title = element_text(colour = "black", size = 15)) +
  coord_flip() +
  theme(aspect.ratio = 0.7)
ggsave("Figures/SCNT_JZ_specific_Track_genes_clusterprofiler_GOBP_barplot.pdf", p5, width = 12, height = 5)


# Monocle2 Analysis ------------------------------------------------------------
# Define Monocle2 processing function
Run_Monocle2 <- function(Seurat_obj) {
  expr_matirx <- as(as.matrix(Seurat_obj@assays$RNA@counts), 'sparseMatrix')    
  p_data <- Seurat_obj@meta.data
  f_Data <- data.frame(gene_short_name = row.names(Seurat_obj), row.names = row.names(Seurat_obj))
  
  pd <- new('AnnotatedDataFrame', data = p_data)
  fd <- new('AnnotatedDataFrame', data = f_Data)
  
  monocle_cds <- newCellDataSet(
    expr_matirx,
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )  
  
  # Estimate size factors and dispersions
  monocle_cds <- estimateSizeFactors(monocle_cds)   
  monocle_cds <- estimateDispersions(monocle_cds)   
  
  # Filter low quality genes
  monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)   
  expressed_genes <- row.names(subset(fData(monocle_cds), num_cells_expressed >= 10))   
  
  # Differential gene test for ordering
  diff_test_res <- differentialGeneTest(
    monocle_cds[expressed_genes, ],
    fullModelFormulaStr = "~Tropho_celltype",
    cores = 6
  )
  
  degs <- subset(diff_test_res, qval < 0.01)
  n_ordering_genes <- min(2000, nrow(degs))
  ordering_genes <- rownames(degs)[order(degs$qval)][1:n_ordering_genes]
  
  monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
  monocle_cds <- reduceDimension(
    monocle_cds,
    max_components = 3,
    num_dim = 20,
    method = 'DDRTree'
  )   
  monocle_cds <- orderCells(monocle_cds)
  
  return(monocle_cds)
}

# Run Monocle2
WT_tropho_monocle2_cds <- Run_Monocle2(WT_tropho_obj)
SCNT_tropho_monocle2_cds <- Run_Monocle2(SCNT_tropho_obj)

plot_Monocle2 <- function(monocle_cds) {
  monocle_cds$Tropho_celltype <- factor(monocle_cds$Tropho_celltype, levels = tropho_cell_levels)
  
  p_traj1 <- plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 1, show_backbone = TRUE) +
    scale_color_gradientn(colours = c("#2c7bb6", "#abd9e9", "#d7191c"))
  
  p_traj2 <- plot_cell_trajectory(monocle_cds, color_by = "Tropho_celltype", cell_size = 1, show_backbone = TRUE) +
    theme(legend.position = 'none', panel.border = element_blank()) +
    scale_color_manual(values = Tropho_cell_colors)
  
  segment_layers <- which(sapply(p_traj2$layers, function(x) inherits(x$geom, "GeomSegment")))
  point_layers <- which(sapply(p_traj2$layers, function(x) inherits(x$geom, "GeomPoint")))
  text_layers <- which(sapply(p_traj2$layers, function(x) inherits(x$geom, "GeomText")))
  new_layers <- c(p_traj2$layers[point_layers], p_traj2$layers[segment_layers], p_traj2$layers[text_layers])
  p_traj2$layers <- new_layers
  
  return(p_traj1 | p_traj2)
}

# Visualize DDRTree trajectories
p6 <- plot_Monocle2(WT_tropho_monocle2_cds)
p7 <- plot_Monocle2(SCNT_tropho_monocle2_cds)
ggsave("Figures/WT_tropho_monocle2.pdf", p6, width = 4.5 * 2, height = 4)
ggsave("Figures/SCNT_tropho_monocle2.pdf", p7, width = 4.5 * 2, height = 4)
