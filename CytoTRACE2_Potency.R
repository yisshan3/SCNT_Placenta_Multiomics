# ==============================================================================
# Description: Run CytoTRACE2 to get Potency Score, align coordinates with WNN UMAP, and visualize correlation with Cell Cycle phases.
# ==============================================================================


# Setup & Load Packages --------------------------------------------------------
set.seed(133)
setwd("/home4/ssyi/Mouse_Placenta")
rm(list = ls())

library(CytoTRACE2)
library(Seurat)
library(Signac)
library(tidyverse)
library(dplyr)
library(gtools)
library(ggplot2)
library(patchwork)


# Load Data --------------------------------------------------------------------
tropho_obj <- readRDS("codes/3.Subtype_anno/tropho_obj.rds.gz")


# Run CytoTRACE2 ---------------------------------------------------------------
cytotrace2_result_sce <- cytotrace2(
  tropho_obj, 
  is_seurat = TRUE, 
  slot_type = "counts", 
  species = 'mouse',
  seed = 133
)
dir.create("codes/16.CytoTRACE2/", recursive = TRUE, showWarnings = FALSE)
saveRDS(cytotrace2_result_sce, "codes/16.CytoTRACE2/Tropho_cytotrace2_result.rds.gz", compress = "gzip")

# Generate annotation plots
annotation <- data.frame(phenotype = tropho_obj@meta.data$Tropho_celltype) %>% 
  set_rownames(., colnames(tropho_obj))
plots <- plotData(
  cytotrace2_result = cytotrace2_result_sce, 
  annotation = annotation, 
  is_seurat = TRUE
)
saveRDS(plots, "codes/16.CytoTRACE2/Tropho_cytotrace2_plots.rds.gz", compress = "gzip")


# Align CytoTRACE2 with WNN UMAP -----------------------------------------------
umap_raw <- as.data.frame(tropho_obj@reductions$wnn.umap@cell.embeddings) 

plot_names <- c("CytoTRACE2_UMAP", "CytoTRACE2_Potency_UMAP", 
                "CytoTRACE2_Relative_UMAP", "Phenotype_UMAP")

# Loop through each plot and update coordinates to match WNN UMAP
for (plot_name in plot_names) {
  if (!is.null(plots[[plot_name]][[1]]$data)) {
    plots[[plot_name]][[1]]$data$umap_1 <- umap_raw$wnnUMAP_1
    plots[[plot_name]][[1]]$data$umap_2 <- umap_raw$wnnUMAP_2
  }
}

# Determine UMAP ranges
x_limits <- range(umap_raw$wnnUMAP_1, na.rm = TRUE)
y_limits <- range(umap_raw$wnnUMAP_2, na.rm = TRUE)

plots$CytoTRACE2_UMAP[[1]] <- plots$CytoTRACE2_UMAP[[1]] + 
  coord_cartesian(xlim = x_limits, ylim = y_limits)
p1_ct <- plots$CytoTRACE2_UMAP

plots$CytoTRACE2_Potency_UMAP[[1]] <- plots$CytoTRACE2_Potency_UMAP[[1]] + 
  coord_cartesian(xlim = x_limits, ylim = y_limits)
p2_ct <- plots$CytoTRACE2_Potency_UMAP

plots$CytoTRACE2_Relative_UMAP[[1]] <- plots$CytoTRACE2_Relative_UMAP[[1]] + 
  coord_cartesian(xlim = x_limits, ylim = y_limits)
p3_ct <- plots$CytoTRACE2_Relative_UMAP

plots$Phenotype_UMAP[[1]] <- plots$Phenotype_UMAP[[1]] + 
  coord_cartesian(xlim = x_limits, ylim = y_limits)
p4_ct <- plots$Phenotype_UMAP

combined_plot <- (p1_ct[[1]] | p2_ct[[1]]) / (p3_ct[[1]] | p4_ct[[1]])
ggsave("Figures/Tropho_CytoTRACE2_combined_plot.pdf", plot = combined_plot, width = 10, height = 10)


# Potency Score FeaturePlots (WT vs SCNT) --------------------------------------
p1 <- FeaturePlot(subset(cytotrace2_result_sce, treatment == 'WT'), 
                  features = "CytoTRACE2_Score", 
                  reduction = "wnn.umap") +
  scale_colour_gradientn(
    colours = c("#5E4FA2", "#66C2A5", "#E6F598", "#FEE08B", "#F46D43", "#9E0142"),
    na.value = "transparent", 
    limits = c(0, 1), 
    breaks = seq(0, 1, by = 0.2), 
    labels = c("0.0 (More diff.)", "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
    name = "Potency\nScore \n",
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  )

p2 <- FeaturePlot(subset(cytotrace2_result_sce, treatment == 'SCNT'),
                  features = "CytoTRACE2_Score",
                  reduction = "wnn.umap") +
  scale_colour_gradientn(
    colours = c("#5E4FA2", "#66C2A5", "#E6F598", "#FEE08B", "#F46D43", "#9E0142"),
    na.value = "transparent", 
    limits = c(0, 1), 
    breaks = seq(0, 1, by = 0.2), 
    labels = c("0.0 (More diff.)", "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
    name = "Potency\nScore \n",
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  )

ggsave("Figures/Tropho_CytoTRACE2_Score_featureplot_WT.pdf", plot = p1, width = 6, height = 4.5)
ggsave("Figures/Tropho_CytoTRACE2_Score_featureplot_SCNT.pdf", plot = p2, width = 6, height = 4.5)


# Potency Statistics & Plotting ------------------------------------------------
# 01.Potency Proportion per Sample
cell_pct <- data.frame(table(cytotrace2_result_sce$orig.ident, cytotrace2_result_sce$CytoTRACE2_Potency))
colnames(cell_pct) <- c("Sample", "Potency", "Num")
cell_pct <- spread(cell_pct, "Potency", "Num")
cell_pct[, 2:ncol(cell_pct)] <- 100 * cell_pct[, 2:ncol(cell_pct)] / rowSums(cell_pct[, 2:ncol(cell_pct)])
cell_pct <- gather(cell_pct, "Potency", "Percentage", -"Sample")

# Extract identifiers
cell_pct$Group <- gsub("\\d+", "", cell_pct$Sample)
cell_pct$Time <- gsub("\\D+", "", cell_pct$Sample)

cell_pct$Group <- factor(cell_pct$Group, levels = c("WT", "SCNT"))
cell_pct$Potency <- factor(cell_pct$Potency, levels = c("Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent", "Totipotent"))
cell_pct$Sample <- factor(cell_pct$Sample, levels = c("SCNT185", "WT185", "SCNT145", "WT145", "SCNT125", "WT125"))

p3 <- ggbarplot(cell_pct, x = "Sample", y = "Percentage", fill = "Potency", color = "Potency",
                palette = c("Differentiated" = "#5E4FA2", "Unipotent" = "#66C2A5", "Oligopotent" = "#E6F598", 
                            "Multipotent" = "#FEE08B", "Pluripotent" = "#F46D43", "Totipotent" = "#9E0142"),
                label = FALSE) +
  ylab("Percentage (%)") +
  theme(aspect.ratio = 0.7) +
  coord_flip()
ggsave("Figures/Tropho_CytoTRACE2_potency_sample_barplot.pdf", plot = p3, width = 6, height = 5)


# 02.Potency Proportion per Celltype
celltype_pct <- data.frame(table(cytotrace2_result_sce$treatment, cytotrace2_result_sce$Tropho_celltype, cytotrace2_result_sce$CytoTRACE2_Potency))
colnames(celltype_pct) <- c("Treatment", "Celltype", "Potency", "Num") 
celltype_pct <- spread(celltype_pct, "Potency", "Num")
celltype_pct[, 3:ncol(celltype_pct)] <- 100 * celltype_pct[, 3:ncol(celltype_pct)] / rowSums(celltype_pct[, 3:ncol(celltype_pct)])
celltype_pct <- gather(celltype_pct, "Potency", "Percentage", -c("Treatment", "Celltype"))

celltype_pct$Treatment <- factor(celltype_pct$Treatment, levels = c("WT", "SCNT"))
celltype_pct$Celltype <- factor(celltype_pct$Celltype, levels = c("LaTP", "LaTP2", "JZP1", "SynTII Precursor", "SynTI Precursor", "JZP2", "S-TGC Precursor", "SpT Precursor", "SynTII", "SynTI", "S-TGC", "SpT", "GlyT"))
celltype_pct$Potency <- factor(celltype_pct$Potency, levels = c("Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent", "Totipotent"))

p4 <- ggbarplot(celltype_pct, x = "Celltype", y = "Percentage", fill = "Potency", color = "Potency",
                palette = c("Differentiated" = "#5E4FA2", "Unipotent" = "#66C2A5", "Oligopotent" = "#E6F598", 
                            "Multipotent" = "#FEE08B", "Pluripotent" = "#F46D43", "Totipotent" = "#9E0142"),
                label = FALSE, facet.by = "Treatment") +
  ylab("Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.8)
ggsave("Figures/Tropho_CytoTRACE2_potency_celltype_barplot.pdf", plot = p4, width = 9, height = 6)


# 03.Potency Score per Celltype (Boxplot)
meta <- cytotrace2_result_sce@meta.data
meta$Tropho_celltype <- factor(meta$Tropho_celltype,
                               levels = c("LaTP", "LaTP2", "JZP1", "SynTII Precursor", "SynTI Precursor", "JZP2", "S-TGC Precursor", "SpT Precursor", "SynTII", "SynTI", "S-TGC", "SpT", "GlyT"))

p5 <- ggplot(meta, aes(x = Tropho_celltype, y = CytoTRACE2_Score, fill = treatment)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 0.5, outlier.color = 'black', position = position_dodge(width = 0.9)) +
  theme_bw() +
  ggpubr::stat_compare_means(aes(group = treatment), method = "wilcox.test", paired = FALSE,
                             method.args = list(alternative = "two.sided"), label = "p.signif", label.y = 0.65) +
  scale_fill_manual(values = c("WT" = "#A6CEE3", "SCNT" = "#FF9FCF")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.6)
ggsave("Figures/Tropho_CytoTRACE2_Score_boxplot.pdf", plot = p5, width = 7, height = 6)


# Joint Analysis: Potency & Cell Cycle -----------------------------------------
tropho_obj$CytoTRACE2_Potency <- as.character(cytotrace2_result_sce@meta.data[rownames(tropho_obj@meta.data), 'CytoTRACE2_Potency'])
tropho_obj$CytoTRACE2_Score <- cytotrace2_result_sce@meta.data[rownames(tropho_obj@meta.data), 'CytoTRACE2_Score']
tropho_obj$New_Phase <- factor(tropho_obj$New_Phase, levels = c('G0', 'G1', 'S', 'G2M'))

cellcycle_colors <- c("G0" = "#9ecae1", "G1" = "#deebf7", "S" = "#fdbe85", "G2M" = "#fee8c8")

# 01.VlnPlot: Potency score by cell cycle phase
p6 <- ggplot(tropho_obj@meta.data, aes(x = New_Phase, y = CytoTRACE2_Score, fill = New_Phase)) + 
  geom_violin(position = position_dodge(width = 0.9), trim = TRUE, scale = "width", width = 0.8) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) + 
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.15, size = 0.3, position = position_dodge(width = 0.9)) +        
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.15, size = 0.3, position = position_dodge(width = 0.9)) +       
  ggpubr::stat_compare_means(
    comparisons = list(c("G1", "S"), c("G1", "G2M"), c("G0", "S"), c("G0", "G2M")),
    method = "wilcox.test", label = "p.signif", paired = FALSE, 
    method.args = list(alternative = "two.sided")
  ) +
  scale_fill_manual(values = cellcycle_colors) +
  theme_bw() +
  labs(x = "Cell Cycle Phase", y = 'Potency Score') +
  theme(aspect.ratio = 1, legend.position = "none") 
ggsave("Figures/Tropho_cellcycle_potency_score_boxplot.pdf", plot = p6, width = 5, height = 4)


# 02.Barplot (Potency by Phase)
pct <- data.frame(table(tropho_obj$CytoTRACE2_Potency, tropho_obj$New_Phase))
colnames(pct) <- c("Potency", "Phase", "Num")
pct <- spread(pct, "Potency", "Num")
pct[, 2:ncol(pct)] <- 100 * pct[, 2:ncol(pct)] / rowSums(pct[, 2:ncol(pct)])
pct <- gather(pct, "Potency", "Percentage", -c("Phase"))

pct$Potency <- factor(pct$Potency, levels = c("Differentiated", 'Unipotent', 'Oligopotent', 'Multipotent'))    
pct$Phase <- factor(pct$Phase, levels = c("G0", 'G1', 'S', 'G2M'))

p7 <- ggbarplot(pct, x = "Phase", y = "Percentage", fill = "Potency", color = "Potency",
                palette = c("Differentiated" = "#5E4FA2", "Unipotent" = "#66C2A5", "Oligopotent" = "#E6F598", "Multipotent" = "#FEE08B"),
                label = FALSE) +
  ylab("Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 1.2)
ggsave("Figures/Tropho_cellcycle_potency_barplot.pdf", plot = p7, width = 5, height = 4)
