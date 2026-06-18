# ==============================================================================
# Description: Cell cycle scoring, AUCell signatures (Hypoxia, Glycolysis, DNA Damage, Lactation), GSVA, GO Enrichment, and Imprinted Genes expression.
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
library(cowplot)
library(ggrastr)
library(clusterProfiler)
library(AUCell)
library(ggridges) 
library(pheatmap)
library(org.Mm.eg.db)


# Load Data --------------------------------------------------------------------
tropho_obj <- readRDS("codes/3.Subtype_anno/tropho_obj.rds.gz")


# Cell Cycle Scoring -----------------------------------------------------------
DefaultAssay(tropho_obj) <- "RNA"
s.genes <- CaseMatch(search = cc.genes$s.genes, match = rownames(tropho_obj))   
g2m.genes <- CaseMatch(search = cc.genes$g2m.genes, match = rownames(tropho_obj))
tropho_obj <- CellCycleScoring(tropho_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

g0_genes <- c("Cdkn1a", "Cdkn1b", "Cdkn2a", "Cdkn2b")  
g1_genes <- c("Mki67", "Ccnd1", "Ccnd2", "Ccnd3", "Cdk4", "Cdk6") 
tropho_obj <- AddModuleScore(tropho_obj, features = list(g0_genes), name = "g0_score")
tropho_obj <- AddModuleScore(tropho_obj, features = list(g1_genes), name = "g1_score")

tropho_obj$New_Phase <- tropho_obj$Phase
old_g1_cells <- WhichCells(tropho_obj, expression = Phase == "G1")
g1_cells <- tropho_obj@meta.data[tropho_obj$g1_score1 > tropho_obj$g0_score1 & tropho_obj$g1_score1 > 0 & tropho_obj$Phase == "G1", "Index"]
g0_cells <- setdiff(old_g1_cells, g1_cells)
tropho_obj@meta.data[g0_cells, "New_Phase"] <- "G0"
tropho_obj$New_Phase <- factor(tropho_obj$New_Phase, levels = c("G0", "G1", "S", "G2M"))

cellcycle_colors <- c("G0" = "#A6CEE3", "G1" = "#1F78B4", "S" = "#B2DF8A", "G2M" = "#33A02C")

# 01.FeaturePlot
p1 <- DimPlot(tropho_obj, group.by = "New_Phase", reduction = "wnn.umap", cols = cellcycle_colors, split.by = "treatment")
ggsave("Figures/Tropho_cellcycle_treatment_UMAP.pdf", p1, width = 5 * 2, height = 5)


# 02.Barplot for each sample
cellcycle_pct <- data.frame(table(tropho_obj$orig.ident, tropho_obj$New_Phase))
colnames(cellcycle_pct) <- c("Sample", "Phase", "Num")
cellcycle_pct <- spread(cellcycle_pct, "Phase", "Num")
cellcycle_pct[, 2:ncol(cellcycle_pct)] <- 100 * cellcycle_pct[, 2:ncol(cellcycle_pct)] / rowSums(cellcycle_pct[, 2:ncol(cellcycle_pct)])
cellcycle_pct <- gather(cellcycle_pct, "Phase", "Percentage", -c("Sample"))

cellcycle_pct$Sample <- factor(cellcycle_pct$Sample, levels = c("SCNT185", "WT185", "SCNT145", "WT145", "SCNT125", "WT125"))   
cellcycle_pct$Phase <- factor(cellcycle_pct$Phase, levels = c("G2M", "S", "G1", "G0"))

p2 <- ggbarplot(cellcycle_pct, x = "Sample", y = "Percentage", fill = "Phase", color = "Phase",
                palette = cellcycle_colors, label = FALSE) +
  ylab("Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.7) +
  coord_flip()
ggsave("Figures/Tropho_cellcycle_stat_sample.pdf", p2, width = 5, height = 4)


# Gene Set Analysis ------------------------------------------------------------
Interested_geneSets <- read.gmt("codes/3.Subtype_anno/Interested_term.gmt")
cells_rankings <- AUCell_buildRankings(tropho_obj@assays$RNA@data)
geneSets <- lapply(unique(Interested_geneSets$term), function(x) { Interested_geneSets$gene[Interested_geneSets$term == x] })
names(geneSets) <- unique(Interested_geneSets$term)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.1)

Precursor_Diff_cell_levels <- c("LaTP", "LaTP2", "SynTII Precursor", "SynTI Precursor", "S-TGC Precursor", "JZP1", "JZP2", "SpT Precursor", "SynTII", "SynTI", "S-TGC", "GlyT", "SpT")   

# 01.Hypoxia
geneSet <- "GOBP_CELLULAR_RESPONSE_TO_HYPOXIA"   
tropho_obj$AUC.Hypoxia <- as.numeric(getAUC(cells_AUC)[geneSet, ])

df <- data.frame(tropho_obj@meta.data, tropho_obj@reductions$wnn.umap@cell.embeddings)
p1 <- ggplot(df, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = AUC.Hypoxia)) + 
  geom_point(size = 0.3) + 
  scale_color_gradientn(colours = c("#4575b4", "#e0f3f8", "#d73027"), values = scales::rescale(c(0, 0.35, 1))) +   
  theme_light(base_size = 10) +
  labs(title = "GOBP_CELLULAR_RESPONSE_TO_HYPOXIA") +
  theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 0.95) +
  facet_wrap(~ treatment)
ggsave("Figures/Lab_Hypoxia_score_umap.pdf", p1, width = 8, height = 5)

AUC_score <- tropho_obj@meta.data[, c("Index", "orig.ident", "treatment", "AUC.Hypoxia", "time", "Tropho_celltype")]
AUC_score$Pre_diff_label <- ifelse(AUC_score$Tropho_celltype %in% c("SynTII", "SynTI", "S-TGC", "GlyT", "SpT"), "Differentiated", "Precursor")
AUC_score$Pre_diff_label <- factor(AUC_score$Pre_diff_label, levels = c('Precursor', 'Differentiated'))
AUC_score$Tropho_celltype <- factor(AUC_score$Tropho_celltype, levels = Precursor_Diff_cell_levels)

p2 <- ggplot(AUC_score, aes(x = Tropho_celltype, y = AUC.Hypoxia, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9)) +     
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.15, size = 0.3, position = position_dodge(width = 0.9)) +        
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.15, size = 0.3, position = position_dodge(width = 0.9)) +      
  ggpubr::stat_compare_means(aes(group = treatment), method = "wilcox.test", paired = FALSE,
                             method.args = list(alternative = "two.sided"), label = "p.signif", label.y = 0.23) +
  scale_fill_manual(values = c("#A6CEE3", "#FF9FCF")) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text.x = element_text(face = "plain", angle = 45, hjust = 1, vjust = 1, size = 12)) +
  ylim(c(0.04, 0.24)) +
  labs(x = "", y = "Hypoxia score", title = "GOBP CELLULAR RESPONSE TO HYPOXIA") +
  facet_grid(~ Pre_diff_label, scales = 'free', space = 'free')
ggsave("Figures/Lab_Hypoxia_score_boxplot.pdf", p2, width = 8, height = 5)


# 02.Glycolysis
geneSet <- "REACTOME_GLYCOLYSIS"   
tropho_obj$AUC.Glycolysis <- as.numeric(getAUC(cells_AUC)[geneSet, ])

AUC_score <- tropho_obj@meta.data[, c("Index", "orig.ident", "treatment", "AUC.Glycolysis", "time", "Tropho_celltype")]
AUC_score$Pre_diff_label <- ifelse(AUC_score$Tropho_celltype %in% c("SynTII", "SynTI", "S-TGC", "GlyT", "SpT"), "Differentiated", "Precursor")
AUC_score$Pre_diff_label <- factor(AUC_score$Pre_diff_label, levels = c('Precursor', 'Differentiated'))
AUC_score$Tropho_celltype <- factor(AUC_score$Tropho_celltype, levels = Precursor_Diff_cell_levels)

p3 <- ggplot(AUC_score, aes(x = Tropho_celltype, y = AUC.Glycolysis, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9)) + 
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.15, size = 0.3, position = position_dodge(width = 0.9)) + 
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.15, size = 0.3, position = position_dodge(width = 0.9)) + 
  ggpubr::stat_compare_means(aes(group = treatment), method = "wilcox.test", paired = FALSE,
                             method.args = list(alternative = "two.sided"), label = "p.signif", label.y = 0.235) +
  scale_fill_manual(values = c("#A6CEE3", "#FF9FCF")) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text.x = element_text(face = "plain", angle = 45, hjust = 1, vjust = 1, size = 12)) +
  ylim(c(0.04, 0.24)) +
  labs(x = "", y = "Glycolysis score", title = "REACTOME GLYCOLYSIS") +
  facet_grid(~ Pre_diff_label, scales = 'free', space = 'free')
ggsave("Figures/Lab_Glycolysis_score_boxplot.pdf", p3, width = 8, height = 5)


# 03.DNA damage
geneSet <- "GOBP_CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS"
tropho_obj$AUC.DNA_Damage <- as.numeric(getAUC(cells_AUC)[geneSet, ])

df <- data.frame(tropho_obj@meta.data, tropho_obj@reductions$wnn.umap@cell.embeddings)
p4 <- ggplot(df, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = AUC.DNA_Damage)) + 
  geom_point(size = 0.3) + 
  scale_color_gradientn(colours = c("#4575b4", "#e0f3f8", "#d73027")) +
  theme_light(base_size = 10) +
  labs(title = "GOBP_CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS") +
  theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 0.9) +
  facet_wrap(~ treatment)
ggsave("Figures/DNA_damage_UMAP.pdf", p4, width = 8, height = 5)

AUC_score <- tropho_obj@meta.data[, c("Index", "orig.ident", "treatment", "AUC.DNA_Damage", "time", "Tropho_celltype")]
AUC_score$Pre_diff_label <- ifelse(AUC_score$Tropho_celltype %in% c("SynTII", "SynTI", "S-TGC", "GlyT", "SpT"), "Differentiated", "Precursor")
AUC_score$Pre_diff_label <- factor(AUC_score$Pre_diff_label, levels = c('Precursor', 'Differentiated'))
AUC_score$Tropho_celltype <- factor(AUC_score$Tropho_celltype, levels = Precursor_Diff_cell_levels)

p5 <- ggplot(AUC_score, aes(x = Tropho_celltype, y = AUC.DNA_Damage, fill = treatment)) + 
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9)) +  
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.15, size = 0.3, position = position_dodge(width = 0.9)) +      
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.15, size = 0.3, position = position_dodge(width = 0.9)) +      
  ggpubr::stat_compare_means(aes(group = treatment), method = "wilcox.test", paired = FALSE,
                             method.args = list(alternative = "two.sided"), label = "p.signif", label.y = 0.225) +
  scale_fill_manual(values = c("#A6CEE3", "#FF9FCF")) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.3), axis.text.x = element_text(face = "plain", angle = 45, hjust = 1, vjust = 1, size = 12)) +
  labs(x = "", y = "AUC score", title = "GOBP_CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS") +
  facet_grid(~ Pre_diff_label, scales = 'free', space = 'free')
ggsave("Figures/DNA_damage_boxplot.pdf", p5, width = 8.5, height = 5)


# 04.Lactation
geneSet <- "GOBP_LACTATION"
tropho_obj$AUC.Lactation <- as.numeric(getAUC(cells_AUC)[geneSet, ])
JZ_obj <- subset(tropho_obj, subset = Tropho_celltype %in% c("SpT", "SpT Precursor", "GlyT", "JZP2", "JZP1"))

AUC_score <- JZ_obj@meta.data[, c("Index", "orig.ident", "treatment", "AUC.Lactation", "time", "Tropho_celltype")]
ylim1 <- boxplot.stats(AUC_score$AUC.Lactation)$stats[c(1, 5)] 
AUC_score$Tropho_celltype <- factor(AUC_score$Tropho_celltype, levels = c("SpT", "SpT Precursor", "GlyT", "JZP2", "JZP1"))

p6 <- ggplot(AUC_score, aes(x = Tropho_celltype, y = AUC.Lactation, color = treatment)) + 
  geom_boxplot(position = position_dodge(width = 0.9)) +     
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.15, size = 0.3, position = position_dodge(width = 0.9)) +        
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.15, size = 0.3, position = position_dodge(width = 0.9)) +        
  coord_cartesian(ylim = ylim1 * 1.5) +
  ggpubr::stat_compare_means(aes(group = treatment), method = "wilcox.test", paired = FALSE,
                             method.args = list(alternative = "two.sided"), label = "p.signif", label.y = 0.37) +
  scale_color_manual(values = c("#A6CEE3", "#FF9FCF")) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text.x = element_text(face = "plain", angle = 45, hjust = 1, vjust = 1, size = 12), aspect.ratio = 0.9) +
  labs(x = "", y = "AUC score", title = "Lactation Score")
ggsave("Figures/JZ_lactation_score_celltype_boxplot.pdf", p6, width = 6, height = 5)


# GSVA (Glycogen Metabolism) ---------------------------------------------------
tropho_obj$treatment_celltype <- paste(tropho_obj$treatment, tropho_obj$Tropho_celltype, sep = "_")

GSEA_glycogen <- read.gmt("codes/3.Subtype_anno/GSEA_glycogen.gmt")
glycogen_geneSets <- lapply(unique(GSEA_glycogen$term), function(x) { GSEA_glycogen$gene[GSEA_glycogen$term == x] })
names(glycogen_geneSets) <- unique(GSEA_glycogen$term)

# Compute mean expression
cal_mean_exp <- function(object, genes_to_cal, column_to_cal) {
  exp_data <- object@assays$RNA@data
  meta <- object@meta.data[, c("Index", column_to_cal)]
  colnames(meta) <- c("Index", "Types")
  
  gene_matrix <- exp_data[rownames(exp_data) %in% genes_to_cal, rownames(meta)] %>% as.matrix()
  gene_matrix <- expm1(gene_matrix)   
  
  ave_gene_matrix <- data.frame(row.names = rownames(gene_matrix))
  for (i in unique(meta$Types)) {
    columns_to_average <- meta[meta$Types == i, "Index"]
    ave_gene_matrix[i] <- rowMeans(gene_matrix[, columns_to_average, drop = FALSE])
  }
  return(ave_gene_matrix)
}
celltype_aveExp <- cal_mean_exp(tropho_obj, rownames(tropho_obj), "treatment_celltype")
celltype_aveExp <- as.matrix(celltype_aveExp)
GSEA_glycogen_gsva_scores <- gsva(celltype_aveExp, glycogen_geneSets, method = "gsva", parallel.sz = 5)

levels <- paste0(rep(c("WT_", "SCNT_"), time = 13), 
                 rep(c("LaTP", "SynTII Precursor", "SynTII", "LaTP2", "SynTI Precursor", "SynTI", "S-TGC Precursor", "S-TGC", "JZP1", "JZP2", "GlyT", "SpT Precursor", "SpT"), each = 2))
GSEA_glycogen_gsva_scores <- GSEA_glycogen_gsva_scores[, levels]

df <- as.data.frame(GSEA_glycogen_gsva_scores['REACTOME_GLYCOGEN_SYNTHESIS', , drop = FALSE])
df1 <- data.frame(WT = t(df)[grepl('WT', colnames(df)), ],
                  SCNT = t(df)[grepl('SCNT', colnames(df)), ],
                  row.names = gsub('WT_', '', colnames(df)[grepl('WT', colnames(df))]))

df1 <- arrange(df1, desc(SCNT))

pdf(file = "Figures/REACTOME_GLYCOGEN_SYNTHESIS_GSVA_heatmap.pdf", width = 3, height = 5)
pheatmap(df1, 
         scale = "none", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         display_numbers = TRUE, 
         main = 'REACTOME_GLYCOGEN_SYNTHESIS')
dev.off()


# GO Enrichment ----------------------------------------------------------------
# 01.SCNT Precursor DEGs
SCNT_Precursor_DEGs_counts <- read_tsv('codes/3.Subtype_anno/Trophoblast/Gene_List/merge_SCNT_Precursor_markers_sort_uniq_count.txt', col_names = FALSE)
colnames(SCNT_Precursor_DEGs_counts) <- c('Gene.ID', 'count')
ego_bp <- enrichGO(
  gene = SCNT_Precursor_DEGs_counts$Gene.ID,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
go_data <- as_tibble(ego_bp@result)
go_data$GeneRatio <- as.numeric(str_split_i(go_data$GeneRatio, '/', 1)) / as.numeric(str_split_i(go_data$GeneRatio, '/', 2))
go_data$BgRatio <- as.numeric(str_split_i(go_data$BgRatio, '/', 1)) / as.numeric(str_split_i(go_data$BgRatio, '/', 2))
go_data$fold_enrichment <- go_data$GeneRatio / go_data$BgRatio
go_data$log10_p <- -log10((go_data$pvalue))

Enriched_term <- c('mitotic cell cycle phase transition', 'chromosome segregation', 'regulation of cell cycle phase transition',
                   'double-strand break repair', 'regulation of mitotic cell cycle phase transition', 'cell cycle checkpoint signaling',
                   'recombinational repair', 'signal transduction in response to DNA damage', 'DNA damage checkpoint signaling',
                   'mitotic G2 DNA damage checkpoint signaling')
go_data <- go_data[go_data$Description %in% Enriched_term, ]

p1 <- ggplot(data = go_data, aes(x = reorder(Description, log10_p), y = log10_p, fill = fold_enrichment)) +
  geom_bar(stat = 'identity', width = 0.9) + theme_classic() +
  scale_fill_gradient(name = "Fold enrichment", low = "#deebf7", high = "#2171b5") +
  geom_text(aes(label = Description), y = 0, hjust = 0, color = "black", size = 5) +
  theme(axis.text.x = element_text(colour = "black", size = 15), axis.text.y = element_text(colour = "black", size = 15),
        axis.title.x = element_text(colour = "black", size = 15), axis.title.y = element_text(colour = "black", size = 15),
        legend.text = element_text(colour = "black", size = 13), legend.title = element_text(colour = "black", size = 15),
        aspect.ratio = 0.7) +
  coord_flip()
ggsave("Figures/SCNT_Precursor_DEGs_clusterprofiler_GOBP_barplot.pdf", p1, width = 12, height = 5)


# 02.SCNT DEGs Overlap
SCNT_DEGs_count <- read_tsv('codes/3.Subtype_anno/Trophoblast/Gene_List/merge_SCNT_subtype_upDEGs_sort_uniq_count.txt', col_names = FALSE)
colnames(SCNT_DEGs_count) <- c('Gene.ID', 'count')
filter_SCNT_DEGs_count <- SCNT_DEGs_count[SCNT_DEGs_count$count > 6, ]  

ego_bp <- enrichGO(gene = filter_SCNT_DEGs_count$Gene.ID, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL',
                   ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

go_data <- as_tibble(ego_bp@result)
go_data$GeneRatio <- as.numeric(str_split_i(go_data$GeneRatio, '/', 1)) / as.numeric(str_split_i(go_data$GeneRatio, '/', 2))
go_data$BgRatio <- as.numeric(str_split_i(go_data$BgRatio, '/', 1)) / as.numeric(str_split_i(go_data$BgRatio, '/', 2))
go_data$fold_enrichment <- go_data$GeneRatio / go_data$BgRatio
go_data$log10_p <- -log10((go_data$pvalue))

Enriched_term <- c('mitotic cell cycle phase transition', "double-strand break repair", 'double-strand break repair via homologous recombination',
                   'G1/S transition of mitotic cell cycle', 'positive regulation of DNA repair', 'DNA damage checkpoint signaling',
                   'signal transduction in response to DNA damage', 'mitotic DNA damage checkpoint signaling',
                   'mitotic G2 DNA damage checkpoint signaling', 'mitotic intra-S DNA damage checkpoint signaling')
go_data <- go_data[go_data$Description %in% Enriched_term, ]

p2 <- ggplot(data = go_data, aes(x = reorder(Description, log10_p), y = log10_p, fill = fold_enrichment)) +
  geom_bar(stat = 'identity', width = 0.9) + theme_classic() +
  scale_fill_gradient(name = "Fold enrichment", low = "#fee5d9", high = "#cb181d") +
  geom_text(aes(label = Description), y = 0, hjust = 0, color = "black", size = 5) +
  theme(axis.text.x = element_text(colour = "black", size = 15), axis.text.y = element_text(colour = "black", size = 15),
        axis.title.x = element_text(colour = "black", size = 15), axis.title.y = element_text(colour = "black", size = 15),
        legend.text = element_text(colour = "black", size = 13), legend.title = element_text(colour = "black", size = 15),
        aspect.ratio = 0.7) +
  coord_flip()
ggsave("Figures/SCNT_cutoff7_upDEGs_clusterprofiler_GOBP_barplot.pdf", p2, width = 13, height = 5)


# Imprinted Genes (Ridge Plots) ------------------------------------------------
celltype_allele_fpkm <- read_tsv("codes/9.sinto_SNPsplit/RNA/Celltype/Stringtie/merged_fpkm.tab")
Precursor_celltypes <- c("LaTP", "LaTP2", "SynTII Precursor", "SynTI Precursor", "S-TGC Precursor", "JZP1", "JZP2", "SpT Precursor")

plot_density_ridges <- function(gene) {
  plot_fpkm <- celltype_allele_fpkm[celltype_allele_fpkm$Gene.ID == gene, ]
  plot_fpkm <- gather(plot_fpkm, key = "Group", value = "FPKM", -"Gene.ID")
  
  plot_fpkm$Genome <- str_split_fixed(plot_fpkm$Group, '\\.', 2)[, 2]
  plot_fpkm$Label <- str_split_fixed(plot_fpkm$Group, '\\.', 2)[, 1]
  plot_fpkm$Treatment <- str_split_fixed(plot_fpkm$Group, '_', 2)[, 1]
  plot_fpkm$Celltype <- str_split_fixed(plot_fpkm$Label, '_', 2)[, 2]
  
  plot_fpkm <- plot_fpkm[plot_fpkm$Celltype %in% Precursor_celltypes, ]
  plot_fpkm$Genome <- ifelse(plot_fpkm$Genome == 'genome1', "Maternal_allele(C57)", "Paternal_allele(PWK)")
  
  plot_fpkm$Genome <- factor(plot_fpkm$Genome, levels = c("Maternal_allele(C57)", "Paternal_allele(PWK)"))
  plot_fpkm$Treatment <- factor(plot_fpkm$Treatment, levels = c("WT", "SCNT"))
  
  p <- ggplot(plot_fpkm, aes(x = FPKM, y = Genome, fill = Treatment)) +  
    geom_density_ridges(scale = .95, rel_min_height = .01, quantiles = c(0.25, 0.5, 0.75)) +
    scale_fill_manual(values = c("#A6CEE3", "#FF9FCF")) +  
    scale_discrete_manual("point_color", values = c("#A6CEE3", "#FF9FCF"), guide = "none") +
    coord_cartesian(clip = "off") +
    guides(fill = guide_legend(override.aes = list(fill = c("#A6CEE3", "#FF9FCF"), color = NA, point_color = NA))) +
    theme_bw() +
    ggtitle(paste0(gene, " density")) +
    theme(aspect.ratio = 0.4)
  return(p)
}

p8 <- plot_density_ridges("Sfmbt2")
p9 <- plot_density_ridges("Smoc1")
p10 <- plot_density_ridges("Slc38a4")
p12 <- plot_density_ridges("Jade1")
ggsave("Figures/Sfmbt2_ridgeplot.pdf", p8, width = 6, height = 4)
ggsave("Figures/Smoc1_ridgeplot.pdf", p9, width = 6, height = 4)
ggsave("Figures/Slc38a4_ridgeplot.pdf", p10, width = 6, height = 4)
ggsave("Figures/Jade1_ridgeplot.pdf", p12, width = 6, height = 4)


# VlnPlot with FDR -------------------------------------------------------------
JZ_celltype <- c("SpT", "SpT Precursor", "GlyT", "JZP2", "JZP1")
Lab_celltype <- c("LaTP", "SynTII Precursor", "SynTII", "LaTP2", "SynTI Precursor", "SynTI", "S-TGC Precursor", "S-TGC")

# load DEGs results
load_DEGs_results <- function(target_cells) {
  results_list <- list()
  for (ct in target_cells) {
    file_path <- paste0("codes/3.Subtype_anno/Trophoblast/", gsub(" ", "_", ct), "_DEGs.csv")
    if(file.exists(file_path)) {
      results_list[[ct]] <- read.csv(file_path)
    }
  }
  return(results_list)
}

Plot_VlnPlot_FDR <- function(seurat_obj = tropho_obj, target_gene, target_celltypes, 
                             FindMarkers_results_DEGs, aspect.ratio = 0.6) {
  
  sub_obj <- subset(seurat_obj, subset = Tropho_celltype %in% target_celltypes)
  expr_matrix <- sub_obj@assays$RNA@data
  meta_data <- sub_obj@meta.data
  
  results_list <- list()
  for (ct in target_celltypes) {
    cells_g1 <- rownames(meta_data)[meta_data$Tropho_celltype == ct & meta_data$treatment == "WT"]
    cells_g2 <- rownames(meta_data)[meta_data$Tropho_celltype == ct & meta_data$treatment == "SCNT"]
    
    SCNT_celltype_DEGs <- FindMarkers_results_DEGs[[ct]]
    p_val_FDR_FindMarkers <- NA_real_ 
    
    if(!is.null(SCNT_celltype_DEGs) && target_gene %in% SCNT_celltype_DEGs$gene) { 
      p_val_FDR_FindMarkers <- SCNT_celltype_DEGs[SCNT_celltype_DEGs$gene == target_gene, 'p_val_BH(FDR)']
    }
    
    results_list[[ct]] <- data.frame(
      Celltype = ct,
      Gene = target_gene,
      n_Group1 = paste0("WT(", length(cells_g1), ")"),
      n_Group2 = paste0("SCNT(", length(cells_g2), ")"),
      p_val_FDR_FindMarkers = p_val_FDR_FindMarkers
    )
  }
  final_results <- do.call(rbind, results_list)
  
  final_results$group1 <- "WT"  
  final_results$group2 <- "SCNT"
  final_results$FDR_label <- format.pval(signif(final_results$p_val_FDR_FindMarkers, digits = 2), digits = 2, eps = 2.2e-16)
  final_results$FDR_label <- paste0("p = ", final_results$FDR_label)
  final_results$FDR_label <- gsub("p = <", "p <", final_results$FDR_label)
  final_results$FDR_label <- ifelse(final_results$p_val_FDR_FindMarkers < 2.2e-16, "p < 2.2e-16", final_results$FDR_label)
  final_results$FDR_label[is.na(final_results$p_val_FDR_FindMarkers)] <- "ns"
  
  max_expression <- max(as.numeric(expr_matrix[target_gene, ]), na.rm = TRUE)
  y_pos_adaptive <- ifelse(max_expression == 0, 0.1, max_expression * 0.995)
  
  p <- VlnPlot(sub_obj, features = target_gene, group.by = "Tropho_celltype",
               split.by = "treatment", split.plot = FALSE, pt.size = 0, combine = TRUE, assay = "RNA") +
    geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
    ggpubr::stat_pvalue_manual(final_results, label = "FDR_label", y.position = y_pos_adaptive, 
                               x = "Celltype", tip.length = 0, inherit.aes = FALSE) +
    scale_fill_manual(values = c("WT" = "#A6CEE3", "SCNT" = "#FF9FCF")) + 
    theme(aspect.ratio = aspect.ratio)
  
  return(p)
}

# Precursor
precursor_genes <- c('Rad51', 'Brca1', 'Mre11a', 'Rbbp8', 'Sfmbt2', 'Smoc1', 'Slc38a4', 'Jade1')
precursor_DEGs_list <- load_DEGs_results(Precursor_celltypes)

for (gene in precursor_genes) {
  p <- Plot_VlnPlot_FDR(target_gene = gene, target_celltypes = Precursor_celltypes, 
                        FindMarkers_results_DEGs = precursor_DEGs_list, aspect.ratio = 0.6)
  ggsave(paste0("Figures/Precursor_", gene, "_vlnplot.pdf"), p, width = 8, height = 6)
}

# GlyT
GlyT <- subset(tropho_obj, Tropho_celltype == "GlyT")
GlyT$treatment <- factor(GlyT$treatment, levels = c("WT", "SCNT"))
glyt_genes <- c('Hk1', 'Pgm2', 'Ugp2', 'Gbe1', 'Slc2a1', 'Slc2a2', 'Slc2a3', 'Slc2a4')
glyt_DEGs_list <- load_DEGs_results("GlyT")

for (gene in glyt_genes) {
  p <- Plot_VlnPlot_FDR(seurat_obj = GlyT, target_gene = gene, target_celltypes = 'GlyT', 
                        FindMarkers_results_DEGs = glyt_DEGs_list, aspect.ratio = 1.6)
  ggsave(paste0("Figures/GlyT_", gene, "_vlnplot.pdf"), p, width = 3, height = 5)
}

# JZ
JZ <- subset(tropho_obj, Tropho_celltype %in% JZ_celltype)
JZ$treatment <- factor(JZ$treatment, levels = c("WT", "SCNT"))
JZ$Tropho_celltype <- factor(JZ$Tropho_celltype, levels = c("SpT", "SpT Precursor", "GlyT", "JZP2", "JZP1"))
JZ_genes <- c('Flt1', 'Tpbpa', 'Prl3b1', 'Prl8a9')
JZ_DEGs_list <- load_DEGs_results(JZ_celltype)

for (gene in JZ_genes) {
  p <- Plot_VlnPlot_FDR(seurat_obj = JZ, target_gene = gene, target_celltypes = JZ_celltype, 
                        FindMarkers_results_DEGs = JZ_DEGs_list, aspect.ratio = 0.6)
  ggsave(paste0("Figures/JZ_", gene, "_celltype_vlnplot.pdf"), p, width = 7, height = 6)
}

# Lab
Lab <- subset(tropho_obj, subset = Tropho_celltype %in% Lab_celltype)
Lab$Tropho_celltype <- factor(Lab$Tropho_celltype, levels = c("LaTP", "SynTII Precursor", "SynTII", "LaTP2", "SynTI Precursor", "SynTI", "S-TGC Precursor", "S-TGC"))
lab_DEGs_list <- load_DEGs_results(Lab_celltype)

p_lab_vegfa <- Plot_VlnPlot_FDR(seurat_obj = Lab, target_gene = 'Vegfa', target_celltypes = Lab_celltype, 
                                FindMarkers_results_DEGs = lab_DEGs_list, aspect.ratio = 0.5)
ggsave("Figures/Lab_Vegfa_celltype_vlnplot.pdf", p_lab_vegfa, width = 7, height = 6)

