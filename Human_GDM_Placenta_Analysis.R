# ==============================================================================
# Description: Preprocessing, annotation, sub-clustering, potency evaluation, and differential analysis of human GDM placenta scRNA-seq data.
# ==============================================================================


# Setup & Load Packages --------------------------------------------------------
set.seed(133)
setwd("/home4/ssyi/Mouse_Placenta/Public_data/GSE173193_GDM/")
rm(list = ls())

library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(harmony)
library(ComplexHeatmap)
library(circlize)
library(CytoTRACE2)
library(scProportionTest)
library(clusterProfiler)
library(org.Hs.eg.db)

options(Seurat.object.assay.version = "v3")


# Load Data --------------------------------------------------------------------
samples <- list.files("/home4/ssyi/Mouse_Placenta/Public_data/GSE173193_GDM/Data/", pattern = "GSM")

create_seurat <- function(sample, dir_path = '/home4/ssyi/Mouse_Placenta/Public_data/GSE173193_GDM/Data/') {
  f1 <- paste(dir_path, sample, "filtered_feature_bc_matrix", sep = "/")
  inputdata.10x <- Read10X(f1)
  seurat.obj <- CreateSeuratObject(counts = inputdata.10x, project = sample, min.cells = 3, min.features = 200)
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")
  return(seurat.obj)
}

for(i in samples) {
  assign(i, create_seurat(i))
}

outdir <- "./ProcessedData/0.PreProcess/raw_seurat_objs"
if(!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

seurat_list <- lapply(samples, function(sample) {
  seurat.obj <- create_seurat(sample)
  saveRDS(seurat.obj, file = paste0(outdir, '/', sample, "-raw_seurat_obj.rds"))
  return(seurat.obj)
})
names(seurat_list) <- samples

merge_seu <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = samples
)

merge_seu@meta.data <- merge_seu@meta.data %>% 
  mutate(
    group1 = case_when(
      orig.ident %in% c('GSM5261695', 'GSM5261696') ~ 'control',
      orig.ident %in% c('GSM5261697', 'GSM5261698') ~ 'gestational diabetes group',
      orig.ident %in% c('GSM5261699', 'GSM5261700') ~ 'preeclampsia',
      orig.ident %in% c('GSM5261701', 'GSM5261702') ~ 'advanced age group'
    ),
    group2 = case_when(
      orig.ident %in% c('GSM5261695', 'GSM5261696') ~ 'Ctrl',
      orig.ident %in% c('GSM5261697', 'GSM5261698') ~ 'GDM',
      orig.ident %in% c('GSM5261699', 'GSM5261700') ~ 'PE',
      orig.ident %in% c('GSM5261701', 'GSM5261702') ~ 'GL'
    ),
    age = case_when(
      orig.ident %in% c('GSM5261695') ~ '29',
      orig.ident %in% c('GSM5261696', 'GSM5261697') ~ '28',
      orig.ident %in% c('GSM5261698') ~ '33',
      orig.ident %in% c('GSM5261699') ~ '34',
      orig.ident %in% c('GSM5261699') ~ '34',
      orig.ident %in% c('GSM5261701') ~ '39',
      orig.ident %in% c('GSM5261702') ~ '38'
    )
  )

merge_seu$tissue <- 'Placenta'
merge_seu <- subset(merge_seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)   
saveRDS(merge_seu, 'ProcessedData/0.PreProcess/merge_8samples.rds.gz', compress = 'gzip')


# GDM scRNA-seq Pipeline -------------------------------------------------------
GDM_obj <- subset(merge_seu, orig.ident %in% c('GSM5261695', 'GSM5261696', 'GSM5261697', 'GSM5261698'))
GDM_obj <- NormalizeData(GDM_obj, normalization.method = "LogNormalize", scale.factor = 10000)
GDM_obj <- FindVariableFeatures(GDM_obj, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(GDM_obj), 10)
all.genes <- rownames(GDM_obj)

GDM_obj <- ScaleData(GDM_obj, features = all.genes)
GDM_obj <- RunPCA(GDM_obj)
ElbowPlot(GDM_obj, ndims = 50, reduction = "pca")

GDM_obj <- FindNeighbors(object = GDM_obj, dims = 1:30)
GDM_obj <- FindClusters(object = GDM_obj)
GDM_obj <- RunUMAP(object = GDM_obj, dims = 1:30)

GDM_obj@meta.data[GDM_obj@meta.data$orig.ident == 'GSM5261695', 'sample'] <- 'C1'
GDM_obj@meta.data[GDM_obj@meta.data$orig.ident == 'GSM5261696', 'sample'] <- 'C2'
GDM_obj@meta.data[GDM_obj@meta.data$orig.ident == 'GSM5261697', 'sample'] <- 'G1'
GDM_obj@meta.data[GDM_obj@meta.data$orig.ident == 'GSM5261698', 'sample'] <- 'G2'

DimPlot(object = GDM_obj, reduction = "umap", group.by = 'sample')
saveRDS(GDM_obj, 'ProcessedData/0.PreProcess/GDM_obj.rds.gz', compress = 'gzip')

# Harmony remove batch effect
GDM_obj <- RunHarmony(
  object = GDM_obj,
  group.by.vars = 'orig.ident',
  reduction = 'pca',
  assay.use = 'RNA',
  project.dim = FALSE,
  reduction.save = "harmony",
  plot_convergence = TRUE
)

ElbowPlot(GDM_obj, ndims = 50, reduction = "harmony")
GDM_obj <- RunUMAP(GDM_obj, reduction = 'harmony', dims = 1:30)
DimPlot(object = GDM_obj, reduction = "umap", group.by = 'sample')
saveRDS(GDM_obj, 'ProcessedData/0.PreProcess/GDM_obj_harmony.rds.gz', compress = 'gzip')


# Annotation -------------------------------------------------------------------
GDM_obj <- FindClusters(object = GDM_obj, resolution = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1))

DimPlot(object = GDM_obj, reduction = "umap", group.by = 'RNA_snn_res.0.5', label = TRUE)
DimPlot(object = GDM_obj, reduction = "umap", group.by = 'RNA_snn_res.0.8', label = TRUE)

marker_list <- list(
  Trophoblast = c("PERP", "KRT7", "GATA3", "GATA2"),
  VCT = c('PARP1', 'MET', 'EGFR', 'CDH1', 'CCNB2', "TP63"),
  EVT = c('HLA-G', 'PAPPA2', 'MMP11', 'MMP2', 'TGFB1', 'CXCR6'),
  SCT = c('ERVFRD-1', 'CYP19A1', 'CGA', 'LGALS13', 'INSL4'),
  Immune = c('PTPRC'),
  TNK = c('CD3D', 'CD3G', 'GZMA', 'XCL2', 'CCL5', 'GZMK', 'IFNG'),
  B = c('CD79A', 'CD79B', 'CD19'), 
  Mono = c('LYZ', 'CD14', 'CD300E', 'CD244', 'HLA-DRA', 'FCN1', 'CLEC12A'),
  Macro = c('CD209', 'CD163', 'AIF1', 'CD68', 'CSF1R'),
  Gran = c('FCGR3B', 'CXCL8', 'MNDA', 'SELL'),
  Myelocyte = c('TCN1', 'CEACAM8', 'MMP8', 'DEFA4', 'CAMP', 'S100A8'),
  Erythroblast = c('HBG1', "GYPC"),
  EC = c('PECAM1', 'VWF', 'ENG', 'ACKR1', 'CA4', 'CLEC4G', 'CD34'),
  Mast = c('HDC', 'CPA3')
)

DotPlot(GDM_obj, features = marker_list, group.by = "RNA_snn_res.0.5") + 
  scale_colour_gradient2(low = "lightblue", mid = "white", high = "firebrick3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

VlnPlot(GDM_obj, features = c('nFeature_RNA', "nCount_RNA", "percent.mt"), group.by = "RNA_snn_res.0.5", pt.size = 0) +
  theme(legend.position = "right")

Idents(GDM_obj) <- GDM_obj$RNA_snn_res.0.5
GDM_obj <- RenameIdents(GDM_obj, c(
  '3' = 'VCT', '4' = 'VCT', '7' = 'VCT', '9' = 'VCT', '13' = 'VCT', '17' = 'VCT',
  '2' = 'EVT',
  '12' = 'SCT',
  '0' = 'Granulocyte', '22' = 'Granulocyte',
  '5' = 'Myelocyte', '11' = 'Myelocyte',
  '8' = 'Monocyte', '14' = 'Monocyte',
  '1' = 'Macrophage', '15' = 'Macrophage',
  '6' = 'Tcell', 
  '10' = 'NK',
  '16' = 'EC',
  '18' = 'Mast',
  '20' = 'Bcell',
  '21' = 'Erythroblast',
  '19' = 'Unknown_immune'
))

GDM_obj$MinorType <- as.character(Idents(GDM_obj))
GDM_obj <- subset(GDM_obj, RNA_snn_res.0.5 != '19')

GDM_obj$MajorType <- GDM_obj$MinorType
GDM_obj@meta.data[GDM_obj$MinorType %in% c('VCT', 'EVT', 'SCT'), 'MajorType'] <- 'Trophoblast'
GDM_obj@meta.data[GDM_obj$MinorType %in% c('Tcell', 'NK'), 'MajorType'] <- 'TNK'

outdir <- "./ProcessedData/1.Annotation/"
if(!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
saveRDS(GDM_obj, 'ProcessedData/1.Annotation/GDM_obj.rds.gz', compress = 'gzip')


# MajorTypes Landscape ---------------------------------------------------------
# 01.DimPlot
GDM_obj$MajorType <- factor(GDM_obj$MajorType, levels = c('Trophoblast', 'Granulocyte', 'Myelocyte', 'Monocyte', 'Macrophage', 'TNK', 'EC', 'Mast', 'Bcell', 'Erythroblast'))

colors <- c("Trophoblast" = '#8dd3c7', "Myelocyte" = '#ffffb3', "Granulocyte" = '#b2df8a',
            "Macrophage" = '#bebada', "Monocyte" = '#80b1d3', "Mast" = '#fb8072',
            "TNK" = '#fdb462', "Bcell" = '#bc80bd', "EC" = '#fccde5', "Erythroblast" = '#ccebc5')

p1 <- DimPlot(GDM_obj, reduction = "umap", group.by = "MajorType", 
              label = FALSE, cols = colors, split.by = "group2", ncol = 3) 
ggsave("Figures/GDM_UMAP_celltype.pdf", p1, width = 5 * 2, height = 5)


# 02.DotPlot
p2 <- DotPlot(GDM_obj, features = c(
  "PERP", "KRT7", "GATA3",
  'PTPRC', 'FCGR3B', 'CSF3R', 'CXCL8', 'SELL', 'MNDA', 'S100A8', 'S100A9',
  'TCN1', 'CEACAM8', 'MMP8', 'DEFA4', 'CAMP',
  'LYZ', 'CD14', 'CD300E', 'HLA-DRA', 'FCN1', 'CLEC12A',
  'CD68', 'CD163', 'CD209', 'AIF1', 'CSF1R',
  'CD3D', 'CD3G', 'GZMA', 'XCL2', 'CCL5', 'GZMK', 'IFNG',
  'PECAM1', 'VWF', 'ENG',
  'HDC', 'CPA3', 'TPSAB1',
  'CD79A', 'CD79B', 'CD19',
  'HBG1', "GYPC"
), group.by = "MajorType") + 
  scale_colour_gradient2(low = "lightblue", mid = "white", high = "firebrick3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Figures/GDM_marker_Dotplot.pdf", p2, width = 13.5, height = 4)


# 03.Cell Proportion per Sample
cell_pct <- data.frame(table(GDM_obj$group2, GDM_obj$MajorType))
colnames(cell_pct) <- c("Group", "Celltype", "Num") 
cell_pct <- spread(cell_pct, "Celltype", "Num")
cell_pct[, 2:ncol(cell_pct)] <- 100 * cell_pct[, 2:ncol(cell_pct)] / rowSums(cell_pct[, 2:ncol(cell_pct)])
cell_pct <- gather(cell_pct, "Celltype", "Percentage", -"Group")
cell_pct$Celltype <- factor(cell_pct$Celltype, levels = c("Trophoblast", "Myelocyte", "Granulocyte", "Macrophage", "Monocyte", "Mast", "TNK", "Bcell", "EC", "Erythroblast"))

p3 <- ggbarplot(cell_pct, x = "Group", y = "Percentage", fill = "Celltype", color = "Celltype", 
                width = 0.6, label = FALSE, palette = colors) +
  ylab("Percentage (%)") + theme(aspect.ratio = 1.6)
ggsave("Figures/GDM_stat_MajorType_ratio.pdf", p3, width = 7, height = 5)


# 04.Sample Origin per Celltype
cell_pct <- data.frame(table(GDM_obj$group2, GDM_obj$MajorType))
colnames(cell_pct) <- c("Sample", "Celltype", "Num") 
sample_pct <- spread(cell_pct, "Sample", "Num") %>% column_to_rownames(var = "Celltype")
sample_pct["col_sum", ] <- colSums(sample_pct)
max_sum <- max(sample_pct["col_sum", ])

scale_sample_pct <- apply(sample_pct, 2, function(col) { col * max_sum / col[nrow(sample_pct)] }) %>% as.data.frame()
scale_sample_pct <- scale_sample_pct[1:(nrow(sample_pct) - 1), ] %>% rownames_to_column(var = "Celltype")

sample_cols <- setdiff(colnames(scale_sample_pct), "Celltype")
scale_sample_pct[, sample_cols] <- 100 * scale_sample_pct[, sample_cols] / rowSums(scale_sample_pct[, sample_cols])
scale_sample_pct <- gather(scale_sample_pct, "Sample", "Percentage", -"Celltype")
scale_sample_pct$Sample <- factor(scale_sample_pct$Sample, levels = c('GDM', 'Ctrl'))

p4 <- ggbarplot(scale_sample_pct, x = "Celltype", y = "Percentage", fill = "Sample", color = "Sample",
                palette = c('Ctrl' = '#92c5de', 'GDM' = '#f4a582'), label = FALSE) +
  ylab("Percentage (%)") + theme(aspect.ratio = 0.8) + coord_flip()
ggsave("Figures/GDM_stat_sample_ratio.pdf", p4, width = 6, height = 7)


# Proportion Analysis ----------------------------------------------------------
# 01.scProportionTest
scu <- sc_utils(GDM_obj)
scu_res <- permutation_test(
  sc_utils_obj = scu,
  cluster_identity = "MajorType",
  sample_identity = "group2",
  sample_1 = "Ctrl",
  sample_2 = "GDM",
  n_permutations = 10000
)

df <- as.data.table(scu_res@results$permutation)
df$obs_log2FD2 <- log2((df$GDM + 1e-3) / (df$Ctrl + 1e-3))

plot_df <- df %>%
  mutate(
    neg_log10_FDR = -log10(FDR + 1e-10),
    Status = case_when(
      FDR < 0.05 & obs_log2FD2 > 0 ~ "Enriched in GDM",
      FDR < 0.05 & obs_log2FD2 < 0 ~ "Depleted in GDM",
      TRUE ~ "Not Significant"
    )
  )

cluster_order <- c('EC', 'Erythroblast', 'Monocyte', 'Granulocyte', 'Myelocyte', 'Bcell', 'Macrophage', 'TNK', 'Trophoblast', 'Mast')
plot_df$clusters <- factor(plot_df$clusters, levels = cluster_order)

p5 <- ggplot(plot_df, aes(x = obs_log2FD2, y = clusters)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = obs_log2FD2, y = clusters, yend = clusters, color = Status), linewidth = 1) +
  geom_point(aes(color = Status, size = neg_log10_FDR)) +
  scale_color_manual(values = c("Enriched in GDM" = "#fb9a99", "Depleted in GDM" = "#92c5de", "Not Significant" = "gray80")) +  
  theme_bw() +
  labs(
    x = bquote("Log"[2]~"Fold Difference (GDM vs Ctrl)"),
    y = "MajorTypes",
    title = "scProportionTest",
    size = bquote("-Log"[10]~"(FDR)"),
    color = "Significance"
  ) +
  theme(
    panel.grid.major.y = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold"), 
    strip.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    aspect.ratio = 1.5
  ) +
  scale_size_continuous(range = c(1, 4), name = bquote("|Log"[2]~"FD|"))
ggsave("Figures/GDM_MajorType_scProportionTest_lollipop.pdf", p5, width = 7, height = 6)


# 02.Ro/e (Ratio of Observed to Expected)
results <- chisq.test(GDM_obj$group2, GDM_obj$MajorType)
observed_data <- results$observed
expected_data <- results$expected
ratio_oe <- observed_data / expected_data

ord <- c('EC', 'Erythroblast', 'Monocyte', 'Granulocyte', 'Myelocyte', 'Bcell', 'Macrophage', 'TNK', 'Trophoblast', 'Mast')
missing_cells <- setdiff(ord, colnames(ratio_oe))
if(length(missing_cells) > 0) {
  pad_mat <- matrix(NA, nrow = nrow(ratio_oe), ncol = length(missing_cells))
  colnames(pad_mat) <- missing_cells
  ratio_oe <- cbind(ratio_oe, pad_mat)
}
ratio_oe <- ratio_oe[, ord, drop = FALSE]

min_val <- min(ratio_oe, na.rm = TRUE)
max_val <- max(ratio_oe, na.rm = TRUE)

pdf("Figures/GDM_MajorType_Roe.pdf", width = 6, height = 2)
ht <- Heatmap(
  ratio_oe,
  col = colorRamp2(breaks = c(min_val, 1, max_val), colors = c("#4393c3", "white", "#d6604d")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = 'left',
  row_names_gp = gpar(fontsize = 13),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = 'Ro/e'),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", ratio_oe[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)
draw(ht)
dev.off()


# Trophoblast Lineage Characterization -----------------------------------------
Tropho_obj <- subset(GDM_obj, subset = MinorType %in% c('VCT', 'EVT', 'SCT'))
saveRDS(Tropho_obj, 'ProcessedData/1.Annotation/Tropho_obj.rds.gz', compress = 'gzip')
# DimPlot
Tropho_colors <- c("VCT" = '#a6cee3', "SCT" = '#fb9a99', "EVT" = '#cab2d6')

p6 <- DimPlot(Tropho_obj, reduction = "umap", group.by = "MinorType", 
              label = TRUE, cols = Tropho_colors, split.by = "group2") 
ggsave("Figures/Hs_GDM_UMAP_Tropho.pdf", p6, width = 5 * 2, height = 5)

# DotPlot
p7 <- DotPlot(Tropho_obj, features = c(
  'PARP1', 'MET', 'CDH1', 'CCNB2', "TP63", 'EGFR',
  'ERVFRD-1', 'CYP19A1', 'CGA', 'LGALS13', 'INSL4',
  'HLA-G', 'PAPPA2', 'MMP11', 'MMP2'
), group.by = "MinorType") + 
  scale_colour_gradient2(low = "lightblue", mid = "white", high = "firebrick3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Figures/GDM_tropho_marker_Dotplot.pdf", p7, width = 6, height = 2.5)


# CytoTRACE2 -------------------------------------------------------------------
cytotrace2_result_sce <- cytotrace2(
  Tropho_obj, 
  is_seurat = TRUE, 
  slot_type = "counts", 
  species = 'human',
  seed = 133
)

# 01.Potency Proportion per Celltype
cellcycle_pct <- data.frame(table(cytotrace2_result_sce$group2, cytotrace2_result_sce$MinorType, cytotrace2_result_sce$CytoTRACE2_Potency))
colnames(cellcycle_pct) <- c("Treatment", "Celltype", "Potency", "Num") 
cellcycle_pct <- spread(cellcycle_pct, "Potency", "Num")
cellcycle_pct[, 3:ncol(cellcycle_pct)] <- 100 * cellcycle_pct[, 3:ncol(cellcycle_pct)] / rowSums(cellcycle_pct[, 3:ncol(cellcycle_pct)])
cellcycle_pct <- gather(cellcycle_pct, "Potency", "Percentage", -c("Treatment", "Celltype"))

cellcycle_pct$Treatment <- factor(cellcycle_pct$Treatment, levels = c("Ctrl", "GDM"))
cellcycle_pct$Celltype <- factor(cellcycle_pct$Celltype, levels = c("VCT", 'SCT', 'EVT'))
cellcycle_pct$Potency <- factor(cellcycle_pct$Potency, levels = c("Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent", "Totipotent"))

p8 <- ggbarplot(cellcycle_pct, x = "Celltype", y = "Percentage", fill = "Potency", color = "Potency",
                palette = c("Differentiated" = "#5E4FA2", "Unipotent" = "#66C2A5", "Oligopotent" = "#E6F598", 
                            "Multipotent" = "#FEE08B", "Pluripotent" = "#F46D43", "Totipotent" = "#9E0142"),
                label = FALSE, facet.by = "Treatment") +
  ylab("Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 1.6)
ggsave("Figures/Hs_GDM_CytoTRACE2_potency_barplot_celltype.pdf", p8, width = 5.5, height = 5)


# 02.Potency Score Comparison
meta <- cytotrace2_result_sce@meta.data
meta$MinorType <- factor(meta$MinorType, levels = c("VCT", 'SCT', 'EVT'))

p9 <- ggplot(meta, aes(x = MinorType, y = CytoTRACE2_Score, fill = group2)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 0.1, position = position_dodge(width = 0.9)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(Ctrl = "#A6CEE3", GDM = "#FF9FCF")) +
  ggpubr::stat_compare_means(aes(group = group2), method = "wilcox.test", paired = FALSE,
                             method.args = list(alternative = "two.sided"), label = "p.signif", label.y = 0.65)
ggsave("Figures/Hs_GDM_CytoTRACE2_potency_score_boxplot.pdf", p9, width = 6, height = 5)


# GO Enrichment ----------------------------------------------------------------
pro_FindMarkers <- function(Seurat_Obj, ident.1, ident.2 = NULL, fdr_cutoff = 0.05) {
  markers <- FindMarkers(Seurat_Obj, ident.1 = ident.1, ident.2 = ident.2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers$gene <- rownames(markers)
  markers$q_val_fdr <- p.adjust(markers$p_val, method = "BH")
  markers <- markers %>% dplyr::filter(q_val_fdr < fdr_cutoff)
  return(markers)
}

Idents(Tropho_obj) <- Tropho_obj$group2
GDM_tropho_marker <- pro_FindMarkers(Tropho_obj, ident.1 = "GDM", ident.2 = "Ctrl", fdr_cutoff = 0.05)

ego_bp <- enrichGO(
  gene = GDM_tropho_marker$gene,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',
  ont = "BP",  
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

go_data <- as_tibble(ego_bp@result)
go_data$GeneRatio <- as.numeric(str_split_i(go_data$GeneRatio, '/', 1)) / as.numeric(str_split_i(go_data$GeneRatio, '/', 2))
go_data$BgRatio <- as.numeric(str_split_i(go_data$BgRatio, '/', 1)) / as.numeric(str_split_i(go_data$BgRatio, '/', 2))
go_data$fold_enrichment <- go_data$GeneRatio / go_data$BgRatio
go_data$log10_p <- -log10((go_data$pvalue))

Enriched_term <- c('wound healing', 'regulation of apoptotic signaling pathway', 'response to decreased oxygen levels',
                   'response to hypoxia', 'regulation of vasculature development', 'cellular response to decreased oxygen levels',
                   'cellular response to hypoxia', 'regulation of intrinsic apoptotic signaling pathway')
go_data <- go_data[go_data$Description %in% Enriched_term, ]

p10 <- ggplot(data = go_data, aes(x = reorder(Description, log10_p), y = log10_p, fill = fold_enrichment)) +
  geom_bar(stat = 'identity', width = 0.9) +
  theme_classic() +
  scale_fill_gradient(name = "Fold enrichment", low = "#EFE8F9", high = "#8B6CB7") +
  geom_text(aes(label = Description), y = 0, hjust = 0, color = "black", size = 5) +
  theme(axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15),
        legend.text = element_text(colour = "black", size = 13),
        legend.title = element_text(colour = "black", size = 15),
        aspect.ratio = 0.7) +
  coord_flip()
ggsave("Figures/Hs_GDM_tropho_DEGs_clusterprofiler_GOBP_barplot.pdf", p10, width = 13, height = 5)

