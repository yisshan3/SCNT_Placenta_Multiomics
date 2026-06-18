# ==============================================================================
# Description: Preprocessing, WNN integration, annotation, and statistics for 10X snRNA + snATAC Multiome data.
# ==============================================================================


# Setup & Load Packages --------------------------------------------------------
set.seed(133)
setwd("/home4/ssyi/Mouse_Placenta")
rm(list = ls())

library(Signac)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(harmony)
library(DoubletFinder)
library(patchwork)
library(ggrastr)
library(SingleR)
library(scater)
library(ggpubr)


# Create Raw Seurat Object -----------------------------------------------------
create_seurat <- function(dir_path) {
  f1 <- paste(dir_path, "filtered_feature_bc_matrix.h5", sep = "/")
  f2 <- paste(dir_path, "atac_fragments.tsv.gz", sep = "/")
  
  inputdata.10x <- Read10X_h5(f1)
  rna_counts <- inputdata.10x$`Gene Expression`
  atac_counts <- inputdata.10x$Peaks
  
  seurat.obj <- CreateSeuratObject(counts = rna_counts)
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")
  
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "mm10"
  
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = f2,
    min.cells = 10,
    annotation = annotations
  )
  
  seurat.obj[["ATAC"]] <- chrom_assay
  return(seurat.obj)
}

samples <- c('WT_batch1_125', 'WT_batch2_145', 'SCNT_batch3_125', 
             'WT_batch4_185', 'SCNT_batch5_145', 'SCNT_batch6_185')
outdir  <- "/home4/ssyi/Mouse_Placenta/codes/0.PreProcess/raw_seurat_objs"
if(!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }

# Save raw objects
lapply(samples, function(sample) {
  seurat.obj <- create_seurat(sample)
  saveRDS(seurat.obj, file = paste0(outdir, "/", sample, "-raw_seurat_obj.rds.gz"), compress = 'gzip')
})


# Quality Control --------------------------------------------------------------
qc_subset <- function(seurat.obj, max_ncount_ATAC = 1e5, min_ncount_ATAC = 2e3, 
                      max_ncount_RNA = 25000, min_ncount_RNA = 1000, 
                      min_feature_RNA = 200, mt.cutoff = 5) {
  obj_sub <- subset(
    x = seurat.obj,
    subset = nCount_ATAC < max_ncount_ATAC &
      nCount_ATAC > min_ncount_ATAC &
      nCount_RNA < max_ncount_RNA &
      nCount_RNA > min_ncount_RNA &
      nFeature_RNA > min_feature_RNA &
      percent.mt < mt.cutoff
  )
  return(obj_sub)
}

names(samples) <- c("WT125", "WT145", "SCNT125", "WT185", "SCNT145", "SCNT185")

filtered_seu_obj_list <- lapply(names(samples), function(i) {
  dir_path <- '/home4/ssyi/Mouse_Placenta/codes/0.PreProcess/raw_seurat_objs/'
  sample_ident <- samples[[i]]
  seu_obj <- readRDS(paste0(dir_path, sample_ident, "-raw_seurat_obj.rds.gz"))
  
  seu_obj$orig.ident <- i
  filtered_seu_obj <- qc_subset(seu_obj)
  
  p1 <- VlnPlot(filtered_seu_obj, features = c("nFeature_RNA", "nCount_RNA", "nCount_ATAC", "percent.mt"), 
                ncol = 4, pt.size = 0, log = TRUE) + NoLegend()
  
  filtered_seu_obj <- NucleosomeSignal(object = filtered_seu_obj, assay = "ATAC")
  filtered_seu_obj$nucleosome_group <- ifelse(filtered_seu_obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  p2 <- FragmentHistogram(filtered_seu_obj, group.by = 'nucleosome_group')
  
  filtered_seu_obj <- TSSEnrichment(object = filtered_seu_obj, fast = FALSE)
  filtered_seu_obj$high.tss <- ifelse(filtered_seu_obj$TSS.enrichment > 2, 'High', 'Low')
  
  p3 <- TSSPlot(filtered_seu_obj, group.by = 'high.tss') + NoLegend()
  p4 <- DensityScatter(filtered_seu_obj, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
  
  combined_plots <- p1 + p2 + p3 + p4
  ggsave(paste0(dir_path, i, "_QC_plots.pdf"), plot = combined_plots, width = 16, height = 10)
  
  filtered_seu_obj <- subset(x = filtered_seu_obj, subset = nucleosome_signal < 2 & TSS.enrichment > 1)
  return(filtered_seu_obj)
})

names(filtered_seu_obj_list) <- names(samples)

# Merge filtered objects
merged_filtered_seu_obj <- merge(
  x = filtered_seu_obj_list[[1]],
  y = filtered_seu_obj_list[[2:length(filtered_seu_obj_list)]],
  add.cell.ids = c("WT125", "WT145", "SCNT125", "WT185", "SCNT145", "SCNT185"),
  project = "Placenta"
)
merged_filtered_seu_obj$Index <- rownames(merged_filtered_seu_obj@meta.data)
merged_filtered_seu_obj@meta.data$treatment <- gsub("\\d+", "", merged_filtered_seu_obj@meta.data$orig.ident)

saveRDS(filtered_seu_obj_list, "codes/0.PreProcess/filtered_seu_obj_list.rds.gz", compress = "gzip")
saveRDS(merged_filtered_seu_obj, "codes/0.PreProcess/merged_filtered_seu_obj.rds.gz", compress = "gzip")


# Remove Doublets --------------------------------------------------------------
seuratlist <- readRDS("codes/0.PreProcess/filtered_seu_obj_list.rds.gz")
dr         <- c(0.038, 0.051, 0.058, 0.039, 0.039, 0.041)

for (i in 1:length(seuratlist)) {
  seuratlist[[i]]$Index <- paste(seuratlist[[i]]$orig.ident, rownames(seuratlist[[i]]@meta.data), sep = "_")
  seuratlist[[i]] <- NormalizeData(seuratlist[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  seuratlist[[i]] <- FindVariableFeatures(seuratlist[[i]], selection.method = "vst", nfeatures = 2000)
  seuratlist[[i]] <- ScaleData(seuratlist[[i]], features = rownames(seuratlist[[i]]))
  seuratlist[[i]] <- RunPCA(seuratlist[[i]], npcs = 50, features = VariableFeatures(seuratlist[[i]]))
  seuratlist[[i]] <- FindNeighbors(seuratlist[[i]], reduction = "pca", dims = 1:30)
  seuratlist[[i]] <- FindClusters(seuratlist[[i]], resolution = 0.3)
  seuratlist[[i]] <- RunUMAP(seuratlist[[i]], dims = 1:30)
  seuratlist[[i]] <- RunTSNE(seuratlist[[i]], dims = 1:30)
  
  # DoubletFinder
  sweep.res.list <- paramSweep(seuratlist[[i]], PCs = 1:30, sct = FALSE) 
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  DoubletRate <- dr[i]    
  homotypic.prop <- modelHomotypic(seuratlist[[i]]$seurat_clusters)   
  nExp_poi <- round(DoubletRate * ncol(seuratlist[[i]])) 
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  seuratlist[[i]] <- doubletFinder(seuratlist[[i]], PCs = 1:30, pN = 0.25, pK = pK_bcmvn, 
                                   nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(seuratlist[[i]]@meta.data)[ncol(seuratlist[[i]]@meta.data)] <- "DoubletFinder"
}

Doublet_Index <- unlist(lapply(seuratlist, function(x) x@meta.data$Index[x@meta.data$DoubletFinder == "Doublet"]))

saveRDS(seuratlist, "codes/0.PreProcess/DoubletFinder_seu_obj_list.rds.gz", compress = 'gzip')
saveRDS(Doublet_Index, "codes/0.PreProcess/Doublet_Index.rds.gz", compress = 'gzip')

# Remove doublets
merged_filtered_seu_obj <- subset(merged_filtered_seu_obj, cells = setdiff(colnames(merged_filtered_seu_obj), Doublet_Index))
saveRDS(merged_filtered_seu_obj, "codes/0.PreProcess/merged_filtered_seu_obj.rds.gz", compress = "gzip")


# WNN Analysis -----------------------------------------------------------------
# RNA processing
DefaultAssay(merged_filtered_seu_obj) <- "RNA"
merged_filtered_seu_obj <- NormalizeData(merged_filtered_seu_obj) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE)

ElbowPlot(merged_filtered_seu_obj, ndims = 50, reduction = "pca")

# ATAC processing
DefaultAssay(merged_filtered_seu_obj) <- "ATAC"
merged_filtered_seu_obj <- RunTFIDF(merged_filtered_seu_obj)   
merged_filtered_seu_obj <- FindTopFeatures(merged_filtered_seu_obj, min.cutoff = 'q0')
merged_filtered_seu_obj <- RunSVD(merged_filtered_seu_obj)  
ElbowPlot(merged_filtered_seu_obj, ndims = 50, reduction = "lsi")

# Harmony RNA
DefaultAssay(merged_filtered_seu_obj) <- "RNA"
integrated_obj <- RunHarmony(
  object = merged_filtered_seu_obj,
  group.by.vars = 'orig.ident', reduction = 'pca', assay.use = 'RNA',
  project.dim = FALSE, reduction.save = "harmony_rna", plot_convergence = TRUE
)

# Harmony ATAC
DefaultAssay(integrated_obj) <- "ATAC"
integrated_obj <- RunHarmony(
  object = integrated_obj,
  group.by.vars = 'orig.ident', reduction = 'lsi', assay.use = 'ATAC',
  project.dim = FALSE, reduction.save = "harmony_atac", plot_convergence = TRUE
)

# WNN Integration
integrated_obj <- FindMultiModalNeighbors(
  object = integrated_obj,
  reduction.list = list("harmony_rna", "harmony_atac"),
  dims.list = list(1:50, 2:50), verbose = TRUE
)   

DefaultAssay(integrated_obj) <- "RNA"
integrated_obj <- RunUMAP(integrated_obj, reduction = 'harmony_rna', dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

DefaultAssay(integrated_obj) <- "ATAC"
integrated_obj <- RunUMAP(integrated_obj, reduction = 'harmony_atac', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

integrated_obj <- RunUMAP(integrated_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", verbose = TRUE)
integrated_obj <- FindClusters(integrated_obj, graph.name = 'wsnn', algorithm = 3, resolution = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1), verbose = FALSE)

dir.create("codes/1.Integration/", recursive = TRUE, showWarnings = FALSE)
saveRDS(integrated_obj, "codes/1.Integration/integrated_obj.rds.gz", compress = 'gzip')


# Add Allele Info (PWK Ratio) --------------------------------------------------
all_allele_info <- data.frame()

load_allele_info <- function(main_path, sample_ident) {
  file_path   <- str_c(main_path, "/outs/SNP/9.calculate_allele_ratio/ATAC_SNPsplit_summary.tab")
  atac_allele <- as.data.frame(read_tsv(file_path))
  atac_allele$pwk_ratio <- atac_allele$genome2_reads / (atac_allele$genome1_reads + atac_allele$genome2_reads)
  rownames(atac_allele) <- paste0(sample_ident, "_", atac_allele$CB, "-1")
  return(atac_allele)
}

for(i in names(samples)) {
  main_path   <- paste0("/home4/ssyi/Mouse_Placenta/2.cellranger-arc/", samples[[i]])
  atac_allele <- load_allele_info(main_path, i)
  all_allele_info <- rbind(all_allele_info, atac_allele)
}

all_allele_info$CB_ident <- rownames(all_allele_info)
all_allele_info <- all_allele_info[all_allele_info$CB_ident %in% rownames(integrated_obj@meta.data), ]
integrated_obj$pwk_ratio <- all_allele_info[rownames(integrated_obj@meta.data), "pwk_ratio"]


# Major Types Annotation -------------------------------------------------------
# 01.SingleR Reference-based Annotation
load("Public_data/elife/AllStages_AllNuclei_obj.Rdata")     
elife_count <- mouse.combined@assays$RNA@counts
elife_cluster <- data.frame(Index = rownames(mouse.combined@meta.data), cluster = mouse.combined@active.ident)   
elife_anno_cluster <- read_tsv("Public_data/elife/elife_anno_cluster.tab")

elife_anno <- merge(elife_cluster, elife_anno_cluster, by = "cluster")
elife_anno_v2 <- column_to_rownames(elife_anno, var = "Index")
elife_anno_v2$cluster <- NULL
colnames(elife_anno_v2) <- "ref_label"
elife_anno_v2 <- elife_anno_v2[colnames(elife_count), , drop = FALSE]     

elife_SE <- SummarizedExperiment(assays = list(counts = elife_count), colData = elife_anno_v2)    
elife_SE <- scater::logNormCounts(elife_SE) 

my_count <- integrated_obj@assays$RNA@counts
my_SE <- SummarizedExperiment(assays = list(counts = my_count))   
my_SE <- scater::logNormCounts(my_SE)

common_gene <- intersect(rownames(my_SE), rownames(elife_SE)) 
elife_SE <- elife_SE[common_gene, ]
my_SE <- my_SE[common_gene, ]

singleR_res <- SingleR(test = my_SE, ref = elife_SE, labels = elife_SE$ref_label)   
anno_df <- data.frame(Index = rownames(singleR_res), ref_label_from_elifeAll = singleR_res$labels)

integrated_obj@meta.data$Index <- rownames(integrated_obj@meta.data)
integrated_obj@meta.data <- integrated_obj@meta.data %>% inner_join(anno_df, by = "Index")     
rownames(integrated_obj@meta.data) <- integrated_obj@meta.data$Index

# Best cluster annotation mapping
stat_celltype_per_cluster <- t(as.data.frame.array(table(integrated_obj$ref_label_from_elifeAll, integrated_obj$wsnn_res.0.3)))     
max_celltype_per_cluster <- data.frame(anno = apply(stat_celltype_per_cluster, 1, function(row) names(row)[which.max(row)]))
max_celltype_per_cluster <- rownames_to_column(max_celltype_per_cluster, var = "cluster")

meta_cluster <- data.frame(Index = rownames(integrated_obj@meta.data), cluster = integrated_obj$wsnn_res.0.3)
meta_cluster_celltype <- merge(meta_cluster, max_celltype_per_cluster, by = "cluster")
rownames(meta_cluster_celltype) <- meta_cluster_celltype$Index
meta_cluster_celltype <- meta_cluster_celltype[rownames(integrated_obj@meta.data), ]
integrated_obj$ref_label_from_elifeAll_cluster <- meta_cluster_celltype$anno


# 02.Manually Define Major Types
DefaultAssay(integrated_obj) <- "RNA"
Idents(integrated_obj) <- integrated_obj$wsnn_res.0.3
DimPlot(integrated_obj, reduction = "wnn.umap", group.by = 'wsnn_res.0.3', label = TRUE, pt.size = 0.5) + NoLegend()
table(integrated_obj$wsnn_res.0.3, integrated_obj$ref_label_from_elifeAll_cluster)
cluster_marker <- FindAllMarkers(integrated_obj, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

new.cluster.ids <- c("0" = "Trophoblast", "1" = "Trophoblast", "2" = "Trophoblast", "3" = "Trophoblast", "4" = "Trophoblast", 
                     "5" = "Fetal_mesenchyme", "6" = "Trophoblast", "7" = "Endothelial", '8' = "Trophoblast", '9' = "Trophoblast", 
                     "10" = "Trophoblast", "11" = "Pericyte", "12" = "Decidual_stroma", "13" = "Trophoblast", "14" = "Blood_cells", 
                     "15" = "Trophoblast", "16" = "Trophoblast", "17" = "Trophoblast", '18' = "Trophoblast", '19' = "Trophoblast", 
                     "20" = "Pericyte", "21" = "Trophoblast", "22" = "Blood_cells", "23" = "Blood_cells")

integrated_obj <- RenameIdents(integrated_obj, new.cluster.ids)
integrated_obj$Major_celltype <- integrated_obj@active.ident

cell_levels <- c("Trophoblast", "Fetal_mesenchyme", "Pericyte", "Endothelial", "Decidual_stroma", "Blood_cells")
Idents(integrated_obj) <- factor(integrated_obj$Major_celltype, levels = cell_levels)

dir.create("codes/2.Celltype_anno/", recursive = TRUE, showWarnings = FALSE)
saveRDS(integrated_obj, "codes/2.Celltype_anno/integrated_obj.rds.gz", compress = 'gzip')


# 03.Annotation Plotting
cell_colors <- c("Trophoblast" = "#66c2a5", "Blood_cells" = "#fc8d62", "Decidual_stroma" = "#8da0cb", 
                 "Fetal_mesenchyme" = "#e78ac3", "Endothelial" = "#a6d854", "Pericyte" = "#ffd92f")  
sample_colors <- c("WT125" = "#d1e5f0", "WT145" = "#67a9cf", "WT185" = "#2166ac", 
                   "SCNT125" = "#fddbc7", "SCNT145" = "#ef8a62", "SCNT185" = "#b2182b")

p1 <- DotPlot(integrated_obj, features = c(
  "Pparg", "Perp", "Lepr", "Ctsq",                               # Trophoblast
  "Gata4", "Kit", "Pdpn",                                        # Fetal_mesenchyme
  "Acta2", "Pdgfrb", "Col1a1",                                   # Pericyte
  "Pecam1", "Kdr", "Tek",                                        # Endothelial
  "Pgr", "Pbx1",                                                 # Decidual_stroma
  "Hbb-y", "Hba-x", "Ptprc", "Adgre1", "Bank1", "Btla", "Gzmc"   # Blood_cells
)) +
  scale_colour_gradient2(low = "#1f78b4", mid = "white", high = "firebrick3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("Figures/Global_marker_Dotplot.pdf", p1, width = 11, height = 3.8)

p2 <- FeaturePlot(integrated_obj, features = c("Pparg", "Gata4", "Col1a1", "Pecam1", "Pgr", "Ptprc"), 
                  cols = c("lightgrey", "firebrick3"), reduction = "wnn.umap", combine = TRUE, ncol = 3) +
  theme(aspect.ratio = 1)
p2 <- rasterize(p2, dpi = 300)
ggsave("Figures/Global_marker_Featureplot.pdf", p2, width = 11, height = 6)

p3 <- DimPlot(integrated_obj, reduction = "wnn.umap", label = TRUE, cols = cell_colors) + theme(aspect.ratio = 1)
p4 <- DimPlot(integrated_obj, reduction = "umap.rna", label = TRUE, cols = cell_colors) + theme(aspect.ratio = 1)
p5 <- DimPlot(integrated_obj, reduction = "umap.atac", label = TRUE, cols = cell_colors) + theme(aspect.ratio = 1)
p6 <- DimPlot(integrated_obj, reduction = "wnn.umap", group.by = 'orig.ident', label = FALSE, cols = sample_colors) + theme(aspect.ratio = 1)
ggsave("Figures/Global_wnnUMAP.pdf", p3, width = 6.5, height = 5)
ggsave("Figures/Global_rnaUMAP.pdf", p4, width = 6.5, height = 5)
ggsave("Figures/Global_atacUMAP.pdf", p5, width = 6.5, height = 5)
ggsave("Figures/Global_wnnUMAP_samples.pdf", p6, width = 6.5, height = 5)

p7 <- FeaturePlot(integrated_obj, features = c("pwk_ratio"), cols = c("#74c476", "#9e9ac8"), reduction = "wnn.umap", 
                  combine = TRUE, min.cutoff = 0, max.cutoff = 0.5) + theme(aspect.ratio = 1)
ggsave("Figures/Global_PWK_ratio_UMAP.pdf", p7, width = 6.5, height = 4)


# 04.FindAllMarkers
celltype_marker <- FindAllMarkers(integrated_obj, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_tsv(celltype_marker, "codes/2.Celltype_anno/MajorTypes_marker.txt")


# Modality Correlation Validation ----------------------------------------------
seu_RNA <- readRDS("codes/2.Celltype_anno/integrated_obj.rds.gz")
DefaultAssay(seu_RNA) <- "RNA"
seu_RNA <- NormalizeData(seu_RNA) %>% FindVariableFeatures() %>% ScaleData()

seu_ATAC <- readRDS("codes/2.Celltype_anno/integrated_obj.rds.gz")
DefaultAssay(seu_ATAC) <- "ATAC"

gene.activities <- GeneActivity(seu_ATAC, features = VariableFeatures(seu_RNA, assay = "RNA"))
seu_ATAC[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)  

DefaultAssay(seu_ATAC) <- "ACTIVITY"
seu_ATAC <- NormalizeData(seu_ATAC) %>% ScaleData(features = rownames(seu_ATAC))

# Cross-modality validation: Predict ATAC labels from RNA anchors to evaluate concordance
transfer.anchors <- FindTransferAnchors(reference = seu_RNA, query = seu_ATAC, features = VariableFeatures(object = seu_RNA, assay = "RNA"), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = seu_RNA$Major_celltype, 
                                     weight.reduction = seu_ATAC[["lsi"]], dims = 2:30)
seu_ATAC <- AddMetaData(seu_ATAC, metadata = celltype.predictions)

seu_ATAC$annotation_correct <- seu_ATAC$predicted.id == seu_ATAC$Major_celltype
seu_ATAC$Major_celltype     <- as.character(seu_ATAC$Major_celltype)

predictions <- table(seu_ATAC$Major_celltype, seu_ATAC$predicted.id)
predictions <- as.data.frame(predictions / rowSums(predictions))
predictions$Var1 <- factor(predictions$Var1, levels = cell_levels)
predictions$Var2 <- factor(predictions$Var2, levels = cell_levels)

p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + 
  geom_tile() + 
  scale_fill_gradient(name = "Fraction of cells", low = "#e5f5f9", high = "#2ca25f") + 
  labs(x = "Cell type annotation (RNA)", y = "Predicted cell type label (ATAC)") + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), aspect.ratio = 1.2)
ggsave("Figures/Global_RNA_ATAC_cor_heatmap.pdf", p1, width = 6, height = 5)


# Statistics & Plotting --------------------------------------------------------
# 01.Cell Origin Stat
integrated_obj$cell_origin <- ifelse(integrated_obj$pwk_ratio > 0.4, 'Fetus', 'Surrogate mother')
cell_pct <- data.frame(table(integrated_obj$Major_celltype, integrated_obj$cell_origin))
colnames(cell_pct) <- c("Celltype", "Origin", "Num") 

cell_pct <- spread(cell_pct, "Origin", "Num")
cell_pct[, 2:ncol(cell_pct)] <- 100 * cell_pct[, 2:ncol(cell_pct)] / rowSums(cell_pct[, 2:ncol(cell_pct)])
cell_pct <- gather(cell_pct, "Origin", "Percentage", -"Celltype")

cell_pct$Origin   <- factor(cell_pct$Origin, levels = c('Fetus', 'Surrogate mother'))
cell_pct$Celltype <- factor(cell_pct$Celltype, levels = cell_levels)

p1 <- ggbarplot(cell_pct, x = "Celltype", y = "Percentage", fill = "Origin", color = "Origin", 
                palette = c('#af8dc3', '#7fbf7b'), label = FALSE) +
  ylab("Percentage (%)") + RotatedAxis() + theme(aspect.ratio = 1.3) 
ggsave("Figures/Global_PWK_stat_cell_origin.pdf", p1, width = 6, height = 6)


# 02. Number of cell Stat
celltype_num <- as.data.frame(table(integrated_obj$Major_celltype))
colnames(celltype_num) <- c("Celltype", "Num")
celltype_num$Num_k    <- celltype_num$Num / 1000
celltype_num$Celltype <- factor(celltype_num$Celltype, levels = cell_levels)

p2 <- ggbarplot(celltype_num, x = "Celltype", y = "Num_k", fill = "Celltype", color = "Celltype", 
                palette = cell_colors, label = TRUE) +
  ylab("Number of cells (k)") + theme(aspect.ratio = 0.8) + coord_flip()
ggsave("Figures/Global_Cell_Num_byCelltype.pdf", p2, width = 8, height = 7)


# 03.Sample Percentage Stat
cell_pct <- data.frame(table(integrated_obj$orig.ident, integrated_obj$Major_celltype))
colnames(cell_pct) <- c("Sample", "Celltype", "Num") 

cell_pct <- spread(cell_pct, "Celltype", "Num")
cell_pct[, 2:ncol(cell_pct)] <- 100 * cell_pct[, 2:ncol(cell_pct)] / rowSums(cell_pct[, 2:ncol(cell_pct)])
cell_pct <- gather(cell_pct, "Celltype", "Percentage", -"Sample")

cell_pct$Sample2  <- gsub("T1", "T_1", cell_pct$Sample)
cell_pct$Group <- str_split_fixed(cell_pct$Sample2, "_", 2)[, 1]
cell_pct$Time <- str_split_fixed(cell_pct$Sample2, "_", 2)[, 2]
cell_pct$Group <- factor(cell_pct$Group, levels = c("WT", "SCNT"))
cell_pct$Sample <- factor(cell_pct$Sample, levels = c("WT125", "WT145", "WT185", "SCNT125", "SCNT145", "SCNT185"))
cell_pct$Celltype <- factor(cell_pct$Celltype, levels = c("Blood_cells", "Decidual_stroma", "Endothelial", "Pericyte", "Fetal_mesenchyme", "Trophoblast"))

p3 <- ggbarplot(cell_pct, x = "Time", y = "Percentage", fill = "Celltype", color = "Celltype", 
                palette = cell_colors, label = FALSE, facet.by = "Group") +
  ylab("Percentage (%)") + theme(aspect.ratio = 1.7) 
ggsave("Figures/Global_stat_celltype.pdf", p3, width = 7, height = 6)

cell_pct$Percentage <- round(cell_pct$Percentage, digits = 1)
Endo_cell_pct <- cell_pct[cell_pct$Celltype == 'Endothelial', ]
Peri_cell_pct <- cell_pct[cell_pct$Celltype == 'Pericyte', ]

p4 <- ggplot(Endo_cell_pct, aes(x = Time, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Percentage), color = "white", size = 4, hjust = 0.5, vjust = 1.5, position = position_dodge(width = 0.9)) +
  theme_bw() + scale_fill_manual(values = c(WT = "#A6CEE3", SCNT = "#FF9FCF")) +
  labs(y = "Percentage (%)", title = 'Endothelial') + theme(aspect.ratio = 1.3)
p5 <- ggplot(Peri_cell_pct, aes(x = Time, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Percentage), color = "white", size = 4, hjust = 0.5, vjust = 1.5, position = position_dodge(width = 0.9)) +
  theme_bw() + scale_fill_manual(values = c(WT = "#A6CEE3", SCNT = "#FF9FCF")) +
  labs(y = "Percentage (%)", title = 'Pericyte') + theme(aspect.ratio = 1.3)
ggsave("Figures/EC_pct.pdf", p4, width = 6, height = 5)
ggsave("Figures/Peri_pct.pdf", p5, width = 6, height = 5)


# 04.Celltype by Sample/Treatment Stat
sample_pct <- data.frame(table(integrated_obj$orig.ident, integrated_obj$Major_celltype))
colnames(sample_pct) <- c("Sample", "Celltype", "Num") 

sample_pct <- spread(sample_pct, "Sample", "Num") %>% column_to_rownames(var = "Celltype")
sample_pct["col_sum", ] <- colSums(sample_pct)
max_sum <- max(sample_pct["col_sum", ])

scale_sample_pct <- apply(sample_pct, 2, function(col) { col * max_sum / col[nrow(sample_pct)] }) %>% as.data.frame()
scale_sample_pct <- scale_sample_pct[1:(nrow(sample_pct) - 1), ] %>% rownames_to_column(var = "Celltype")

scale_sample_pct[, 2:7] <- 100 * scale_sample_pct[, 2:7] / rowSums(scale_sample_pct[, 2:7])
scale_sample_pct <- gather(scale_sample_pct, "Sample", "Percentage", -"Celltype")

scale_sample_pct$Sample2 <- gsub("T1", "T_1", scale_sample_pct$Sample)
scale_sample_pct$Group   <- str_split_fixed(scale_sample_pct$Sample2, "_", 2)[, 1]
scale_sample_pct$Time    <- str_split_fixed(scale_sample_pct$Sample2, "_", 2)[, 2]

scale_sample_pct$Group  <- factor(scale_sample_pct$Group, levels = c("WT", "SCNT"))
scale_sample_pct$Sample <- factor(scale_sample_pct$Sample, levels = c("SCNT185", "SCNT145", "SCNT125", "WT185", "WT145", "WT125"))

p6 <- ggbarplot(scale_sample_pct, x = "Celltype", y = "Percentage", fill = "Sample", color = "Sample", 
                palette = sample_colors, label = FALSE) +
  ylab("Percentage (%)") + theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.8) + coord_flip()
ggsave("Figures/Global_stat_sample.pdf", p6, width = 6, height = 7)

treatment_colors <- c("WT" = "#A6CEE3", "SCNT" = "#FF9FCF")
p7 <- ggbarplot(scale_sample_pct, x = "Celltype", y = "Percentage", fill = "Group", color = "Group", 
                palette = treatment_colors, label = FALSE) +
  ylab("Percentage (%)") + theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.8) + coord_flip()
ggsave("Figures/Global_stat_treatment.pdf", p7, width = 6, height = 7)
