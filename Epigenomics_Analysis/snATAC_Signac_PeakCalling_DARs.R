# ==============================================================================
# Description: Peak calling, Differential Accessibility (DARs), integration with imprinted gene lists, and metadata export for Sinto.
# ==============================================================================


# Setup & Load Packages --------------------------------------------------------
set.seed(133)
setwd("/home4/ssyi/Mouse_Placenta")
rm(list = ls())

library(chromVAR)
library(cicero)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(patchwork)
library(tidyverse)
library(ggVennDiagram)
library(future)

plan("multicore", workers = 10)
options(future.globals.maxSize = 50 * 1024^3)


# Load Data --------------------------------------------------------------------
tropho_obj <- readRDS("codes/3.Subtype_anno/tropho_obj.rds.gz")
tropho_obj$treatment_celltype <- paste0(tropho_obj$treatment, "_", tropho_obj$Tropho_celltype)


# Peak Calling -----------------------------------------------------------------
peaks <- CallPeaks(tropho_obj, group.by = "treatment_celltype")

# Quality control on peaks
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

# Quantify counts in peaks
macs2_counts <- FeatureMatrix(
  fragments = Fragments(tropho_obj),
  features = peaks,
  cells = colnames(tropho_obj)
)

tropho_obj[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(tropho_obj),
  annotation = GetGRangesFromEnsDb(EnsDb.Mmusculus.v79)
)

dir.create("codes/8.Signac/", recursive = TRUE, showWarnings = FALSE)
saveRDS(tropho_obj, "codes/8.Signac/tropho_obj_signac.rds.gz", compress = "gzip")


# Differential Accessibility Analysis (DARs) -----------------------------------
Idents(tropho_obj) <- tropho_obj$treatment_celltype
cell_levels <- c("LaTP", "SynTII Precursor", "SynTII", "LaTP2", "SynTI Precursor", 
                 "SynTI", "S-TGC Precursor", "S-TGC", "JZP1", "JZP2", "GlyT", "SpT Precursor", "SpT")

SCNT_pst_celltype_DARs.list <- list()
WT_pst_celltype_DARs.list <- list()

# Identify DARs (SCNT vs WT)
for (celltype in cell_levels) {
  SCNT_pst_celltype_DARs.list[[celltype]] <- FindMarkers(
    object = tropho_obj,
    group.by = "treatment_celltype",
    ident.1 = paste0("SCNT_", celltype),
    ident.2 = paste0("WT_", celltype),
    only.pos = TRUE,
    test.use = 'LR',
    latent.vars = 'nCount_peaks'
  )
  
  WT_pst_celltype_DARs.list[[celltype]] <- FindMarkers(
    object = tropho_obj,
    group.by = "treatment_celltype",
    ident.1 = paste0("WT_", celltype),
    ident.2 = paste0("SCNT_", celltype),
    only.pos = TRUE,
    test.use = 'LR',
    latent.vars = 'nCount_peaks'
  )
}
save(SCNT_pst_celltype_DARs.list, WT_pst_celltype_DARs.list, file = "codes/8.Signac/DARs_list.RData")

# filter DARs
precursor_celltypes <- c("LaTP", "LaTP2", "SynTII Precursor", "SynTI Precursor", 
                         "S-TGC Precursor", "JZP1", "JZP2", "SpT Precursor")

filter_SCNT_precursor_DARs <- lapply(SCNT_pst_celltype_DARs.list[precursor_celltypes], function(x) {
  x %>% dplyr::filter(p_val < 0.05 & abs(avg_log2FC) > 0.25 & (pct.1 > 0.1 | pct.2 > 0.1))
})
names(filter_SCNT_precursor_DARs) <- precursor_celltypes

outdir  <- "./codes/12.Modifications/DARs/"
if(!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }

lapply(names(filter_SCNT_precursor_DARs), function(x){
  celltype <- gsub(" ",'_',x)
  DARs <- filter_SCNT_precursor_DARs[[x]]
  bed <- str_split_fixed(rownames(DARs), "-", 3) %>% as.data.frame()
  write_tsv(bed, paste0("codes/12.Modifications/DARs/SCNT_",celltype,'_DARs.bed'), col_names = F)
})

# merge SCNT precursor DARs as DDARs
SCNT_precursor_DARs_peaks <- Reduce(union, lapply(filter_SCNT_precursor_DARs, rownames))
write_tsv(as.data.frame(str_split_fixed(SCNT_precursor_DARs_peaks, "-", 3)), 
          "codes/12.Modifications/DARs/DDARs.bed", col_names = FALSE)

# Generate Random background regions
DefaultAssay(tropho_obj) <- "peaks"
NC_Non_DARs <- sample(setdiff(rownames(tropho_obj), SCNT_precursor_DARs_peaks), size = length(SCNT_precursor_DARs_peaks))
write_tsv(as.data.frame(str_split_fixed(NC_Non_DARs, "-", 3)), 
          "codes/12.Modifications/DARs/Random_regions.bed", col_names = FALSE)


# Imprinting & Venn Diagrams ---------------------------------------------------
close_genes <- ClosestFeature(tropho_obj, regions = SCNT_precursor_DARs_peaks)
up_DEGs <- read_tsv("codes/3.Subtype_anno/Trophoblast/Gene_List/merge_SCNT_subtype_upDEGs_sort_uniq_count.txt", col_names = FALSE)$X1

# Define Imprinting Gene Sets
imprinting_raw <- read_tsv("codes/10.Imprinting/mouse_imprinting_genes.txt")
imprinting_v2 <- imprinting_raw %>% filter(`Expressed Allele` %in% c("Paternal", "Maternal"))
k27me3_genes <- c('Sfmbt2', 'Gab1', 'Smoc1', 'Slc38a4', 'Jade1', 'Xist')
dna_methyl_imprinted <- setdiff(imprinting_v2$Gene, k27me3_genes)

# Venn Visualization
gene_list1 <- list(DDARs = unique(close_genes$gene_name), Up_DEGs = up_DEGs, K27me3_Imprinted = k27me3_genes)
p1 <- ggVennDiagram(gene_list1) + scale_fill_gradient(low = "white", high = "#67a9cf")

gene_list2 <- list(DDARs = unique(close_genes$gene_name), Up_DEGs = up_DEGs, DNAm_Imprinted = dna_methyl_imprinted)
p2 <- ggVennDiagram(gene_list2) + scale_fill_gradient(low = "white", high = "#c994c7")

ggsave("Figures/DARs_DEGs_H3K27me3_imprinting_VennDiagram.pdf", p1, width = 5, height = 5)
ggsave("Figures/DARs_DEGs_DNAmtheyl_imprinting_VennDiagram.pdf", p2, width = 5, height = 5)


# Export Metadata for Sinto ----------------------------------------------------
sinto_base <- "codes/9.sinto/ATAC/Batch"
meta <- tropho_obj@meta.data
meta$CB <- str_split_i(meta$Index, "_", 2)
meta$Tropho_celltype <- gsub(" ", "_", meta$Tropho_celltype)

samples_map <- list(
  "SCNT125" = "SCNT_batch3_125_PWK", "SCNT145" = "SCNT_batch5_145_PWK", 
  "SCNT185" = "SCNT_batch6_185_PWK", "WT125" = "WT_batch1_125_PWK", 
  "WT145" = "WT_batch2_145_PWK", "WT185" = "WT_batch4_185_PWK"
)

for (sample_id in names(samples_map)) {
  out_path <- file.path(sinto_base, samples_map[[sample_id]])
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  
  meta_subset <- meta[meta$orig.ident == sample_id, c("CB", "Tropho_celltype")]
  write_tsv(meta_subset, file.path(out_path, "Cell_meta.txt"), col_names = FALSE)
}

# Export cell types
celltype_dir <- "codes/9.sinto/ATAC/Celltype"
dir.create(celltype_dir, recursive = TRUE, showWarnings = FALSE)

tropho_celltypes <- unique(meta$Tropho_celltype)
writeLines(tropho_celltypes, file.path(celltype_dir, "Tropho_celltypes.txt"))

precursor_celltypes <- c("LaTP", "LaTP2", "SynTII_Precursor", "SynTI_Precursor", 
                         "S-TGC_Precursor", "JZP1", "JZP2", "SpT_Precursor")
writeLines(precursor_celltypes, file.path(celltype_dir, "Precursor_celltypes.txt"))

