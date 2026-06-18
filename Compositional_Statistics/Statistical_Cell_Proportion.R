# ==============================================================================
# Description: Subpopulation proportion analysis using scProportionTest, Ro/e, and data preparation/visualization for scCODA.
# ==============================================================================


# Setup & Load Packages --------------------------------------------------------
set.seed(133)
setwd("/home4/ssyi/Mouse_Placenta")
rm(list = ls())

library(scProportionTest)
library(ggrepel)
library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(circlize)


# Load Data ----------------------------------------------------------------
tropho_obj <- readRDS("codes/3.Subtype_anno/tropho_obj.rds.gz")


# scProportionTest ---------------------------------------------------------
obj_125 <- subset(tropho_obj, subset = time == "125")
obj_145 <- subset(tropho_obj, subset = time == "145")
obj_185 <- subset(tropho_obj, subset = time == "185")

# Permutation test for each timepoint
scu_125 <- sc_utils(obj_125)
scu_145 <- sc_utils(obj_145)
scu_185 <- sc_utils(obj_185)

E125_res <- permutation_test(
  sc_utils_obj = scu_125,
  cluster_identity = "Tropho_celltype",
  sample_identity = "treatment",
  sample_1 = "WT",
  sample_2 = "SCNT",
  n_permutations = 10000
)

E145_res <- permutation_test(
  sc_utils_obj = scu_145,
  cluster_identity = "Tropho_celltype",
  sample_identity = "treatment",
  sample_1 = "WT",
  sample_2 = "SCNT",
  n_permutations = 10000
)

E185_res <- permutation_test(
  sc_utils_obj = scu_185,
  cluster_identity = "Tropho_celltype",
  sample_identity = "treatment",
  sample_1 = "WT",
  sample_2 = "SCNT",
  n_permutations = 10000
)

# Merge and visualize results
extract_perm <- function(res, stage_label) {
  dt <- as.data.table(res@results$permutation)
  dt[, stage := stage_label]
  dt[, sig := FDR < 0.05]
  dt[, neglogFDR := -log10(FDR + 1e-300)]
  return(dt)
}

dt_all <- rbindlist(list(
  extract_perm(E125_res, "E12.5"),
  extract_perm(E145_res, "E14.5"),
  extract_perm(E185_res, "E18.5")
), use.names = TRUE, fill = TRUE)

dt_all[, stage := factor(stage, levels = c("E12.5", "E14.5", "E18.5"))]

# Lollipop Visualization
df <- dt_all
df$obs_log2FD2 <- log2((df$SCNT + 1e-3) / (df$WT + 1e-3))

dir.create("codes/Review/3.1.scProportionTest/", recursive = TRUE, showWarnings = FALSE)
write_tsv(df, 'codes/Review/3.1.scProportionTest/Tropho_scProportionTest_result.tsv')

plot_df <- df %>%
  filter(stage != "Meta") %>%
  mutate(
    neg_log10_FDR = -log10(FDR + 1e-10),
    Status = case_when(
      FDR < 0.05 & obs_log2FD2 > 0 ~ "Enriched in SCNT",
      FDR < 0.05 & obs_log2FD2 < 0 ~ "Depleted in SCNT",
      TRUE ~ "Not Significant"
    )
  )

plot_df$stage <- factor(plot_df$stage, levels = c("E12.5", "E14.5", "E18.5"))
cluster_order <- c("SynTII", "SpT", "S-TGC", "SynTI", "GlyT", "SynTII Precursor", "SpT Precursor",
                   "S-TGC Precursor", "JZP2", "JZP1", "SynTI Precursor", "LaTP2", "LaTP")
plot_df$clusters <- factor(plot_df$clusters, levels = cluster_order)

p1 <- ggplot(plot_df, aes(x = obs_log2FD2, y = clusters)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = obs_log2FD2, y = clusters, yend = clusters, color = Status), linewidth = 1) +
  geom_point(aes(color = Status, size = neg_log10_FDR)) +
  facet_wrap(~ stage, ncol = 3) +
  scale_color_manual(values = c("Enriched in SCNT" = "#fb9a99", "Depleted in SCNT" = "#92c5de", "Not Significant" = "gray80")) +  
  theme_bw() +
  labs(
    x = bquote("Log"[2]~"Fold Difference (SCNT vs WT)"),
    y = "Trophoblast Subtypes",
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
    aspect.ratio = 2
  ) +
  scale_size_continuous(range = c(1, 4), name = bquote("|Log"[2]~"FD|"))
ggsave("Figures/Tropho_scProportionTest_lollipop.pdf", p1, width = 10, height = 9)


# Ratio of Observed to Expected (Ro/e) -----------------------------------------
Tropho_obj_list <- list(obj_125, obj_145, obj_185)
names(Tropho_obj_list) <- c('E125', 'E145', 'E185')

dir.create("codes/Review/3.1.Roe/", recursive = TRUE, showWarnings = FALSE)
lapply(c('E125', 'E145', 'E185'), function(i) {
  obj <- Tropho_obj_list[[i]]
  # Chi-square test
  results <- chisq.test(obj$treatment, obj$Tropho_celltype)
  observed_data <- results$observed
  expected_data <- results$expected
  ratio_oe <- observed_data / expected_data
  
  ord <- c("LaTP", "JZP1", "JZP2", "LaTP2", "SynTI Precursor", "S-TGC Precursor", "SpT Precursor", "SynTII Precursor", "GlyT", "SynTI", "SynTII", "S-TGC", "SpT")
  
  missing_cells <- setdiff(ord, colnames(ratio_oe))
  if(length(missing_cells) > 0) {
    pad_mat <- matrix(NA, nrow = nrow(ratio_oe), ncol = length(missing_cells))
    colnames(pad_mat) <- missing_cells
    ratio_oe <- cbind(ratio_oe, pad_mat)
  }
  ratio_oe <- ratio_oe[, ord, drop = FALSE]
  
  # Generate Ro/e heatmaps
  min_val <- min(ratio_oe, na.rm = TRUE)
  max_val <- max(ratio_oe, na.rm = TRUE)
  
  pdf(paste0("codes/Review/3.1.Roe/",'Tropho_',i,'_Roe.pdf'),width=7, height = 2)
  ht <- Heatmap(
    ratio_oe,
    col = colorRamp2(breaks = c(min_val, 1, max_val), colors = c("#4393c3", "white", "#d6604d")),
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_side = 'left',
    row_names_gp = gpar(fontsize = 13),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(title = 'Ro/e'),
    cell_fun = function(j, k, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", ratio_oe[k, j]), x, y, gp = gpar(fontsize = 10))
    }
  )
  draw(ht)
  dev.off()
})


# scCODA -------------------------------------------------------------------
# Prepare Input metadata for scCODA
E145_4samples_tropho <- readRDS("codes/Review/E145_4samples_tropho.rds.gz")
E145_Tropho_meta <- E145_4samples_tropho@meta.data[, c("orig.ident", "time", "treatment", "Tropho_celltype")]
E145_Tropho_meta$time <- '145'
E145_Tropho_meta$cell_id <- rownames(E145_Tropho_meta)
colnames(E145_Tropho_meta) <- c("sample_id", "stage", "condition", "cell_type", "cell_id")

dir.create("codes/Review/3.1.scCODA/", recursive = TRUE, showWarnings = FALSE)
write.csv(E145_Tropho_meta, "codes/Review/3.1.scCODA/E145_Tropho_cell_metadata_for_sccoda.csv", quote = FALSE, row.names = FALSE)


# Visualization from scCODA Results
# Note: scCODA statistical model was executed via Python (scCODA_Analysis.py)
df <- read.csv("codes/Review/3.1.scCODA/scCODA_result.csv", stringsAsFactors = FALSE)
colnames(df) <- c("Covariate", "CellType", "FinalParameter", "HDI_3", "HDI_97", "SD", "Inclusion_prob", "Expected_Sample", "Log2FC")

df <- df %>%
  arrange(Inclusion_prob) %>%
  mutate(CellType = factor(CellType, levels = CellType)) 

p2 <- ggplot(df, aes(x = Inclusion_prob, y = CellType)) +
  geom_segment(aes(x = 0, xend = Inclusion_prob, y = CellType, yend = CellType), color = "grey70", size = 1.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#E64B35", size = 0.5) +
  geom_point(aes(fill = Log2FC), shape = 21, size = 5, color = "black", stroke = 0.5) +
  scale_fill_gradient2(
    low = "#4393c3", mid = "white", high = "#d6604d", midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    name = "Log2 Fold Change\n(SCNT vs WT)"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_classic() +
  labs(
    x = "Inclusion Probability",
    y = "Cell Types",
    title = "scCODA (E14.5)"
  ) +
  theme(
    axis.text.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 11),
    axis.title = element_text(face = "bold", size = 13),
    title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    aspect.ratio = 0.9
  )
ggsave('Figures/scCODA_E145_lollipop.pdf', p2, width = 7, height = 7)
