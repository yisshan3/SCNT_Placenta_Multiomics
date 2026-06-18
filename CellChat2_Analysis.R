# ==============================================================================
# Description: Infer and visualize cell-cell communication networks between labyrinthine trophoblast lineages and Endothelial/Pericyte cells.
# ==============================================================================


# Setup & Load Packages --------------------------------------------------------
set.seed(133)
setwd("/home4/ssyi/Mouse_Placenta")
rm(list = ls())

library(Seurat)
library(Signac)
library(tidyverse)
library(dplyr)
library(CellChat)
library(patchwork)
library(ggpubr)
library(ComplexHeatmap)


# Load Data --------------------------------------------------------------------
integrated_obj <- readRDS("codes/2.Celltype_anno/integrated_obj.rds.gz")
tropho_obj <- readRDS("codes/3.Subtype_anno/tropho_obj.rds.gz")

# Unify the Celltype column before merging
Lab <- subset(tropho_obj, Tropho_celltype %in% c("LaTP", "SynTII Precursor", "SynTII", "LaTP2", "SynTI Precursor", "SynTI", "S-TGC Precursor", "S-TGC"))
Lab$Celltype <- as.character(Lab$Tropho_celltype)
EC_Peri <- subset(integrated_obj, Major_celltype %in% c("Endothelial", "Pericyte"))
EC_Peri$Celltype <- as.character(EC_Peri$Major_celltype)

Lab_EC_Peri <- merge(x = Lab, y = EC_Peri)

Lab_EC_Peri_WT <- subset(Lab_EC_Peri, subset = treatment == "WT")
Lab_EC_Peri_SCNT <- subset(Lab_EC_Peri, subset = treatment == "SCNT")


# Define CellChat Function -----------------------------------------------------
Run_Cellchat <- function(Seurat_obj, sample_name) {
  cellchat <- createCellChat(
    object = Seurat_obj,
    group.by = "Celltype",
    assay = "RNA"
  )  
  
  cellchat@DB <- CellChat::CellChatDB.mouse
  
  cellchat <- subsetData(cellchat)   
  future::plan("multisession", workers = 6)   
  cellchat <- identifyOverExpressedGenes(cellchat)   
  cellchat <- identifyOverExpressedInteractions(cellchat)    
  cellchat <- projectData(cellchat, CellChat::PPI.mouse)    
  
  # Infer communication network (Ligand-Receptor level)
  cellchat <- computeCommunProb(
    cellchat,
    population.size = TRUE,  
    raw.use = FALSE
  )   
  cellchat <- filterCommunication(cellchat, min.cells = 10) 
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, paste0("codes/6.Cellchat/net_lr_", sample_name, ".csv"), row.names = FALSE)
  
  # Infer communication network (Signaling Pathway level)
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp <- subsetCommunication(cellchat, slot.name = "netP")
  write.csv(df.netp, paste0("codes/6.Cellchat/net_pathway_", sample_name, ".csv"), row.names = FALSE)
  
  cellchat <- aggregateNet(cellchat)  
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  return(cellchat)
}


# Run CellChat -----------------------------------------------------------------
cellchat_Lab_EC_Peri_WT <- Run_Cellchat(Lab_EC_Peri_WT, "WT")    
cellchat_Lab_EC_Peri_SCNT <- Run_Cellchat(Lab_EC_Peri_SCNT, "SCNT") 
dir.create("codes/6.Cellchat/", recursive = TRUE, showWarnings = FALSE)
saveRDS(cellchat_Lab_EC_Peri_WT, "codes/6.Cellchat/cellchat_Lab_EC_Peri_WT.rds.gz", compress = "gzip")
saveRDS(cellchat_Lab_EC_Peri_SCNT,"codes/6.Cellchat/cellchat_Lab_EC_Peri_SCNT.rds.gz", compress = "gzip")


# Visualization ----------------------------------------------------------------
# Update and merge CellChat objects
cellchat_WT <- updateCellChat(cellchat_Lab_EC_Peri_WT)
cellchat_SCNT <- updateCellChat(cellchat_Lab_EC_Peri_SCNT)

object.list <- list(WT = cellchat_WT, SCNT = cellchat_SCNT)
merge_cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Differential interactions (Circle Diagram)
pdf(file = "Figures/Lab_cellchat2_netVisual_diffInteraction.pdf", width = 10, height = 7)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(merge_cellchat, weight.scale = TRUE)
netVisual_diffInteraction(merge_cellchat, weight.scale = TRUE, measure = "weight")
dev.off()

# Differential interactions (Heatmap)
p1 <- netVisual_heatmap(merge_cellchat)
p2 <- netVisual_heatmap(merge_cellchat, measure = "weight")
pdf("Figures/Lab_cellchat2_diffInteraction_heatmap.pdf", width = 10, height = 6)
print(p1 + p2)
dev.off()


# Specific Signaling Interaction Pathways --------------------------------------
source_cells <- c("S-TGC", "SynTI", "SynTII")
target_cells_v1 <- c("Endothelial", "Pericyte")
target_cells_v2 <- c("Endothelial", "S-TGC", "SynTI", "SynTII")

# 01.Conserved and Specific Signaling Pathways
p3 <- rankNet(
  merge_cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE,
  sources.use = source_cells,
  targets.use = target_cells_v1,
  color.use = c("#A6CEE3", "#FF9FCF"), 
  cutoff.pvalue = 0.05
)
pdf(file = "Figures/Lab_cellchat2_signaling_pathway_STGC_SynTI_SynTII_2Endo_Peri.pdf", width = 5, height = 8.5)
print(p3)
dev.off()


# 02.Interaction Pairs (Bubble Plots)
# Angiogenesis related pathways
p4 <- netVisual_bubble(
  merge_cellchat,
  sources.use = source_cells,
  targets.use = target_cells_v1,
  signaling = c("VEGF", "NOTCH", "PDGF", 'PECAM1', 'CDH5', "TGFb"),   
  comparison = c(1, 2),
  angle.x = 45,
  remove.isolate = FALSE
)

# LAMININ pathway
p5 <- netVisual_bubble(
  merge_cellchat, 
  sources.use = source_cells, 
  targets.use = target_cells_v2,
  comparison = c(1, 2), 
  signaling = c("LAMININ"),
  angle.x = 45, 
  remove.isolate = FALSE
)
ggsave("Figures/Lab_cellchat2_STGC_SynTI_SynTII_Endo_Peri_LR.pdf", p4, width = 5.5, height = 7)
ggsave("Figures/Lab_cellchat2_STGC_SynTI_SynTII_Endo_LAMININ_LR.pdf", p5, width = 7, height = 9)


# 03.Specified Pathway Networks (Circle Plots)
# VEGF
pathways.show <- c("VEGF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
pdf(file = "Figures/Lab_cellchat2_netVisual_VEGF_pathway.pdf", width = 10, height = 7)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(
    object.list[[i]], 
    signaling = pathways.show, 
    layout = "circle", 
    edge.weight.max = weight.max[1], 
    edge.width.max = 10, 
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}
dev.off()

# LAMININ
pathways.show <- c("LAMININ") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
pdf(file = "Figures/Lab_cellchat2_netVisual_LAMININ_pathway.pdf", width = 10, height = 7)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(
    object.list[[i]], 
    signaling = pathways.show, 
    layout = "circle", 
    edge.weight.max = weight.max[1], 
    edge.width.max = 10, 
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}
dev.off()


# Pairwise Pathway Strength Visualization --------------------------------------
Pathway.df <- subsetCommunication(cellchat_WT, slot.name = "netP")

plot_pathway_dotplot <- function(df, src, tgt, title_text) {
  ggplot(df %>% dplyr::filter(source %in% src & target %in% tgt),
         aes(target, pathway_name, size = -log10(pval + 0.001), color = prob)) +
    geom_point() +
    scale_color_gradient2(high = "#d6604d", mid = '#92c5de', low = "white") + 
    theme_bw() + 
    scale_size_continuous(range = c(1, 3)) +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      aspect.ratio = 6
    ) +
    ggtitle(title_text)
}

p6 <- plot_pathway_dotplot(Pathway.df, "SynTII", "Endothelial", "SynTII to Endo Pathway")
p7 <- plot_pathway_dotplot(Pathway.df, "Endothelial", "SynTII", "Endo to SynTII Pathway")
p8 <- plot_pathway_dotplot(Pathway.df, "SynTI", "SynTII", "SynTI to SynTII Pathway")
p9 <- plot_pathway_dotplot(Pathway.df, "SynTII", "SynTI", "SynTII to SynTI Pathway")
p10 <- plot_pathway_dotplot(Pathway.df, "S-TGC", "SynTI", "S-TGC to SynTI Pathway")
p11 <- plot_pathway_dotplot(Pathway.df, "SynTI", "S-TGC", "SynTI to S-TGC Pathway")

ggsave("Figures/Lab_cellchat2_WT_SynTII2Endo_pathway.pdf", p6, width = 5, height = 4)
ggsave("Figures/Lab_cellchat2_WT_Endo2SynTII_pathway.pdf", p7, width = 5, height = 4)
ggsave("Figures/Lab_cellchat2_WT_SynTI2SynTII_pathway.pdf", p8, width = 5, height = 4)
ggsave("Figures/Lab_cellchat2_WT_SynTII2SynTI_pathway.pdf", p9, width = 5, height = 4)
ggsave("Figures/Lab_cellchat2_WT_STGC2SynTI_pathway.pdf", p10, width = 5, height = 4)
ggsave("Figures/Lab_cellchat2_WT_SynTI2STGC_pathway.pdf", p11, width = 5, height = 4)
