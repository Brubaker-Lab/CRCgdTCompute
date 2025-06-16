library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(car)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(stringr)
# BiocManager::install("GenomicRanges")
# BiocManager::install("rtracklayer")
library(GenomicRanges)
library(rtracklayer)

# I am doing the NicheNet on the cells from all -- normal, cancer, adenoma
# But I am going to use the cells in the cancer tissue for enrichment analysis?

setwd('F:/gdT_aim2/NicheNet')
seuratObj <- readRDS('F:/gdT_aim2/organized_data/unimputed/gdT_integrated.rds')
# metadata = read.table('F:/gdT_aim2/organized_data/unimputed/final_annotation.csv', sep = ',',header = 1)
general_type = seuratObj@meta.data$cell.type
general_type <- str_replace_all(general_type, "IL7R.* TRM", "TRM")
general_type <- str_replace_all(general_type, "^.*Teff.*", "Teff")
general_type <- str_replace_all(general_type, "^Tex.*", "Tex")

seuratObj@meta.data$general_type = general_type
# Idents(seuratObj) <- seuratObj@meta.data$cell.type
Idents(seuratObj) <- seuratObj@meta.data$general_type

# seuratObj@meta.data %>% head()
# seuratObj@meta.data$cell.type %>% table()
# seuratObj <- NormalizeData(seuratObj)
# saveRDS(seuratObj , 'AA_Done/data_organized/gdT_integrated_treatmentnaive_count.rds')

seuratObj_subset <- subset(seuratObj, subset = tissue == "Carcinoma")
# we only take the subset here because fibroblast likely blocks off the interaction of normal gdT to cancer cells to some extent

lr_network <- readRDS('lr_network_human_21122021.rds')
ligand_target_matrix <- readRDS('ligand_target_matrix_nsga2r_final.rds')
weighted_networks <- readRDS('weighted_networks_nsga2r_final.rds')
gr_network = readRDS("gr_network_human_21122021.rds")
sig_network = readRDS("signaling_network_human_21122021.rds")
combine_lr = readRDS("combined_lr.rds")
combined_grn = readRDS("combined_grn.rds")
new_ligand_target_matrix=readRDS('new_ligand_target_matrix_cancer_only.rds')
# upregulated_genes_list = readRDS('deg_results_unimputed_1vsALL')

# receiver_pairs <- combn(unique(Idents(seuratObj)), 2, simplify = FALSE)
# receiver_pairs = list(c('Tex','Teff'),c('Teff','TRM')) #two transitions that we are interested in: Teff ->Tex and TRM -> Teff
# 
# deg_results <- lapply(receiver_pairs, function(pair) {
# 
#   seurat_obj_receiver <- subset(seuratObj_subset, idents = pair)
#   seurat_obj_receiver = PrepSCTFindMarkers(seurat_obj_receiver)
#   FindMarkers(object = seurat_obj_receiver, ident.1 = pair[1], ident.2 = pair[2], group.by = "general_type", min.pct = 0.05, assay = "SCT") %>%
#     rownames_to_column("gene") %>%
#     mutate(group1 = pair[1], group2 = pair[2])
# })
# 
# deg_results <- bind_rows(deg_results)
# saveRDS(deg_results,'deg_results_unimputed_interest')

# markers <- FindAllMarkers(
#   seuratObj_subset, 
#   only.pos = TRUE,         # Only return upregulated genes
#   min.pct = 0.05,           # Minimum 10% of cells expressing the gene
#   logfc.threshold = 1,     # log2FC > 1
#   test.use = "wilcox"      # Wilcoxon Rank Sum test (default)
# )
# filtered_markers <- markers %>%
#   filter(p_val_adj < 0.05, avg_log2FC > 1)
# upregulated_genes_list <- split(filtered_markers, filtered_markers$cluster)
# saveRDS(upregulated_genes_list,'deg_results_unimputed_1vsALL_canceronly')

#after doing quick DEG in python, see other script
TRM2Teff = read.csv('TRM2Teff_up.csv',row.names = 1)
Teff2Tex = read.csv('Teff2Tex_up.csv',row.names = 1)
upregulated_genes_list <- list()
upregulated_genes_list$TRM2Teff = row.names(TRM2Teff)
upregulated_genes_list$Teff2Tex = row.names(Teff2Tex)
receivers = list()
receivers$TRM2Teff = 'TRM'
receivers$Teff2Tex = 'Teff'
goals = list()
goals$TRM2Teff = 'Teff'
goals$Teff2Tex = 'Tex'
# 
# for (celltype in names(upregulated_genes_list)) {
# geneset_oi <- upregulated_genes_list[[celltype]]
# volcano_labels <- ifelse(geneset_oi$pval_adj <= 0.05 & abs(geneset_oi$log2fold) > 2, rownames(geneset_oi), '')
# volcano_labels[is.na(volcano_labels)] <- ""
# volcano_labels <- ifelse(!startsWith(string_array, 'ENSG'), volcano_labels, '')
# 
# volcanoplot <- EnhancedVolcano(
#   res_for_vp,
#   lab = '',  # No labels here to avoid duplication
#   x = 'log2FoldChange',
#   y = 'padj',
#   ylab = bquote(~-Log[10] ~ italic(Padj)),
#   pCutoff = 0.05,
#   FCcutoff = 1,
#   pointSize = 1.0,
#   labSize = 0.5,
#   drawConnectors = TRUE,
#   widthConnectors = 0.5,
#   max.overlaps = 30
# )
# 
# volcanoplot <- volcanoplot + geom_text_repel(
#   aes(label = volcano_labels),
#   size = 2,  # Adjust size as needed
#   segment.color = 'grey50',
#   color = 'black',
#   bg.color = 'cyan',  # Background color for the outline
#   bg.r = 0.1,  # Radius for the outline
#   max.overlaps = 70,  # Allow more overlaps
#   min.segment.length = 0.5  # Minimum segment length for label connectors
# )
# 
# volcano_filename <- paste0( donor,  "_", res_name_sanitized, "EPIi_DesignCompare_volcanoplot.png")
# ggsave(volcano_filename, plot = volcanoplot, width = 5, height = 7)
# }
# with SCENIC -------------------------------------------------------------

ligand_activities_list <- list()
for (celltype in names(upregulated_genes_list)) {
  geneset_oi <- upregulated_genes_list[[celltype]]#$gene
  expressed_genes_receiver <- get_expressed_genes(receivers[[celltype]],seuratObj_subset, pct = 0.01, assay_oi = 'SCT')
  expressed_receptors <- intersect(unique(combined_lr$to), expressed_genes_receiver)

  potential_ligands <- combined_lr %>%
    filter(to %in% expressed_receptors) %>%
    pull(from) %>%
    unique()

  background_expressed_genes <- expressed_genes_receiver %>%
    .[. %in% rownames(new_ligand_target_matrix)]

  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = new_ligand_target_matrix,
    potential_ligands = potential_ligands
  )

  ligand_activities_list[[celltype]] <- ligand_activities
}

lapply(names(ligand_activities_list), function(celltype) {
  ligand_activities_list[[celltype]] = ligand_activities_list[[celltype]][order(ligand_activities_list[[celltype]]$auroc,decreasing = TRUE), ]
  write.csv(ligand_activities_list[[celltype]], paste0("withSCENIC/",celltype, "_ligand_activities_cancer_only.csv"), row.names = FALSE)
})

sorted_table <- ligand_activities_list[[celltype]][order(ligand_activities_list[[celltype]]$auroc,decreasing = TRUE), ]

# without SCENIC ----------------------------------------------------------
ligand_activities_list_withoutSCENIC <- list()
for (celltype in names(upregulated_genes_list)) {
  geneset_oi <- upregulated_genes_list[[celltype]]#$gene
  expressed_genes_receiver <- get_expressed_genes(receivers[[celltype]],seuratObj_subset, pct = 0.01, assay_oi = 'SCT')
  expressed_receptors <- intersect(unique(lr_network$to), expressed_genes_receiver)
  
  potential_ligands <- lr_network %>%
    filter(to %in% expressed_receptors) %>%
    pull(from) %>%
    unique()
  
  background_expressed_genes <- expressed_genes_receiver %>%
    .[. %in% rownames(new_ligand_target_matrix)]
  
  ligand_activities_withoutSCENIC <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  
  ligand_activities_list_withoutSCENIC[[celltype]] <- ligand_activities_withoutSCENIC
}

lapply(names(ligand_activities_list), function(celltype) {
  ligand_activities_list_withoutSCENIC[[celltype]] = ligand_activities_list_withoutSCENIC[[celltype]][order(ligand_activities_list_withoutSCENIC[[celltype]]$auroc,decreasing = TRUE), ]
  write.csv(ligand_activities_list_withoutSCENIC[[celltype]], paste0("withoutSCENIC/",celltype, "_ligand_activities_cancer_only.csv"), row.names = FALSE)
})

# Print out the target of inferred ligand ---------------------------------

for (celltype in names(upregulated_genes_list)) {
  ligand_activities = read.csv(paste0("withSCENIC/",celltype, "_ligand_activities_cancer_only.csv"))
  top_15_ligands = ligand_activities$test_ligand[1:15]
  target_gene_union <- as.data.frame( matrix(NA, nrow = 5, ncol = 15))
  colnames(target_gene_union) <- top_15_ligands
  for (ligand in top_15_ligands) {
    top_5_targets =new_ligand_target_matrix[,ligand][order(new_ligand_target_matrix[,ligand],decreasing = TRUE) ][1:5]
    target_gene_union[ligand] = names(top_5_targets)
  }
  all_targets <- unique(as.vector(as.matrix(target_gene_union)))
  write.table(all_targets,paste0(celltype,"_ligand_target.csv"),sep = ',')
  power_mat_to_plot = new_ligand_target_matrix[all_targets,top_15_ligands]
  
  # subset_obj_cancer_celltype <- subset(seuratObj_subset, idents = goals[[celltype]])

  # dotplot = DotPlot(subset_obj_cancer_celltype, features = all_targets, assay = 'SCT') +
    RotatedAxis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # ggsave(paste0('vis/',celltype,"_dotplot.jpg"), plot = dotplot, width = 12, height = 4,dpi = 300)
  
  # divisors <- unlist(apply(power_mat_to_plot, 2, max))
  # normalized_power_mat_to_plot <- sweep(power_mat_to_plot, 2, divisors, "/")
  
  # htmap = pheatmap(t(normalized_power_mat_to_plot),# Scale by row (optional, for better visualization)
           # cluster_rows = FALSE,   # Cluster rows
           # cluster_cols = TRUE,   # Cluster columns
           # color = colorRampPalette(c("blue", "white", "red"))(100),
           # show_rownames = TRUE,
           # show_colnames = TRUE)
  # ggsave(paste0('vis/',celltype,"_powermap.jpg"), plot = htmap, width = 12, height = 5,dpi = 300)
}
