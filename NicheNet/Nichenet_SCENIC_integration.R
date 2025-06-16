library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(car)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
# BiocManager::install("GenomicRanges")
# BiocManager::install("rtracklayer")
library(GenomicRanges)
library(rtracklayer)
library(ggbeeswarm)
library(cowplot)

setwd('F:/gdT_aim2/NicheNet')
seuratObj <- readRDS('F:/gdT_aim2/organized_data/unimputed/gdT_integrated.rds')
Idents(seuratObj) <- seuratObj@meta.data$cell.type

seuratObj@meta.data %>% head()
seuratObj@meta.data$cell.type %>% table()

# seuratObj <- NormalizeData(seuratObj)
# saveRDS(seuratObj , 'AA_Done/data_organized/gdT_integrated_treatmentnaive_count.rds')

lr_network <- readRDS('lr_network_human_21122021.rds')
ligand_target_matrix <- readRDS('ligand_target_matrix_nsga2r_final.rds')
weighted_networks <- readRDS('weighted_networks_nsga2r_final.rds')
gr_network = readRDS("gr_network_human_21122021.rds")
sig_network = readRDS("signaling_network_human_21122021.rds")

combined_lr= readRDS("combined_lr_cancer_only.rds")
combined_grn = readRDS("combined_grn_cancer_only.rds")
# volcano -----------------------------------------------------------------
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.5) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(new_ligand_target_matrix)]

# only part of the labels are shown
volcano_labels <- ifelse(DE_table_receiver$p_val_adj <= 0.05 & abs(DE_table_receiver$avg_log2FC) > 3.2, DE_table_receiver$gene, '')
volcano_labels[is.na(volcano_labels)] <- ""
volcano_labels <- ifelse(!startsWith(DE_table_receiver$gene, 'IGB'), volcano_labels, '')
volcano_labels <- ifelse(!startsWith(DE_table_receiver$gene, 'RP'), volcano_labels, '')

volcanoplot <- EnhancedVolcano(
  DE_table_receiver,
  lab = '',  # No labels here to avoid duplication
  x = 'avg_log2FC',
  y = 'p_val_adj',
  ylab = bquote(~-Log[10] ~ italic(Padj)),
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 1.0,
  labSize = 0.5,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  max.overlaps = 30
)

volcanoplot <- volcanoplot + geom_text_repel(
  aes(label = volcano_labels),
  size = 2,  # Adjust size as needed
  segment.color = 'grey50',
  color = 'black',
  bg.color = 'cyan',  # Background color for the outline
  bg.r = 0.1,  # Radius for the outline
  max.overlaps = 70,  # Allow more overlaps
  min.segment.length = 0.5  # Minimum segment length for label connectors
)
volcanoplot

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(new_ligand_target_matrix)]

length(background_expressed_genes)
length(geneset_oi)

ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = new_ligand_target_matrix,
                                               potential_ligands = potential_ligands)

#ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

sorted_table <- ligand_activities[order(ligand_activities$aupr_corrected,decreasing = TRUE), ]
sorted_table

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange",bins = 60)  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity


all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
ligand_activities_no_scenic <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

sorted_table_no_scenic <- ligand_activities_no_scenic[order(ligand_activities_no_scenic$aupr_corrected,decreasing = TRUE), ]
sorted_table_no_scenic

# if (!require("pheatmap")) install.packages("pheatmap")
# if (!require("ComplexHeatmap")) install.packages("ComplexHeatmap", dependencies = TRUE)

library(pheatmap)
library(ComplexHeatmap)

get_top_genes <- function(matrix, ligand) {
  gene_scores <- matrix[ligand, ]
  top_genes <- names(sort(gene_scores, decreasing = TRUE)[1:5])
  return(top_genes)
}

ligand_list = sorted_table$test_ligand[1:5]
top_genes_per_ligand <- sapply(ligand_list, function(l) get_top_genes(ligand_target_matrix, l))
DotPlot(seuratObj, features = top_genes_per_ligand[,'CDHR2'], group.by = "cell.type") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

unique_genes <- unique(as.vector(top_genes_per_ligand))
heatmap_matrix <- ligand_target_matrix[ligand_list,unique_genes]

png("Teff2Tex_Nichenet.png", width = 1000, height = 300)
Heatmap(heatmap_matrix, name = "Interaction Strength", 
        cluster_rows = TRUE, cluster_columns = TRUE, 
        show_row_names = TRUE, show_column_names = TRUE,
        col = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# Trainig new nichenet with SCENIC ----------------------------------------
library(nichenetr)
library(tidyverse)

# tf_unique = read.table("predicted_GRN.csv", sep = ",")
tf_unique = read.table("F:/gdT_aim2/SCENIC/reg_pairs_cancer_only.csv", sep = ",", header = 1,row.names = 1)
tf_unique <- tf_unique %>% distinct(tf, target)
tf_unique$enrichment = NULL
tf_unique$source = 'SCENIC'
tf_unique$database = 'SCENIC'
colnames(tf_unique) = c('from','to','source','database')
combined_grn <- bind_rows(gr_network, tf_unique)
saveRDS(combined_grn,'combined_grn_cancer_only.rds')

lit_lr = read.table("lr_add.csv", sep = ",",header = 1)
combined_lr =  bind_rows(lr_network, lit_lr)
saveRDS(combined_lr,'combined_lr_cancer_only.rds')

new_network_weights_df = tibble(source = "SCENIC", avg_weight = 1, median_weight = 1)
new_source_weights_df = optimized_source_weights_df %>% bind_rows(new_network_weights_df)
new_network_weights_df = tibble(source = "literature", avg_weight = 1, median_weight = 1)
new_source_weights_df = new_source_weights_df %>% bind_rows(new_network_weights_df)

# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
source_weights_df_input = new_source_weights_df %>% select(source, avg_weight)
colnames(source_weights_df_input) = c("source","weight")
new_weighted_networks = construct_weighted_networks(lr_network = combined_lr, sig_network = sig_network, gr_network = combined_grn,
                                                source_weights_df = source_weights_df_input)

# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
new_weighted_networks = apply_hub_corrections(weighted_networks = new_weighted_networks,
                                          lr_sig_hub = hyperparameter_list %>% filter(parameter == "lr_sig_hub") %>% pull(avg_weight),
                                          gr_hub = hyperparameter_list %>% filter(parameter == "gr_hub") %>% pull(avg_weight))

# Infer ligand-target regulatory potential scores based on the weighted integrated networks
ligands = as.list(unique(combined_lr$from))
new_ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = new_weighted_networks, ligands = ligands, algorithm = "PPR",
                                                      damping_factor = hyperparameter_list %>% filter(parameter == "damping_factor") %>% pull(avg_weight),
                                                      ltf_cutoff = hyperparameter_list %>% filter(parameter == "ltf_cutoff") %>% pull(avg_weight))
saveRDS(new_ligand_target_matrix,'new_ligand_target_matrix_cancer_only.rds')
saved_new_ligand_target_matrix = readRDS('new_ligand_target_matrix_cancer_only.rds')

# check if saved = calculated
saved_new_ligand_target_matrix['WNT5A',1:3]
new_ligand_target_matrix['WNT5A',1:3]

# Validate the withSCENIC performance -------------------------------------
expression_settings_validation = readRDS("evaluation_final/evaluation/ligand_treatment_datasets/expression_settings")

# source_names_all = new_source_weights_df$source %>% unique()
# source_names = c(combined_lr$source, sig_network$source, combined_grn$source) %>% unique()
#source_weights_df = tibble(source = source_names, weight = rep(1, times = length(source_names))) # default for unoptimized model!

# > load("F:/CRC/Nichenet/hyperparameter_list.rda")
# Warning message:
#   R graphics engine version 16 is not supported by this version of RStudio. The Plots tab will be disabled until a newer version of RStudio is installed. 
# > hyperparameter_list
# # A tibble: 4 × 3
# parameter      avg_weight median_weight
# <chr>               <dbl>         <dbl>
#   1 damping_factor     0.789         0.830 
# 2 gr_hub             0.0803        0.0525
# 3 lr_sig_hub         0.115         0.0850
# 4 ltf_cutoff         0.926         0.920 

load("hyperparameter_list.rda")

parameters_setting = list()
parameters_setting$model_name = "with_scenic_prediction"
parameters_setting$source_weights = new_source_weights_df$avg_weight %>% magrittr::set_names(new_source_weights_df$source) 
parameters_setting$lr_sig_hub =  as.numeric(hyperparameter_list[hyperparameter_list['parameter'] =='lr_sig_hub','avg_weight'])
parameters_setting$gr_hub = as.numeric(hyperparameter_list[hyperparameter_list['parameter'] =='gr_hub','avg_weight'])
parameters_setting$algorithm = "PPR"
parameters_setting$correct_topology = FALSE
parameters_setting$damping_factor =  as.numeric(hyperparameter_list[hyperparameter_list['parameter'] =='damping_factor','avg_weight'])
parameters_setting$ltf_cutoff =   as.numeric(hyperparameter_list[hyperparameter_list['parameter'] =='ltf_cutoff','avg_weight'])

model_performances = evaluate_model(parameters_setting = parameters_setting,
                                       lr_network = combined_lr,
                                       sig_network = sig_network,
                                       gr_network = combined_grn,
                                       settings = lapply(expression_settings_validation,convert_expression_settings_evaluation),
                                       calculate_popularity_bias_ligand_prediction = FALSE, 
                                       calculate_popularity_bias_target_prediction = FALSE,
                                       ncitations = ncitations, 
                                       secondary_targets = FALSE, 
                                       remove_direct_links = "no")

saveRDS(model_performances,"model_performances_with_scenic_prediction")
saved_model_performances = readRDS("model_performances_with_scenic_prediction")

model_performances$performances_target_prediction 

saved_model_performances$performances_target_prediction[1,]
model_performances$performances_target_prediction[1,]

model_performances$performances_ligand_prediction_single

#BTW, source code https://github.com/saeyslab/nichenetr/blob/master/R/evaluate_model_target_prediction.R

# Test out the unchanged ----------------------------------------------------
# list_performances_random = readRDS("evaluation_final/evaluation/evaluation/results/100_random_networks_performance")
# #list_performances_unoptimized_ovo = readRDS("../characterization/results/list_performances_ovo_unoptimized_df05")
# list_performances_unoptimized_all = readRDS("evaluation_final/evaluation/evaluation/results/model_performances_unoptimized_all_sources")
# list_performances = readRDS("evaluation_final/evaluation/evaluation/results/list_performances_cv")
#list_performances_ipa = readRDS("results/fis_res_ipa")
#list_performances_ccce = readRDS("results/fis_res_ccce")
source_weights_df_input_old = optimized_source_weights_df %>% select(source, avg_weight)
colnames(source_weights_df_input_old) = c("source","weight")

weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network,
                                                source_weights_df = source_weights_df_input_old)
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks,
                                          lr_sig_hub = hyperparameter_list %>% filter(parameter == "lr_sig_hub") %>% pull(avg_weight),
                                          gr_hub = hyperparameter_list %>% filter(parameter == "gr_hub") %>% pull(avg_weight)) 
ligands = as.list(unique(lr_network$from))
ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR",
                                                      damping_factor = hyperparameter_list %>% filter(parameter == "damping_factor") %>% pull(avg_weight),
                                                      ltf_cutoff = hyperparameter_list %>% filter(parameter == "ltf_cutoff") %>% pull(avg_weight))

saved_ligand_target_matrix <- readRDS('ligand_target_matrix_nsga2r_final.rds')
#confirmed that saved ligand target matrix is the one with the current hyperparameters
saved_ligand_target_matrix['IL17A',1:3]
ligand_target_matrix['IL17A',1:3]
new_ligand_target_matrix['IL17A',1:3]

# source_names_all = source_weights_df$source %>% unique()
# source_names = c(lr_network$source, sig_network$source, gr_network$source) %>% unique()
source_weights_df = optimized_source_weights_df

parameters_setting_old = list()
parameters_setting_old$model_name = "without_scenic_prediction"
parameters_setting_old$source_weights = source_weights_df$avg_weight %>% magrittr::set_names(source_weights_df$source) 
parameters_setting_old$lr_sig_hub =  as.numeric(hyperparameter_list[hyperparameter_list['parameter'] =='lr_sig_hub','avg_weight'])
parameters_setting_old$gr_hub = as.numeric(hyperparameter_list[hyperparameter_list['parameter'] =='gr_hub','avg_weight'])
parameters_setting_old$algorithm = "PPR"
parameters_setting_old$correct_topology = FALSE
parameters_setting_old$damping_factor =  as.numeric(hyperparameter_list[hyperparameter_list['parameter'] =='damping_factor','avg_weight'])
parameters_setting_old$ltf_cutoff =   as.numeric(hyperparameter_list[hyperparameter_list['parameter'] =='ltf_cutoff','avg_weight'])

old_model_performances = evaluate_model(parameters_setting = parameters_setting_old,
                                           lr_network = lr_network,
                                           sig_network = sig_network,
                                           gr_network = gr_network,
                                           settings = lapply(expression_settings_validation,convert_expression_settings_evaluation),
                                           calculate_popularity_bias_ligand_prediction = FALSE,
                                           calculate_popularity_bias_target_prediction = FALSE,
                                           ncitations = ncitations, 
                                           secondary_targets = FALSE,
                                           remove_direct_links = "no",
                                           mc.cores =4)
saveRDS(old_model_performances,"model_performances_without_scenic_prediction")
old_model_performances = readRDS("model_performances_without_scenic_prediction")

# visualize the performance -----------------------------------------------
selected_performances = bind_rows(
  old_model_performances$performances_target_prediction %>% select(-ligand) %>% mutate(model = "Vanilla Nichnet")
) %>% bind_rows(
  model_performances$performances_target_prediction %>% select(-ligand) %>% mutate(model = "Added SCENIC")
)

selected_performances = selected_performances %>% mutate(model = factor(model, levels = c("Vanilla Nichnet",
                                                                                          "Added SCENIC")))

selected_performances_dataset_median = selected_performances %>% group_by(model, setting) %>% summarise(auroc = median(auroc), aupr_corrected = median(aupr_corrected), auc_iregulon_corrected = median(auc_iregulon_corrected), pearson = median(pearson), spearman = median(spearman), mean_rank_GST_log_pval = median(mean_rank_GST_log_pval))
selected_performances_model_median = selected_performances %>% group_by(model) %>% summarise(auroc = median(auroc), aupr_corrected = median(aupr_corrected), auc_iregulon_corrected = median(auc_iregulon_corrected), pearson = median(pearson), spearman = median(spearman), mean_rank_GST_log_pval = median(mean_rank_GST_log_pval))

performances_dataset_median = selected_performances_dataset_median  %>% gather(key = scorename, value = scorevalue, auroc:mean_rank_GST_log_pval)
scorelabels = c(auroc="AUROC", aupr_corrected="AUPR (corrected)", auc_iregulon_corrected = "AUC-iRegulon (corrected)",pearson = "Pearson correlation", spearman = "Spearman's rank correlation",mean_rank_GST_log_pval = "Mean-rank gene-set enrichment")
scorerandom = c(auroc=0.5, aupr_corrected=0, auc_iregulon_corrected = 0, pearson = 0, spearman = 0,mean_rank_GST_log_pval = 0) %>% data.frame(scorevalue=.) %>% rownames_to_column("scorename")
performances_model_median = selected_performances_model_median %>% gather(key = scorename, value = scorevalue, auroc:mean_rank_GST_log_pval)
scorelabels = c(auroc="AUROC", aupr_corrected="AUPR (corrected)", auc_iregulon_corrected = "AUC-iRegulon (corrected)",pearson = "Pearson correlation", spearman = "Spearman's rank correlation",mean_rank_GST_log_pval = "Mean-rank gene-set enrichment")
scorerandom = c(auroc=0.5, aupr_corrected=0, auc_iregulon_corrected = 0, pearson = 0, spearman = 0,mean_rank_GST_log_pval = 0) %>% data.frame(scorevalue=.) %>% rownames_to_column("scorename")

target_prediction_performances_cv = performances_dataset_median %>% ggplot() + 
  geom_quasirandom(aes(model, scorevalue, group=model, color = model), alpha = 0.7) + 
  geom_point(data = performances_model_median, aes(model,scorevalue), shape = "_", size = 15, color = "black") + 
  scale_color_manual(values = c("darkgrey","lightblue","steelblue1","royalblue")) +
  facet_wrap(~scorename, scales = "free", labeller=as_labeller(scorelabels),nrow = 3) +
  scale_y_continuous("Score target prediction") + 
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())+ 
  theme(legend.position="bottom", legend.title = element_blank()) + 
  geom_hline(aes(yintercept=scorevalue), data=scorerandom, linetype = 2, color = "red") 

target_prediction_performances_cv