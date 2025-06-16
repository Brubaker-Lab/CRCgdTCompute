library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(car)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
# 
# BiocManager::install("GenomicRanges")
# BiocManager::install("rtracklayer")
library(GenomicRanges)
library(rtracklayer)

setwd('F:/gdT_aim2/NicheNet')
seuratObj <- readRDS('F:/gdT_aim2/organized_data/unimputed/gdT_integrated.rds')
# metadata = read.table('F:/gdT_aim2/organized_data/unimputed/final_annotation.csv', sep = ',',header = 1)
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

predicted_gr_network = read.table('F:/gdT_aim2/SCENIC/reg_pairs_cancer_only.csv',sep = ',', header = 1,row.names = 1)
gr_network <- rename(gr_network,  tf = 'from', target = 'to')
result <- left_join(predicted_gr_network, gr_network, by = c("tf", "target")) # this is the step to see which of the predicted network has NicheNet source
print(result)
df_unique <- result %>% distinct(tf, target, .keep_all = TRUE) # to see which distinct tf-target has source. This is exactly the unique(predicted_gr_network) but with source
print(df_unique)
dim(df_unique)
no_source_predicted_grn = df_unique[is.na(df_unique$source),]
dim(no_source_predicted_grn)

write.table(df_unique, file = "F:/gdT_aim2/SCENIC/predicted_GRN.csv", sep = ",")
write.table(no_source_predicted_grn, file = "F:/gdT_aim2/SCENIC/predicted_GRN_no_source.csv", sep = ",")

# Calculate the proportions of NA and non-NA
source_data <- df_unique %>%
  mutate(has_source = ifelse(is.na(source), "No Source", "Has Source")) %>%
  group_by(has_source) %>%
  summarise(count = n()) %>%
  ungroup() %>% 
  mutate(perc = count/sum(count)) %>% 
  mutate(labels = scales::percent(perc))

# Create the pie chart
piechart = ggplot(source_data, aes(x = "", y = perc, fill = has_source)) +
  geom_bar(stat = "identity", width = 1) +
  geom_col() +
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5), cex = 9) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(fill = "Source Status", 
       title = "Proportion of TF-Target pairs with and without Source Information")
ggsave("vis/SCENIC_compare.jpg", plot = piechart, width = 5, height = 7,dpi = 300)