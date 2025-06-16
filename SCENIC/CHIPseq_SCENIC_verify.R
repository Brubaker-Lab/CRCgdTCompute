library(GenomicRanges)
library(rtracklayer)

setwd("F:/gdT_aim2/SCENIC/")
bed_dir    <- "F:/gdT_aim2/SCENIC/SCENIC_database"
bed_files  <- list.files(bed_dir, pattern="\\.bed$", full.names=TRUE)

gr_list <- lapply(bed_files, function(f){
  # 1) celltype from the file name
  celltype <- sub("\\.bed$", "", basename(f))
  celltype <- sub("^Oth\\.Bld\\.05\\.AllAg\\.", "", celltype)  
  
  # 2) import the peaks
  gr <- import(f, format="BED")
  
  # 3) split that big metadata string on “;”,
  #    grab the element that starts with “Name=”
  attrs     <- strsplit(mcols(gr)$name, ";")
  name_field <- sapply(attrs, function(x){
    nm <- grep("^Name=", x, value=TRUE)
    sub("^Name=", "", nm)
  })
  
  # 4) URL‑decode it
  name_decoded <- URLdecode(name_field)
  # 5) TF is the very first token (everything before space or “(”)
  mcols(gr)$TF       <- sub("([^( ]*).*", "\\1", name_decoded)
  # 6) tag cell type
  mcols(gr)$celltype <- celltype
  gr
})

# combine & sort
all_gr <- do.call(c, gr_list)
all_gr <- sort(all_gr)

# export(all_gr,file="all_celltypes_sorted.bed",format="BED")
saveRDS(all_gr,"all_celltypes_sorted_gr_bed")

target_tss_dir <- "F:/gdT_aim2/SCENIC/SCENIC_database/TF_target_genes_tss"
tss_files      <- list.files(target_tss_dir,
                             pattern="no_source_tss_5kb_for_.*\\.bed$",
                             full.names=TRUE)

# get the full list of cell‑types in your ChIP data
all_celltypes <- unique(mcols(all_gr)$celltype)

results_list <- list()
verified_target_percentage = list()

for(tss_file in tss_files){
  tf_name <- sub("\\.bed$","",
                 sub("no_source_tss_5kb_for_","",
                     basename(tss_file)))
  message("TF = ", tf_name)
  
  ## 1) import TSS and get the **unique** gene list
  tss_gr      <- import(tss_file, format="BED")
  gene_levels <- sort(unique(mcols(tss_gr)$name))
  
  ## 2) extract that TF's peaks
  tf_peaks <- all_gr[mcols(all_gr)$TF == tf_name]
  if(length(tf_peaks)==0){
    warning("  no ChIP peaks for ", tf_name, "; skipping\n")
    next
  }
  
  ## 3) find overlaps
  ov      <- findOverlaps(tf_peaks, tss_gr)
  df_hits <- data.frame(
    celltype = mcols(tf_peaks)$celltype[queryHits(ov)],
    gene     = mcols(tss_gr)$name     [subjectHits(ov)]
  )
  
  ## 4) force your factors to use **all** levels
  df_hits$celltype <- factor(df_hits$celltype,
                             levels = all_celltypes)
  df_hits$gene     <- factor(df_hits$gene,
                             levels = gene_levels)
  
  ## 5) tabulate into a full matrix (zeros everywhere else)
  tab    <- table(df_hits$celltype, df_hits$gene)
  df_tab <- as.data.frame.matrix(tab)        # rows=celltypes, cols=genes
  
  ## 6) (optional) push the rownames into a proper column
  # df_tab$celltype <- rownames(df_tab)
  # rownames(df_tab) <- NULL
  verified_target_percentage[[tf_name]] <-sum(colSums(df_tab)!=0)/dim(df_tab)[2]

  ## 7) save it
  results_list[[tf_name]] <- df_tab
  write.table(df_tab,
            file = paste0(tf_name, "_celltype_gene_peakCounts.csv"),
            sep = ',')
}

write.table(as.data.frame(verified_target_percentage), 
            sep = ',',
            file = 'verified_target_percentage.csv')
