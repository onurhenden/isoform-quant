if(!require(rtracklayer))
{
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rtracklayer")
}
if(!require(dplyr)) {
  install.packages("dplyr")
}
library(rtracklayer)
library(dplyr)

if(!file.exists("data/gencode.v27.annotation.gtf.gz")){
  download.file("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz", destfile = "data/gencode.v27.annotation.gtf.gz")
}

if (!interactive()) {
  gtf_df <- as.data.frame(rtracklayer::import('data/gencode.v27.annotation.gtf.gz'))
  tumor_genes <- c("GGCT", "H3F3A", "NME1", "PABPC1", "MYL12B", "CHCHD2", "HSP90AB1", "NR4A1", "TSTA3", "ZNF706", "SET", "SOD1", "ALDOA")
  myeloid_genes <- c("ATP6V0E1", "HLA-A", "SEP15", "RHOA", "B2M")
  tcell_genes <- c("SH3BGRL3", "PSMC2", "GABARAP", "WAC", "CFL1", "CLIC1", "TMEM126B", "CTSB", "RPL23A", "HSP90AA1", "LAP3", "DDX5", "RHOA", "B2M")
  stromal_genes <- c("TMEM59", "TSC22D1")
  bcell_genes <- c("HNRNPK", "MMADHC", "UBB", "ARPC2", "HNRNPA1", "EIF4A2", "CLK1", "ACTG1", "YWHAZ", "FKBP1A", "SUMO2", "BRK1", "EIF2S1", "GSPT1", "GPBP1L1", "CSNK1A1", "PFN1", "LTV1", "SNAP23", "ACTB", "TAF9", "DDX5", "LAP3")
  
  gene_list <- c(tumor_genes, myeloid_genes, tcell_genes, stromal_genes, bcell_genes)
  class_labels <- c( rep("tumor", length(tumor_genes)),
                     rep("myeloid", length(myeloid_genes)),
                     rep("tcell", length(tcell_genes)),
                     rep("stromal", length(stromal_genes)),
                     rep("bcell", length(bcell_genes)))
  
  genes_df <- data.frame(gene_list, class_labels)
  
  rm(tumor_genes,myeloid_genes,tcell_genes,stromal_genes,bcell_genes,gene_list,class_labels )
  
  
  gtf_df <- gtf_df %>% filter(gene_name %in% genes_df$gene_list)
  
  unique_transcript_ids <- gtf_df %>% distinct(transcript_id)
  
  expression_reads <- read.table("data/expression.matrix.tx.numreads.tsv", sep = '\t' , header = TRUE , stringsAsFactors = FALSE)
  
  expression_reads <- expression_reads  %>% 
    filter(Gene_or_Transcript_ID %in% c(unique_transcript_ids$transcript_id)) %>% 
    left_join(unique_transcript_ids, by = c("Gene_or_Transcript_ID" = "transcript_id"))
  
  expression_reads <- rbind(colnames(expression_reads) , expression_reads)
  expression_reads <- as.data.frame(t(expression_reads))
  names(expression_reads) <- as.matrix(expression_reads[1,])
  expression_reads <- expression_reads[-1,] 
  colnames(expression_reads)[1] = "Sample ID"
  
  #convert characters  to doubles
  expression_reads[,2:614] <- lapply(expression_reads[,2:614], function(x) as.numeric(replace(x, is.na(x), 0)))
  
  library("readxl")
  sampleid_disease_mapping <- read_excel("../bio-pipeline/Data/appended_Excel_corrected.xlsx",col_names = FALSE)
  sampleid_disease_mapping <- as.data.frame(t(sampleid_disease_mapping))
  names(sampleid_disease_mapping) <- as.matrix(sampleid_disease_mapping[1,])
  sampleid_disease_mapping <- sampleid_disease_mapping[-1,] 
  sampleid_disease_mapping <- sampleid_disease_mapping %>% select(`Sample ID`, `Patient ID`, `Cell Type`, Subtype)
  
  joined_data <- left_join(expression_reads, sampleid_disease_mapping, by="Sample ID")
  joined_data_subset <- joined_data %>% mutate(classLabel = ifelse(`Cell Type` == "Tumor", "Tumor", "NonTumor" )) %>% 
    select(-`Patient ID`,-`Cell Type`, -Subtype , -`Sample ID`) %>% 
    filter(!is.na(classLabel))#filtering out bulk samples
    
  
  mean_transcript_values <- joined_data_subset %>% group_by(classLabel) %>% summarise_all( mean)
  # median_transcript_values <- joined_data_subset %>% group_by(classLabel) %>% summarise_all( median)
  
  transcript_id_gene_mapping <- gtf_df %>% filter(type == "transcript") %>%  distinct(transcript_id, gene_name, gene_type,width)
  
  mean_transcript_values_t <- as.data.frame(t(rbind(colnames(mean_transcript_values) , mean_transcript_values)))
  names(mean_transcript_values_t) <- as.matrix(mean_transcript_values_t[1,])
  mean_transcript_values_t <- mean_transcript_values_t[-1,] 
  
  mean_transcript_values_t <- mean_transcript_values_t %>% dplyr::rename(transcript_id = classLabel)
  joined_mean_data <- left_join(mean_transcript_values_t, transcript_id_gene_mapping)

  ratio_mean_data <- joined_mean_data %>% 
    mutate(ratio = as.numeric(Tumor) /  as.numeric(NonTumor)) %>% 
    filter(as.numeric(NonTumor) >10 ) %>%  # Filter the low count
    filter(as.numeric(Tumor) >10 ) %>%  # Filter the low count
    group_by(gene_name) %>% mutate(avg_ratio = mean(ratio)) %>% ungroup() %>% 
    mutate(candidate = ifelse( avg_ratio > 1,
                               ifelse( ratio < avg_ratio / 2, TRUE, FALSE ),
                               ifelse( ratio > avg_ratio * 2 , ifelse(ratio>1, TRUE, FALSE), FALSE)))
    
  
  write.csv(ratio_mean_data %>% select(gene_name, transcript_id, gene_type, width, Tumor, NonTumor, ratio, avg_ratio, candidate), "ratio_mean_data.csv")
            
}
