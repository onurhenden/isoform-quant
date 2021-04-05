if(!require(rtracklayer))
{
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rtracklayer")
}
library(rtracklayer)

if(!require(dplyr)) {
  install.packages("dplyr")
}
library(dplyr)

if(!file.exists("data/gencode.v27.annotation.gtf.gz")){
  download.file("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz", destfile = "data/gencode.v27.annotation.gtf.gz")
}

if (!interactive()) {
  gtf_df <- as.data.frame(rtracklayer::import('data/gencode.v27.annotation.gtf.gz'))
  expression_reads <- read.csv("data/expression.matrix.tx.numreads.tsv", sep = '\t' , header = TRUE)
  
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
  
  
  
  subset_gtf <- gtf_df %>% filter(gene_name %in% genes_df$gene_list)
  unique_transcript_ids <- subset_gtf %>% distinct(transcript_id,gene_name)
  subset_expression_reads <- expression_reads  %>% 
    filter(Gene_or_Transcript_ID %in% c(unique_transcript_ids$transcript_id)) %>% 
    left_join(unique_transcript_ids, by = c("Gene_or_Transcript_ID" = "transcript_id"))
}
