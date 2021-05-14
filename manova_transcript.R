library(rtracklayer)
library(dplyr)
library("readxl")
library(tidyverse)
library(tidyfst)
library(tibble)
library(micompr)


if (!interactive()) 
{
  gtf_df <- as.data.frame(rtracklayer::import('data/gencode.v27.annotation.gtf.gz'))
  
  phenotypes <- read.csv("data/phenotype.csv", sep = ",", header= TRUE)
  classLabels <- phenotypes$phenotype
  selected_samples <- phenotypes$Sample.ID
  
  expression_reads <- read.table("data/expression.matrix.tx.numreads.tsv", sep = '\t' , header = TRUE , stringsAsFactors = FALSE)
  rownames(expression_reads) <- expression_reads[,1] # add transcriptIDs as rownames
  
  # get only samples given in phenotypeCSV
  expression_reads_filtered <- expression_reads %>% select(Gene_or_Transcript_ID,all_of(selected_samples))
  # take out transcripts with zero read count on all samples
  # expression_reads_filtered <- expression_reads_filtered[rowSums(expression_reads_filtered[,-1])>0,] 
  
  sampleSize <- ncol(expression_reads_filtered) -1 # take out gene_id column
  
  expression_reads_filtered <- expression_reads_filtered %>% 
    mutate(ZeroCounts =  rowSums(.[,-1]== 0) ) %>% 
    filter(ZeroCounts < sampleSize * 0.9) %>% 
    select(-ZeroCounts)
  
  transcripts_genes <- expression_reads_filtered %>% 
    select(Gene_or_Transcript_ID) %>% #get the transcript column
    left_join(gtf_df %>% 
                select(gene_id, transcript_id ) %>% 
                distinct(gene_id,transcript_id) , by = c("Gene_or_Transcript_ID" = "transcript_id")) %>% 
    group_by(gene_id) %>% 
    summarise(transcript_ids = list(Gene_or_Transcript_ID)) %>%  # add list of transcript ids for each gene
    filter(!is.na(gene_id)) # filter out the na group
  # There are two transcripts without gene_id - Check for later
  # gene_id transcript_ids   
  # <chr>   <chr>            
  # 1 NA      ENST00000360403.2
  # 2 NA      ENST00000372183.3
  
  
  # transpose matrix and keep transcriptIDs as colnames and add seperate sampleID column
  expression_reads_filtered <- expression_reads_filtered[,-1] # Take out geneID column
  transpose_tpms <- t_dt(expression_reads_filtered)

  output <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("gene_id", "manova_p_value"))
  
  for(i in 1:nrow(transcripts_genes))
  {
    selected_row <- transcripts_genes[i,]
    selected_gene_id<- selected_row$gene_id
    selected_tids<- selected_row$transcript_ids[[1]]
    selected_tid_tpms <- transpose_tpms %>% select(all_of(selected_tids))
    
    transcript_size_gene <- length(selected_tids)
    # Perform manova on it
    cmp <- cmpoutput("scRNA-seq", transcript_size_gene , selected_tid_tpms, as.factor(classLabels))
    output[i,] = c(selected_gene_id, cmp$p.values$manova)
  }
  
  output[,2] <- sapply(output[,2] , as.numeric) # convert p values to numeric
  
  write.csv(output_converted %>% arrange(manova_p_value)  , "manova-outputs.csv")
}
