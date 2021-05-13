library(rtracklayer)
library(dplyr)
library("readxl")
library(tidyverse)
library(tidyfst)
library(tibble)

if (!interactive()) 
{
  gtf_df <- as.data.frame(rtracklayer::import('data/gencode.v27.annotation.gtf.gz'))
  
  expression_reads <- read.table("data/expression.matrix.tx.numreads.tsv", sep = '\t' , header = FALSE , stringsAsFactors = FALSE)
  
  phenotypes <- read.csv("data/phenotype.csv", sep = ",", header= TRUE)
  
  # phenotype_groups <- phenotypes %>% group_by(phenotype) %>% summarise(groups = list(Sample.ID))
  # TODO : filter only selected two groups  and make a proper class labels
  
  transcripts_genes <- expression_reads %>% 
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
  expression_reads_flipped <- expression_reads
  rownames(expression_reads_flipped) <- expression_reads_flipped[,1] # add transcriptIDs as rownames
  expression_reads_flipped <- expression_reads_flipped[,-1]
  transpose_tpms <- t_dt(expression_reads_flipped)
  
  for(i in 1:nrow(transcripts_genes))
  {
    selected_row <- transcripts_genes[i,]
    selected_gene_id<- selected_row$gene_id
    selected_tids<- selected_row$transcript_ids[[1]]
    selected_tid_tpms <- transpose_tpms %>% select(all_of(selected_tids))
    
    # Perform manova on it
  }
  
}
