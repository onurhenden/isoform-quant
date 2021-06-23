library(rtracklayer)
library(dplyr)
library("readxl")
library(tidyverse)
library(tidyfst)
library(tibble)
library(micompr)

boxplot_for_transcripts = function( gene_id ,tid_data , phenotype , sample_size)
{
  transposed_data <- as.data.frame(t(tid_data))
  # Log (TPM+1) transformation
  transposed_data <- log(transposed_data+1)
  
  pivoted_data <- transposed_data %>% 
    mutate(transcript_id = row.names(.)) %>% 
    tidyr::gather("Sample.ID", "log(TPM+1)", 1:sample_size)
  pivoted_data <- pivoted_data %>% left_join(phenotype)
  ggplot(pivoted_data, aes(x=transcript_id, y=TPM, fill=phenotype)) + 
    geom_boxplot() + ggtitle(gene_id) + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(filename = paste0(gene_id, ".png"), path= paste0("experiment-outputs/", experiment_name, "/boxplots") ) 
}

if (TRUE) 
{
  experimentList <- c(
    "bc03ln_tumor_vs_non_tumor"
    # "bc07ln_tumor_vs_non_tumor"
    # "er_positive_tumor_vs_non_tumor",
    # "her2_positive_tumor_vs_non_tumor"
  )
  
  
  gtf_df <- as.data.frame(rtracklayer::import('data/gencode.v27.annotation.gtf.gz'))
  expression_reads <- read.table("data/expression.matrix.tx.numreads.tsv", sep = '\t' , header = TRUE , stringsAsFactors = FALSE)
  rownames(expression_reads) <- expression_reads[,1] # add transcriptIDs as rownames
  
  for(i in 1:length(experimentList))
  {
    experiment_name <- experimentList[i]
    
    dir.create(paste0("experiment-outputs/", experiment_name))
    dir.create(paste0("experiment-outputs/", experiment_name, "/gene-tpms"))
    dir.create(paste0("experiment-outputs/", experiment_name, "/boxplots"))
    
    phenotypes <- read.csv(paste0("data/", experiment_name, ".csv"), sep = ",", header= TRUE)
    classLabels <- phenotypes$phenotype
    selected_samples <- phenotypes$Sample.ID
    
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
    nRowgenes  <- nrow(transcripts_genes)
    for(i in 1:nRowgenes)
    {
      # selected_row <- transcripts_genes %>% filter(gene_id == "ENSG00000115884.10")
      selected_row <- transcripts_genes[i,]
      selected_gene_id<- selected_row$gene_id
      selected_tids<- selected_row$transcript_ids[[1]]
      selected_tid_tpms <- transpose_tpms %>% select(all_of(selected_tids))
      
      #DTU FILTERING
      tpms_w_classlabels <- selected_tid_tpms %>% mutate(cls_lbl = classLabels)
      
      clusterCenters <- tpms_w_classlabels %>% group_by(cls_lbl) %>% summarise_all(list(median))
      
      clusterCenters <- rbind ( clusterCenters, data.frame(cls_lbl = "ratio" ,  clusterCenters[1,-1] / clusterCenters[2,-1] ))
      
      ratioList <- as.numeric( clusterCenters[1,-1] / clusterCenters[2,-1] )
      
      ratioList[is.nan(ratioList)] = -1 # replace NANS
      
      isoformSwitchFound = FALSE
      if(any(ratioList>1 ) && any(ratioList>0 && ratioList<1))
        isoformSwitchFound  = TRUE
      clusterCenters <- clusterCenters %>% mutate(totalSum = rowSums(across(where(is.numeric))))
      
      dtu_found <- (clusterCenters$totalSum[1] > 10 && clusterCenters$totalSum[2] > 10)
      # print(paste0("Selected gene -> ", selected_gene_id))
      if(dtu_found && isoformSwitchFound)
      {
        print(paste0("DTU found for gene -> ", selected_gene_id))
        transcript_size_gene <- length(selected_tids)
        # Perform manova on it
        cmp <- cmpoutput("scRNA-seq", transcript_size_gene , selected_tid_tpms, as.factor(classLabels))
        # print( paste0(i, " - ", nRowgenes))
        # rbindlist(list(output,c(selected_gene_id, cmp$p.values$manova) ))
        output[i,] = c(selected_gene_id, cmp$p.values$manova)
        print(paste0(cmp$p.values$manova))
        
        if(!is.na(cmp$p.values$manova) && cmp$p.values$manova < 0.05 && isoformSwitchFound)
        {
          boxplot_for_transcripts(selected_gene_id,selected_tid_tpms,phenotypes,sampleSize)
          write.csv( selected_tid_tpms %>% mutate(cls_lbl = classLabels) , paste0("experiment-outputs/", experiment_name,"/gene-tpms/" , selected_gene_id,".csv"))
          print(clusterCenters)
        }
        
      }
      
    }
    
    output[,2] <- sapply(output[,2] , as.numeric) # convert p values to numeric
    
    output <- output %>% arrange(manova_p_value)
    
    # dir.create(paste0("experiment-outputs/", experiment_name))
    # dir.create(paste0("experiment-outputs/", experiment_name, "/gene_tpms"))
    
    write.csv( output , paste0("experiment-outputs/", experiment_name, "/manova_outputs.csv"))
    
    # class_colors <- phenotypes %>% mutate(cLr = ifelse(phenotype == "Tumor", "red","blue"))
    # class_colors <- class_colors$cLr
    # library(Rtsne) # Load package
    # topValues = 10
    # 
    # for ( i in 1:topValues )
    # {
    #   gene_to_plot <- output[i,1]
    #   png(paste0("experiment-outputs/", experiment_name,"/",i,"th-" ,gene_to_plot ,".png"))
    #   txs <- transcripts_genes %>% filter(gene_id == gene_to_plot)
    #   selected_tids<- txs$transcript_ids[[1]]
    #   selected_tid_tpms <- transpose_tpms %>% select(all_of(selected_tids))
    #   tsne_out <- Rtsne(as.matrix(selected_tid_tpms), perplexity = 7, check_duplicates = FALSE) # Run TSNE
    #   plot(tsne_out$Y,col=class_colors,asp=1 , main = gene_to_plot)
    #   dev.off()
    # }
    # stripplot(tpm$ENST00000641190.1, groups = tpm$classL)
  }
  
  
}




