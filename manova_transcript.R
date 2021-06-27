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
  transposed_data <- log2(transposed_data+1)
  # print(transposed_data)
  pivoted_data <- transposed_data %>% 
    mutate(transcript_id = row.names(.)) %>% 
    tidyr::gather("Sample.ID", "TPM", 1:sample_size)
  pivoted_data <- pivoted_data %>% left_join(phenotype)
  # print(pivoted_data)
  ggplot(pivoted_data[which(pivoted_data$TPM>0),], aes(x=transcript_id, y=TPM, fill=phenotype)) + 
    geom_boxplot(position=position_dodge(1)) + 
    geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1), dotsize = 0.3, binwidth= 1/3) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    ggtitle(gene_id) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(y = "log(TPM+1)")
  ggsave(filename = paste0(gene_id, ".png"), path= paste0("experiment-outputs/", experiment_name, "/boxplots") ) 
}

if (TRUE) 
{
  experimentList <- c(
    "bc03ln_tumor_vs_non_tumor",
    "bc07ln_tumor_vs_non_tumor",
    "er_positive_tumor_vs_non_tumor",
    "her2_positive_tumor_vs_non_tumor"
  )
  
  
  gtf_df <- as.data.frame(rtracklayer::import('data/gencode.v27.annotation.gtf.gz'))
  expression_reads <- read.table("data/expression.matrix.tx.numreads.tsv", sep = '\t' , header = TRUE , stringsAsFactors = FALSE)
  rownames(expression_reads) <- expression_reads[,1] # add transcriptIDs as rownames
  transcript_types <- read.table("data/transcript_types.txt", sep = '\t', header= TRUE, stringsAsFactors = FALSE)
  
  for(i in 1:length(experimentList))
  {
    experiment_name <- experimentList[i]
    print(experiment_name)
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
    output_transcript_csv <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), 
                                      c("gene_id", "transcript_1", "transcript1_type", "transcript_2", "transcript2_type",
                                        "manova-p-value", "fold_change_difference"))
    
    nRowgenes  <- nrow(transcripts_genes)
    for(i in 1:nRowgenes) # if debug start at 415
    {
      # selected_row <- transcripts_genes %>% filter(gene_id == "ENSG00000115884.10")
      print(paste0(i, " - ", nRowgenes))
      selected_row <- transcripts_genes[i,]
      selected_gene_id<- selected_row$gene_id
      selected_tids<- selected_row$transcript_ids[[1]]
      selected_tid_tpms <- transpose_tpms %>% select(all_of(selected_tids))
      
      #DTU FILTERING
      tpms_w_classlabels <- selected_tid_tpms %>% mutate(cls_lbl = classLabels)
      
      clusterCenters <- tpms_w_classlabels %>% group_by(cls_lbl) %>% summarise_all(list(median))
      
      clusterCenters <- rbind ( clusterCenters, data.frame(cls_lbl = "ratio" ,  clusterCenters[1,-1] / clusterCenters[2,-1] ))
      
      # ratioList <- as.numeric( clusterCenters[1,-1] / clusterCenters[2,-1] )
      
      # ratioList = ratioList[!is.nan(ratioList)]  #1  replace NANS
      
      # clusterCenters <- clusterCenters %>% filter(is.na())
      
      # isoformSwitchFound = FALSE
      # if(any(cl>1 ) && any(ratioList>0 && ratioList<1))
      #   isoformSwitchFound  = TRUE
      # clusterCenters <- clusterCenters %>% mutate(totalSum = rowSums(across(where(is.numeric))))
      clusterCenters <- t_dt(column_to_rownames(clusterCenters, var = "cls_lbl"))
      clusterCenters <- clusterCenters %>% filter(`Non-Tumor`>10 & Tumor>10  & !is.nan(ratio))
      significant_transcripts_found <- (nrow(clusterCenters) >=1)
      isoformSwitchFound = FALSE
      ratioList = as.vector(clusterCenters[['ratio']])
      if(any(ratioList>1 ) && any(ratioList>0 && ratioList<1))
        isoformSwitchFound  = TRUE
      # print(paste0("Selected gene -> ", selected_gene_id))
      if(significant_transcripts_found && isoformSwitchFound)
      {
        print(paste0("Significant switch found for gene -> ", selected_gene_id))
        transcript_size_gene <- length(selected_tids)
        # Perform manova on it
        cmp <- cmpoutput("scRNA-seq", transcript_size_gene , selected_tid_tpms, as.factor(classLabels))
        # print( paste0(i, " - ", nRowgenes))
        # rbindlist(list(output,c(selected_gene_id, cmp$p.values$manova) ))
        # output[i,] = c(selected_gene_id, cmp$p.values$manova)
        output <- rbind( output, data.frame(selected_gene_id, cmp$p.values$manova))
        print(paste0(cmp$p.values$manova))
        
        if(!is.na(cmp$p.values$manova) && cmp$p.values$manova < 0.05)
        {
          # 2 transcriptSelection and foldchangeDifference
          significant_transcript_ids <- row.names(clusterCenters)
          significant_transcript_ids_subsets <- combn(significant_transcript_ids,2)
          for(t_subset_id in 1:ncol(significant_transcript_ids_subsets))
          {
            selected_two_transcript <- significant_transcript_ids_subsets[,t_subset_id]
            tmp_clusterCenters <- t_dt(clusterCenters) %>% select(all_of(selected_two_transcript)) # selected two transcript only
            two_isoform_ratio_list <- t_dt(tmp_clusterCenters)[['ratio']]
            fold_change_difference <- diff(log2(t_dt(tmp_clusterCenters)[['ratio']]))
            if(any(two_isoform_ratio_list>1 ) && any(two_isoform_ratio_list>0 && two_isoform_ratio_list<1))
            {
              boxplot_for_transcripts(paste0(selected_gene_id, "_", t_subset_id),selected_tid_tpms %>% select(all_of(selected_two_transcript)),phenotypes,sampleSize)
              # write.csv( selected_tid_tpms %>% mutate(cls_lbl = classLabels) , paste0("experiment-outputs/", experiment_name,"/gene-tpms/" , selected_gene_id,".csv"))
              # print(clusterCenters)
              
              Transcript1_Type = (transcript_types %>% filter(Transcript.ID == selected_two_transcript[1]))[["Type"]]
              Transcript2_Type = (transcript_types %>% filter(Transcript.ID == selected_two_transcript[2]))[["Type"]]
              
              output_transcript_csv <- rbind(output_transcript_csv, 
                                             data.frame(selected_gene_id, 
                                                        selected_two_transcript[1],
                                                        Transcript1_Type,
                                                        selected_two_transcript[2],
                                                        Transcript2_Type,
                                                        cmp$p.values$manova,
                                                        fold_change_difference))
              # print(paste0("Gene_id->", selected_gene_id))
              # print(paste0("Selected transcript->", selected_two_transcript))
              # print(paste0("Fold change difference->", fold_change_difference))
              # print(paste0("Manova p-value->", cmp$p.values$manova))
              # print(clusterCenters)
            }
          }

          # fold change end
          write.csv( selected_tid_tpms %>% mutate(cls_lbl = classLabels) , paste0("experiment-outputs/", experiment_name,"/gene-tpms/" , selected_gene_id,".csv"))
        }
        
      }
      
    }
    
    colnames(output) <- c("gene_id", "manova_p_value")
    colnames(output_transcript_csv) <- c("gene_id", "transcript_1", "transcript1_type", "transcript_2", "transcript2_type",
                                         "manova_p_value", "fold_change_difference")
    
    output[,2] <- sapply(output[,2] , as.numeric) # convert p values to numeric
    
    output <- output %>% arrange(manova_p_value)
    
    write.csv( output , paste0("experiment-outputs/", experiment_name, "/manova_outputs.csv"))
    write.csv( output_transcript_csv %>% arrange(manova_p_value) , paste0("experiment-outputs/", experiment_name, "/transcript_fold_change.csv"))
    
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




