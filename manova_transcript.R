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
  # transposed_data <- log2(transposed_data+1)
  # print(transposed_data)
  colnames(transposed_data) <- phenotype$Sample.ID
  pivoted_data <- transposed_data %>% 
    mutate(transcript_id = row.names(.)) %>% 
    tidyr::gather("Sample.ID", "TPM", 1:sample_size)
  pivoted_data <- pivoted_data %>% left_join(phenotype)
  # print(pivoted_data)
  ggplot(pivoted_data, aes(x=transcript_id, y=TPM, fill=phenotype)) + 
    geom_boxplot(position=position_dodge(1)) + 
    geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1), dotsize = 0.3, binwidth= 1/3) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    ggtitle(gene_id) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(y = "log(TPM+1)")
  ggsave(filename = paste0(gene_id, ".jpeg"), path= paste0("experiment-outputs/", experiment_name, "/boxplots") ) 
}

violing_plot_for_transcripts = function(gene_id, tid_data, phenotype)
{
  
  with_groups <- left_join(tid_data %>% mutate(Sample.ID=phenotype$Sample.ID), phenotype)
  data <- with_groups %>%  pivot_longer(!c(phenotype,Sample.ID), names_to = "transcripts", values_to = "TPM")
  ggstatsplot::grouped_ggbetweenstats(data, 
                                      x = phenotype,
                                      y = TPM,
                                      grouping.var = transcripts,
                                      ylab = "TPM (log scale)",
                                      pairwise.display = "significant",
                                      p.adjust.method = "fdr",
                                      results.subtitle = FALSE)
  ggsave(filename = paste0(gene_id, ".jpeg"), path= paste0("experiment-outputs/", experiment_name, "/boxplots") ) 
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  x_median <- median(x)
  x_mean <- mean(x)
  x[x == 0] <- median(x)
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y_median <- median(y ,na.rm = TRUE)
  y_mean <- mean(y ,na.rm = TRUE)
  # y[is.na(y)] <- x_median
  # y[is.na(y)] <- x_mean
  y[is.na(y)] <- y_median
  # y[is.na(y)] <- y_mean
  y
}

if (TRUE) 
{
  experimentList <- c(
    "bc07_vs_bc07ln_tumor"
    # "bc03_vs_bc03ln_tumor"
    #25	ENSG00000143549.19	ENST00000611659.4	Intron_9	ENST00000341485.9	3UTR	0.0277421482451465	1.24215773503332
    # "er_tumor_vs_her2_tumor",
    # "er_tumor_vs_tnbc_tumor",
    # "er_tumor_vs_er_her2_tumor",
    # "her2_tumor_vs_tnbc_tumor",
    # "her2_tumor_vs_er_her2_tumor",
    # "tnbc_tumor_vs_er_her2_tumor"
  )
  
  setwd("/Users/vlmedia/Thesis/R-Projects/isoform-quant")
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
    
    output <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("gene_id", "manova_p_value", "bartlett_t_test_p_value"))
    output_transcript_csv <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), 
                                      c("gene_id", "transcript_1", "transcript1_type", "transcript_2", "transcript2_type",
                                        "manova-p-value", "fold_change_difference"))
    
    phenotypes <- read.csv(paste0("data/", experiment_name, ".csv"), sep = ",", header= TRUE)
    classLabels <- phenotypes$phenotype
    selected_samples <- phenotypes$Sample.ID
    
    # get only samples given in phenotypeCSV
    expression_reads_filtered <- expression_reads %>% select(Gene_or_Transcript_ID,all_of(selected_samples))
    # take out transcripts with zero read count on all samples
    # expression_reads_filtered <- expression_reads_filtered[rowSums(expression_reads_filtered[,-1])>0,] 
    
    sampleSize <- ncol(expression_reads_filtered) -1 # take out gene_id column
    
    expression_reads_filtered <- expression_reads_filtered %>% 
      mutate(SignificantCounts =  rowSums(.[,-1] > 10) ) %>% 
      filter(SignificantCounts > sampleSize * 0.9) %>% 
      select(-SignificantCounts)
    
    
    expression_reads_filtered_grouped <- t_dt(expression_reads_filtered %>% select(-Gene_or_Transcript_ID)) %>%  mutate(cls_lbl = classLabels)
    clusterCenters <- expression_reads_filtered_grouped %>% group_by(cls_lbl) %>% summarise_all(list(median))
    clusterCenters_ratio <- rbind ( clusterCenters, data.frame(cls_lbl = "ratio" ,  clusterCenters[1,-1] / clusterCenters[2,-1] ))
    clusterCenters_ratio <- t_dt(clusterCenters_ratio)
    colnames(clusterCenters_ratio) <- clusterCenters_ratio[1,]
    clusterCenters_ratio <- clusterCenters_ratio[-1,]
    group_1 <- colnames(clusterCenters_ratio)[1]
    group_2 <- colnames(clusterCenters_ratio)[2]
    clusterCenters_ratio <- mutate_all(clusterCenters_ratio, function(x) as.numeric(as.character(x)))
    
    clusterCenters_filtered <- clusterCenters_ratio %>% filter((!!as.name(group_1)>10 | !!as.name(group_2)>10) & !is.na(ratio)) # & !is.infinite(ratio) & ratio != 0)
    clusterCenters_filtered[clusterCenters_filtered == 0] <-1
    clusterCenters_filtered <- clusterCenters_filtered %>% mutate(ratio = !!as.name(group_1) / !!as.name(group_2))
    
    # clusterCenters_filtered2 <- clusterCenters_filtered %>% filter(ratio != "Inf" & ratio != 0)
    
    isoform_switch_genes <- clusterCenters_filtered %>% 
      rownames_to_column("Gene_or_Transcript_ID") %>% 
      left_join(gtf_df %>% select(gene_id, transcript_id ) %>% distinct(gene_id,transcript_id) , by = c("Gene_or_Transcript_ID" = "transcript_id")) %>% 
      group_by(gene_id) %>% 
      mutate(isoformSwitch = ifelse(min(ratio)<1 && max(ratio)>1, TRUE,FALSE)) %>% 
      filter(isoformSwitch == TRUE) %>% 
      filter(n()>1) %>% 
      summarise(transcript_ids = list(Gene_or_Transcript_ID)) %>%  # add list of transcript ids for each gene
      filter(!is.na(gene_id)) # filter out the na group
    
    col_selectable_clusterCenteres <- t_dt(clusterCenters_filtered)
    
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
    
    transcripts_genes_filtered <- transcripts_genes %>% filter(gene_id %in% isoform_switch_genes$gene_id)
    
    # transpose matrix and keep transcriptIDs as colnames and add seperate sampleID column
    expression_reads_filtered <- expression_reads_filtered[,-1] # Take out geneID column
    transpose_tpms <- t_dt(expression_reads_filtered)

    nRowgenes  <- nrow(isoform_switch_genes)
    for(i in 1:nRowgenes) # if debug start at 415
    {
      # selected_row <- isoform_switch_genes %>% filter(gene_id == "ENSG00000075624.13")
      print(paste0(i, " - ", nRowgenes))
      selected_row <- isoform_switch_genes[i,]
      selected_gene_id<- selected_row$gene_id
      selected_tids<- selected_row$transcript_ids[[1]]
      selected_tid_tpms <- transpose_tpms %>% select(all_of(selected_tids))
      
      isoform_switch_selected_tids <- (isoform_switch_genes %>% filter(gene_id == !!selected_gene_id ))$transcript_ids[[1]]
      
      selected_clusterCenter <- t_dt(col_selectable_clusterCenteres %>% select(all_of(isoform_switch_selected_tids)))
      
      #DTU FILTERING
      # tpms_w_classlabels <- selected_tid_tpms %>% mutate(cls_lbl = classLabels)
      
      # selected_clusterCenter <- tpms_w_classlabels %>% group_by(cls_lbl) %>% summarise_all(list(median))
      
      # selected_clusterCenter <- rbind ( selected_clusterCenter, data.frame(cls_lbl = "ratio" ,  clusterCenters[1,-1] / clusterCenters[2,-1] ))
      
      
      
      # ratioList <- as.numeric( clusterCenters[1,-1] / clusterCenters[2,-1] )
      
      # ratioList = ratioList[!is.nan(ratioList)]  #1  replace NANS
      
      # clusterCenters <- clusterCenters %>% filter(is.na())
      
      # isoformSwitchFound = FALSE
      # if(any(cl>1 ) && any(ratioList>0 && ratioList<1))
      #   isoformSwitchFound  = TRUE
      # clusterCenters <- clusterCenters %>% mutate(totalSum = rowSums(across(where(is.numeric))))
      # selected_clusterCenter <- t_dt(column_to_rownames(selected_clusterCenter, var = "cls_lbl"))
      # group_1 <- colnames(clusterCenters)[1]
      # group_2 <- colnames(clusterCenters)[2]
      # clusterCenters <- clusterCenters %>% filter(!!group_1>10 & !!group_2>10  & !is.nan(ratio) & !is.infinite(ratio))
      # significant_transcripts_found <- (nrow(clusterCenters) >=1)
      # isoformSwitchFound = FALSE
      # ratioList = as.vector(clusterCenters[['ratio']])
      # if(any(ratioList>1 ) && any(ratioList>0 && ratioList<1))
      #   isoformSwitchFound  = TRUE
      # # print(paste0("Selected gene -> ", selected_gene_id))
      # if(significant_transcripts_found && isoformSwitchFound)
      # {
        # print(paste0("Significant switch found for gene -> ", selected_gene_id))
      tid_tpms_groups <- left_join(log2(selected_tid_tpms+1) %>% mutate(Sample.ID=rownames(selected_tid_tpms)), phenotypes, by ="Sample.ID")

      outlier_corrected_tpms <- tid_tpms_groups %>% group_by(phenotype) %>%
        mutate_if(is.numeric, remove_outliers) %>% ungroup()  %>% select(where(is.numeric))
      # rownames(outlier_corrected_tpms) <- phenotypes$Sample.ID
      
      # tid_corrected_groups <- left_join(outlier_corrected_tpms %>% mutate(Sample.ID=phenotypes$Sample.ID), phenotypes)
      # ggbetweenstats(tid_corrected_groups,phenotype, ENST00000464611.1, outlier.tagging = TRUE)
      # ggbetweenstats(tid_corrected_groups,phenotype, ENST00000331789.9, outlier.tagging = TRUE)
      
        transcript_size_gene <- length(selected_tids)
        # Perform manova on it
        # cmp_non_corrected <- cmpoutput("scRNA-seq", transcript_size_gene , selected_tid_tpms, as.factor(classLabels))
        cmp <- cmpoutput("scRNA-seq", transcript_size_gene , outlier_corrected_tpms, as.factor(classLabels))
        t_test_assumptions <- assumptions(cmp)
        # print( paste0(i, " - ", nRowgenes))
        # rbindlist(list(output,c(selected_gene_id, cmp$p.values$manova) ))
        # output[i,] = c(selected_gene_id, cmp$p.values$manova)
        output <- rbind( output, data.frame(selected_gene_id, cmp$p.values$manova , t_test_assumptions$ttest$vartest[[1]]$p.value ))
        print(paste0(cmp$p.values$manova))
        
        if(!is.na(cmp$p.values$manova) && cmp$p.values$manova < 0.05)
        {
          # 2 transcriptSelection and foldchangeDifference
          significant_transcript_ids <- row.names(selected_clusterCenter)
          significant_transcript_ids_subsets <- combn(significant_transcript_ids,2)
          for(t_subset_id in 1:ncol(significant_transcript_ids_subsets))
          {
            selected_two_transcript <- significant_transcript_ids_subsets[,t_subset_id]
            tmp_clusterCenters <- t_dt(selected_clusterCenter) %>% select(all_of(selected_two_transcript)) # selected two transcript only
            two_isoform_ratio_list <- t_dt(tmp_clusterCenters)[['ratio']]
            fold_change_difference <- abs(diff(log2(as.numeric(t_dt(tmp_clusterCenters)[['ratio']]))))
            if(any(two_isoform_ratio_list>=1 ) && any(two_isoform_ratio_list<1))
            {
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
              # boxplot_for_transcripts(paste0(selected_gene_id,  "_", Transcript1_Type, "_",Transcript2_Type , "_", t_subset_id ),outlier_corrected_tpms %>% select(all_of(selected_two_transcript)),phenotypes,sampleSize)
              violing_plot_for_transcripts(paste0(selected_gene_id,  "_", Transcript1_Type, "_",Transcript2_Type , "_", t_subset_id ),outlier_corrected_tpms %>% select(all_of(selected_two_transcript)),phenotypes)
              
              # print(paste0("Gene_id->", selected_gene_id))
              # print(paste0("Selected transcript->", selected_two_transcript))
              # print(paste0("Fold change difference->", fold_change_difference))
              # print(paste0("Manova p-value->", cmp$p.values$manova))
              # print(selected_clusterCenter)
            }
          }

          # fold change end
          write.csv( selected_tid_tpms %>% mutate(cls_lbl = classLabels) , paste0("experiment-outputs/", experiment_name,"/gene-tpms/" , selected_gene_id,".csv"))
        }
        
      }
      
    
    colnames(output) <- c("gene_id", "manova_p_value", "bartlett_t_test_p_value")
    colnames(output_transcript_csv) <- c("gene_id", "transcript_1", "transcript1_type", "transcript_2", "transcript2_type",
                                         "manova_p_value", "fold_change_difference")
    
    output[,2] <- sapply(output[,2] , as.numeric) # convert p values to numeric
    
    output <- output %>% arrange(manova_p_value)
    
    write.csv( output , paste0("experiment-outputs/", experiment_name, "/manova_outputs.csv"))
    write.csv( output %>% filter(bartlett_t_test_p_value > 0.05) , paste0("experiment-outputs/", experiment_name, "/manova_outputs_t_test_passed.csv"))
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

# tid_tpms_groups <- left_join(selected_tid_tpms %>% mutate(Sample.ID=rownames(selected_tid_tpms)), phenotypes)
# 
# corrected <- tid_tpms_groups %>% group_by(phenotype) %>% 
#   mutate_if(is.numeric, remove_outliers) %>% ungroup()  %>% select(where(is.numeric))
# 
# cmp <- cmpoutput("scRNA-seq", transcript_size_gene , corrected, as.factor(classLabels))
# 
# shapiro.test(corrected$ENST00000509541.5)
# shapiro.test(corrected$ENST00000515580.1)
