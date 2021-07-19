library(dplyr)
library(ggvenn)  
# if (!require(devtools)) install.packages("devtools")
#   devtools::install_github("yanlinlin82/ggvenn")
library("ggVennDiagram")
# if (!require(devtools)) install.packages("devtools")
#   devtools::install_github("gaospecial/ggVennDiagram")
library("VennDiagram")

# experimentList <- c(
#   "er_tumor_vs_her2_tumor",
#   "er_tumor_vs_tnbc_tumor",
#   "er_tumor_vs_er_her2_tumor",
#   "her2_tumor_vs_tnbc_tumor",
#   "her2_tumor_vs_er_her2_tumor",
#   "tnbc_tumor_vs_er_her2_tumor"
# )
er_experiments <- c("er_tumor_vs_her2_tumor",
                    "er_tumor_vs_tnbc_tumor",
                    "er_tumor_vs_er_her2_tumor")

her2_experiments <- c("er_tumor_vs_her2_tumor",
                      "her2_tumor_vs_tnbc_tumor",
                      "her2_tumor_vs_er_her2_tumor")

er_her2_experiments <- c("er_tumor_vs_er_her2_tumor",
                         "her2_tumor_vs_er_her2_tumor",
                         "tnbc_tumor_vs_er_her2_tumor")

tnbc_experiments <- c("er_tumor_vs_tnbc_tumor",
                      "her2_tumor_vs_tnbc_tumor",
                      "tnbc_tumor_vs_er_her2_tumor")

all_exp_groups <- list(er_experiments, her2_experiments, er_her2_experiments, tnbc_experiments )
all_exp_names <- c("er", "her2", "er_her2", "tnbc" )


if(TRUE)
{
  for(k in 1:length(all_exp_groups))
  {
    experimentList <- all_exp_groups[[k]]
    manova_output_holder <- list()
    for(i in 1:length(experimentList))
    {
      experiment_name <- experimentList[i]
      exp_manova_outputs <- read.csv(paste0("experiment-outputs/", experiment_name, "/manova_outputs.csv"), sep = ",", header= TRUE)
      manova_output_holder[[as.name(experiment_name)]] <- unique(exp_manova_outputs[["gene_id"]])
    }

    ggvenn(
      manova_output_holder,
      fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
      stroke_size = 0.1, set_name_size = 3
    )
    ggsave(file=paste0(all_exp_names[k],"_gene.jpeg"))
    # ggVennDiagram(manova_output_holder, label_alpha = 0)
  }
  all_common <- list()
  for(k in 1:length(all_exp_groups))
  {
    experimentList <- all_exp_groups[[k]]
    isoform_switch_holder <- list()
    for(i in 1:length(experimentList))
    {
      experiment_name <- experimentList[i]
      exp_isoform_switches <- read.csv(paste0("experiment-outputs/", experiment_name, "/transcript_fold_change.csv"), sep = ",", header= TRUE)
      tmp_isoform_switch_list <- c()
      for(j in 1:nrow(exp_isoform_switches))
      {
        tx_to_sort <- c(exp_isoform_switches[j,][["transcript_1"]], exp_isoform_switches[j,][["transcript_2"]])
        tx_sorted <- str_sort(tx_to_sort)
        tmp_isoform_switch_list <- append(tmp_isoform_switch_list, paste0(exp_isoform_switches[j,]["gene_id"],"-",
                                                                          tx_sorted[1],"-",
                                                                          tx_sorted[2]))
      }
      isoform_switch_holder[[as.name(experiment_name)]] <- tmp_isoform_switch_list
    }
    ol <- calculate.overlap(x =isoform_switch_holder )
    all_common[[all_exp_names[k]]]<- ol$a5
    sink(paste0(all_exp_names[k], "_common_isoform_switches.txt"))
    writeLines(unlist(lapply(ol$a5, paste, collapse=" ")))
    sink()
    # dput(ol$a5 , paste0(all_exp_names[k], "_common_isoform_switches.txt"))
    ol_size=sapply(ol, length)
    ggvenn(
      isoform_switch_holder,
      fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
      stroke_size = 0.1, set_name_size = 3
    )
    ggsave(file=paste0(all_exp_names[k],"_isoformswitch.jpeg"))
    # ggVennDiagram(manova_output_holder, label_alpha = 0)
  }
  ggvenn(
    all_common,
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 4
  )
  ggsave(file=paste0("all_common_isoform_switch.jpeg"))
  
}
