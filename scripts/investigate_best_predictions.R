#investigate_best_predictions.r


# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(coin, quietly = T))
suppressPackageStartupMessages (library(see, quietly = T))
suppressPackageStartupMessages (library(patchwork, quietly = T))
suppressPackageStartupMessages (library(svglite, quietly = T))
# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="DawnRank;OncoImpact;PNC;PRODIGY;PersonaDrive;SCS;sysSVM2;PhenoDriverR;CSN_NCUA", 
              help="algorithms to include in comparison separated by semicolons, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Threads")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


network_choice <- opt$network
algorithms <- opt$algorithms
algorithms <- str_split(algorithms,";") %>% unlist()
threads <- opt$threads

# SL-partner score filtering threshold (default 0.5)
SL_th <- 0.5
# pvalue filtering threshold for mann-whitney test for rare essential genes
pval <- 0.05
# max number of SL partners to consider for each gene
max_SL <- 5


if(threads>1){
  cl <- makeCluster(threads, outfile = "log/investigate_best_predictions.log")
  registerDoParallel(cl)
}



#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv"))
samples <- sample_info$cell_ID %>% unique() %>% sort()

#############################
# Read In Driver Predictions
#############################

source("scripts/read_driver_results.R")

driver_gene_predictions <- read_driver_results("benchmark",algorithms,"ALL")

##############################
# Quantify Driver Predictions
##############################

count_driver_predictions_per_tissue <-  driver_gene_predictions %>%
  #count cohort sizes
  group_by(cancer_type) %>%
  mutate(cohort_size=length(unique(sample_ID))
  ) %>%
  ungroup() %>%
  #count driver occurrences across algorithms
  group_by(cancer_type,driver) %>%
  mutate(driver_occurrences=length(unique(sample_ID))) %>%
  ungroup() %>%
  #get mean rank per algorithm
  group_by(algorithm,cancer_type,driver) %>%
  mutate(mean_rank=round(mean(rank,na.rm = T),2)) %>%
  ungroup() %>%
  dplyr::select(algorithm,cancer_type,cohort_size,driver,driver_occurrences,mean_rank) %>%
  unique() %>%
  pivot_wider(names_from = "algorithm", values_from = "mean_rank", names_prefix = "mean_rank_")

count_driver_predictions <-  driver_gene_predictions %>%
  mutate(cohort_size=length(unique(sample_ID))) %>%
  #count driver occurrences across algorithms
  group_by(driver) %>%
  mutate(driver_occurrences=length(unique(sample_ID))) %>%
  ungroup() %>%
  #get mean rank per algorithm
  group_by(algorithm,driver) %>%
  mutate(mean_rank=round(mean(rank,na.rm = T),2)) %>%
  ungroup() %>%
  dplyr::select(algorithm,driver,cohort_size,driver_occurrences,mean_rank) %>%
  unique() %>%
  pivot_wider(names_from = "algorithm", values_from = "mean_rank", names_prefix = "mean_rank_")


#####################################
# Read In Predicted Essential Genes #
#####################################

message(paste0("Reading in rare essential gene predictions"))
rare_essential_genes <- fread(paste0("results/benchmark/network_",network_choice,"/rare_essential_genes.csv")) %>%
  filter(sample_ID %in% samples, algorithm%in% unique(driver_gene_predictions$algorithm))

all_essential_genes <- fread(paste0("results/benchmark/network_",network_choice,"/all_essential_genes.csv")) %>%
  filter(sample_ID %in% samples, algorithm%in% unique(driver_gene_predictions$algorithm))

gene_effect_z_scores <- fread("benchmark_data/gene_effect_z_scores.csv") %>%
  dplyr::select(sample_ID=cell_ID, gene_ID, UIg=weighted_average)

count_essential_gene_predictions_per_tissue <- all_essential_genes %>%
  dplyr::select(cancer_type,sample_ID,driver,annotation,target) %>%
  unique() %>%
  # combine predictions with ground truth data
  left_join(gene_effect_z_scores, by = c("sample_ID","target"="gene_ID")) %>%
  
  #summarise driver info
  group_by(cancer_type,driver,target) %>%
  mutate(target_occurrences=length(unique(sample_ID)),
         mean_UIg=mean(UIg, na.rm = T),
         gene_effect_ground_truth_n = length(which(!is.na(UIg)))
         ) %>%
  dplyr::select(-c(sample_ID,UIg)) %>% unique()

count_essential_gene_predictions <- all_essential_genes %>%
  dplyr::select(sample_ID,driver,annotation,target) %>%
  unique() %>%
  # combine predictions with ground truth data
  left_join(gene_effect_z_scores, by = c("sample_ID","target"="gene_ID")) %>%
  #summarise driver info
  group_by(driver,target) %>%
  mutate(target_occurrences=length(unique(sample_ID)),
         mean_UIg=mean(UIg, na.rm = T),
         gene_effect_ground_truth_n = length(which(!is.na(UIg)))
         ) %>%
  dplyr::select(-c(sample_ID,UIg)) %>% unique()
  

###########################
# Read In Predicted Drugs #
###########################

message(paste0("Reading in rare drug predictions"))
rare_drug_predictions <- fread(paste0("results/benchmark/network_",network_choice,"/rare_drug_predictions.csv"))%>%
  filter(sample_ID %in% samples, algorithm%in% unique(driver_gene_predictions$algorithm))

all_drug_predictions <- fread(paste0("results/benchmark/network_",network_choice,"/all_drug_predictions.csv"))%>%
  filter(sample_ID %in% samples, algorithm%in% unique(driver_gene_predictions$algorithm))

drug_sensitivity <- fread("benchmark_data/drug_sensitivity.csv") %>%
  dplyr::select(sample_ID=cell_ID, drug_ID, UId=weighted_average)

count_drug_predictions_per_tissue <- all_drug_predictions %>%
  dplyr::select(cancer_type,sample_ID,driver,target,drug_ID) %>%
  unique() %>%
  # combine predictions with ground truth data
  left_join(drug_sensitivity, by = c("sample_ID","drug_ID"="drug_ID")) %>%
  #summarise driver info
  group_by(cancer_type,driver,target,drug_ID) %>%
  mutate(drug_occurrences=length(unique(sample_ID)),
         mean_UId=mean(UId, na.rm = T),
         drug_ground_truth_n = length(which(!is.na(UId)))
         ) %>%
  dplyr::select(-c(sample_ID,UId)) %>% unique()


count_drug_predictions <- all_drug_predictions %>%
  dplyr::select(sample_ID,driver,target,drug_ID) %>%
  unique() %>%
  # combine predictions with ground truth data
  left_join(drug_sensitivity, by = c("sample_ID","drug_ID"="drug_ID")) %>%
  #summarise driver info
  group_by(target,driver,drug_ID) %>%
  mutate(drug_occurrences=length(unique(sample_ID)),
         mean_UId=mean(UId, na.rm = T),
         drug_ground_truth_n = length(which(!is.na(UId)))
         ) %>%
  dplyr::select(-c(sample_ID,UId)) %>% unique()

# Combine all


combined_result_per_tissue <- count_driver_predictions_per_tissue %>%
  left_join(count_essential_gene_predictions_per_tissue, by= c("cancer_type","driver")) %>%
  left_join(count_drug_predictions_per_tissue, by= c("cancer_type","driver","target")) %>%
  ungroup() %>%
  # calculate frequencies
  mutate(driver_frequency=driver_occurrences/cohort_size,
         target_gene_frequency=target_occurrences/cohort_size,
         drug_frequency=drug_occurrences/cohort_size
         ) %>%
  #filters
  filter(
    # driver occurs in >1% of population
    driver_frequency > 0.01,
    # essential gene predicted in > 1% of population
    #target_gene_frequency > 0.01,
    # drug predicted in > 1% of population
    #drug_frequency > 0.01,
    # uniqueness <-1
    mean_UIg < -1,
    mean_UId < -0.5,
    gene_effect_ground_truth_n > 2
  ) %>%
  # At least one algorithm has ranked the driver < 10
  #pivot_wider(names_from = algorithm, values_from = mean_rank) %>%
  rowwise() %>%
  filter(min(mean_rank_CSN_NCUA,
             mean_rank_DawnRank,
             mean_rank_OncoImpact,
             mean_rank_PNC,
             mean_rank_PRODIGY,
             mean_rank_PersonaDrive,
             mean_rank_PhenoDriverR,
             mean_rank_SCS,
             mean_rank_sysSVM2,
             na.rm = T) <=10)
  
  

write_csv(combined_result_per_tissue %>%
            dplyr::select(
              cancer_type,driver,annotation,driver_frequency,target,drug_ID,mean_UIg,mean_UId, starts_with("mean_rank")
            ) %>%
            arrange(mean_UIg)
            ,"best_predictions_per_tissue.csv")



combined_result <- count_driver_predictions %>%
  left_join(count_essential_gene_predictions, by= c("driver")) %>%
  left_join(count_drug_predictions, by= c("driver","target")) %>%
  ungroup() %>%
  # calculate frequencies
  mutate(driver_frequency=driver_occurrences/cohort_size,
         target_gene_frequency=target_occurrences/cohort_size,
         drug_frequency=drug_occurrences/cohort_size
  ) %>%
  #filters
  filter(
    # driver occurs in >1% of population
    driver_frequency > 0.01,
    # essential gene predicted in > 1% of population
    #target_gene_frequency > 0.01,
    # drug predicted in > 1% of population
    #drug_frequency > 0.01,
    # uniqueness <-1
    mean_UIg < -1,
    #mean_UId < -0.5,
    gene_effect_ground_truth_n > 5
  ) %>%
  # At least one algorithm has ranked the driver < 10
  rowwise() %>%
  filter(min(mean_rank_CSN_NCUA,
             mean_rank_DawnRank,
             mean_rank_OncoImpact,
             mean_rank_PNC,
             mean_rank_PRODIGY,
             mean_rank_PersonaDrive,
             mean_rank_PhenoDriverR,
             mean_rank_SCS,
             mean_rank_sysSVM2,
             na.rm = T) <=10)


write_csv(combined_result %>%
            dplyr::select(
              driver,annotation,driver_frequency,target,drug_ID,mean_UIg,mean_UId, starts_with("mean_rank")
            ) %>%
            arrange(mean_UIg)
            ,"best_predictions.csv")


# Plots

plot_indiv_predictions <- function(
    plot_cancer_type = "Lung",
    plot_driver = "CDKN1B",
    plot_target = "CASP8",
    plot_drug = "BRD-K00023644-001-01-9",
    max_rank = 10,
    show_x_axis=F,
    show_y_axis=F,
    show_guides=F,
    driver_type=ifelse(plot_driver==plot_target,"GOF","LOF")
){
  
  target_samples_gene <- fread(paste0("results/benchmark/network_",network_choice,"/all_essential_genes.csv")) %>%
    filter(driver==plot_driver, target==plot_target, driver_rank <= max_rank,cancer_type==plot_cancer_type) %>%
    pull(sample_ID) %>%
    unique()
  
  target_samples_drug <- fread(paste0("results/benchmark/network_",network_choice,"/all_drug_predictions.csv")) %>%
    filter(driver==plot_driver,driver_rank <= max_rank, target==plot_target, drug_ID==plot_drug,cancer_type==plot_cancer_type) %>%
    pull(sample_ID) %>%
    unique()
  
  target_samples <- union(target_samples_drug,target_samples_gene)
  
  same_cancer_type_samples <- sample_info %>%
    filter(lineage==plot_cancer_type, !cell_ID%in%target_samples) %>%
    pull(cell_ID) %>%
    unique()
  
  other_samples <- sample_info %>%
    filter(!cell_ID%in%c(target_samples,same_cancer_type_samples)) %>%
    pull(cell_ID) %>%
    unique()
  
  plot_data <- full_join(
    fread("benchmark_data/gene_effect_z_scores.csv") %>%
      dplyr::select(sample_ID=cell_ID, gene_ID, gene_effect) %>%
      filter(gene_ID==plot_target),
    
    fread("benchmark_data/drug_sensitivity.csv") %>%
      dplyr::select(sample_ID=cell_ID, drug_ID, sensitivity_value, sensitivity_source) %>%
      filter(drug_ID==plot_drug), 
    by = c("sample_ID")
  ) %>%
    mutate(group= ifelse(sample_ID %in% target_samples, "target",
                         ifelse(sample_ID%in%same_cancer_type_samples, "same_cell_type", "other")))
  
  pos <- position_jitter(width = 0.1, seed = 2)
  
  
  p1 <- ggplot(plot_data %>% filter(!is.na(gene_ID)), aes(x=gene_ID,y=gene_effect)) +
    geom_violinhalf(data=plot_data %>% filter(!is.na(gene_ID),group=="other"),aes(fill=group), alpha = 0.75) +
    geom_violinhalf(data=plot_data %>% filter(!is.na(gene_ID),group=="same_cell_type"),aes(fill=group), flip = T, alpha = 0.75) +
    geom_point(data = plot_data %>% filter(!is.na(gene_ID),group=="target"), aes(fill=group), colour="black", shape=21, alpha=1, size=2, position = pos, show.legend = T) +
    
    #scale_x_discrete(limits=rev) +
    scale_fill_manual(values=c("#8CB2E5","#61B550","red"), breaks = c("other","same_cell_type","target"), labels = c("All Cells", "Same Cell Type","Predicted Target Cells"), name = "") +
    ylab("") +
    xlab(paste0(plot_cancer_type,"\n",driver_type,"-",plot_driver)) +
    #geom_text(data = plot_data %>% filter(group==loi,cell_ID== coi), label=coi, colour = "red", angle=90, hjust = 1, nudge_x = 0.6) +
    #geom_text_repel(data = plot_data %>% filter(group==loi, cell_ID!=coi), aes(label=label), box.padding = unit(1,"cm"), position = pos, max.overlaps = Inf) +
    coord_flip() +
    ggtitle(plot_target) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 8), axis.text.y = element_blank(), axis.title.y = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
    guides(fill="none", colour="none")
  
  
  p2 <- ggplot(plot_data %>% filter(!is.na(drug_ID)), aes(x=drug_ID,y=sensitivity_value)) +
    geom_violinhalf(data=plot_data %>% filter(!is.na(drug_ID),group=="other"),aes(fill=group), alpha = 0.75) +
    geom_violinhalf(data=plot_data %>% filter(!is.na(drug_ID),group=="same_cell_type"),aes(fill=group), flip = T, alpha = 0.75) +
    geom_point(data = plot_data %>% filter(!is.na(drug_ID),group=="target"), aes(fill=group), colour="black", shape=21, alpha=1, size=2, position = pos, show.legend = T) +
    
    #scale_x_discrete(limits=rev) +
    scale_fill_manual(values=c("#8CB2E5","#61B550","red"), breaks = c("other","same_cell_type","target"), labels = c("All Cells", "Same Cell Type","Predicted Target Cells"), name = "") +
    ylab("") +
    xlab("") +
    #geom_text(data = plot_data %>% filter(group==loi,cell_ID== coi), label=coi, colour = "red", angle=90, hjust = 1, nudge_x = 0.6) +
    #geom_text_repel(data = plot_data %>% filter(group==loi, cell_ID!=coi), aes(label=label), box.padding = unit(1,"cm"), position = pos, max.overlaps = Inf) +
    coord_flip() +
    ggtitle(plot_drug) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size = 8), axis.text.y = element_blank(), axis.title.y = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm" ))
  
  if(show_x_axis){
    p1 <- p1 + ylab("Gene Effect")
    p2 <- p2 + ylab("Drug Sensitivity")
  }
  
  if(!show_guides){
    #p1 <- p1 + guides(fill="none", colour="none")
    p2 <- p2 + guides(fill="none", colour="none")
  }
  
  
  
  p1 | p2 + plot_layout(nrow=1, guides = "collect", widths = 1) & theme(legend.position = "bottom")
}

p1.1 <- plot_indiv_predictions(
  plot_cancer_type = "Breast",
  plot_driver = "KRAS",
  plot_target = "KRAS",
  plot_drug = "BRD-K00003260-001-01-9",
  max_rank = Inf
)

p1.2 <- plot_indiv_predictions(
  plot_cancer_type = "Myeloid",
  plot_driver = "BCR",
  plot_target = "BCR",
  plot_drug = "Dasatinib",
  max_rank = Inf
)

p1.3 <- plot_indiv_predictions(
  plot_cancer_type = "Liver",
  plot_driver = "SMARCA4",
  plot_target = "SMARCA2",
  plot_drug = "PFI-3",
  max_rank = Inf
)

p1.4 <- plot_indiv_predictions(
  plot_cancer_type = "Biliary_Tract",
  plot_driver = "BRAF",
  plot_target = "BRAF",
  plot_drug = "BRD-K00003576-001-01-9",
  max_rank = Inf
)

p1.5 <- plot_indiv_predictions(
  plot_cancer_type = "Skin",
  plot_driver = "PPM1D",
  plot_target = "PPM1D",
  plot_drug = "GSK2830371A",
  max_rank = Inf
)

p1.6 <- plot_indiv_predictions(
  plot_cancer_type = "Myeloid",
  plot_driver = "ABL1",
  plot_target = "ABL1",
  plot_drug = "Nilotinib",
  max_rank = Inf
)

p1.7 <- plot_indiv_predictions(
  plot_cancer_type = "Skin",
  plot_driver = "BIRC2",
  plot_target = "BIRC2",
  plot_drug = "BRD-K83030136-001-02-9",
  max_rank = Inf
)

p1.8 <- plot_indiv_predictions(
  plot_cancer_type = "Esophagus_Stomach",
  plot_driver = "ERBB2",
  plot_target = "ERBB2",
  plot_drug = "Lapatinib",
  max_rank = 10,
  show_x_axis = T,
  show_guides = T
)

p1.1/p1.2/p1.3/p1.4/p1.5/p1.6/p1.7/p1.8 + plot_layout(axis_titles = "collect", widths = 1)

ggsave("plots/benchmark/best_predictions_tissue_specific.svg", device = svglite, width = 20, height = 20, units = "cm")









plot_indiv_predictions2 <- function(
    plot_driver = "CDKN1B",
    plot_target = "CASP8",
    plot_drug = "BRD-K00023644-001-01-9",
    max_rank = 10,
    show_x_axis=F,
    show_guides=F,
    driver_type=ifelse(plot_driver==plot_target,"GOF","LOF")
){
  
  target_samples_gene <- fread(paste0("results/benchmark/network_",network_choice,"/all_essential_genes.csv")) %>%
    filter(driver==plot_driver, target==plot_target, driver_rank <= max_rank) %>%
    pull(sample_ID) %>%
    unique()
  print(length(target_samples_gene))
  target_samples_drug <- fread(paste0("results/benchmark/network_",network_choice,"/all_drug_predictions.csv")) %>%
    filter(driver==plot_driver,driver_rank <= max_rank, target==plot_target, drug_ID==plot_drug) %>%
    pull(sample_ID) %>%
    unique()
  
  print(length(target_samples_drug))
  
  target_samples <- union(target_samples_drug,target_samples_gene)
  
  other_samples <- fread(paste0("benchmark_data/network_",network_choice,"/sample_info.csv")) %>%
    filter(!cell_ID%in%c(target_samples)) %>%
    pull(cell_ID) %>%
    unique()
  
  plot_data <- left_join(
    fread("benchmark_data/gene_effect_z_scores.csv") %>%
      dplyr::select(sample_ID=cell_ID, gene_ID, gene_effect) %>%
      filter(gene_ID==plot_target),
    
    fread("benchmark_data/drug_sensitivity.csv") %>%
      dplyr::select(sample_ID=cell_ID, drug_ID, sensitivity_value, sensitivity_source) %>%
      filter(drug_ID==plot_drug), 
    by = c("sample_ID")
  ) %>%
    mutate(group= ifelse(sample_ID %in% target_samples, "target", "other"))

  
  pos <- position_jitter(width = 0.1, seed = 2)
  
  pp1 <- ggplot(plot_data %>% filter(!is.na(gene_ID)), aes(x=gene_ID,y=gene_effect)) +
    geom_violinhalf(data=plot_data %>% filter(!is.na(gene_ID),group=="other"),aes(fill=group), alpha = 0.75) +
    geom_violinhalf(data=plot_data %>% filter(!is.na(gene_ID),group=="target"),aes(fill=group), flip = T, alpha = 0.75) +
    geom_point(data = plot_data %>% filter(!is.na(gene_ID),group=="target"), aes(fill=group), colour="black", shape=21, alpha=1, size=2, position = pos, show.legend = T) +
    
    #scale_x_discrete(limits=rev) +
    scale_fill_manual(values=c("#8CB2E5","red"), breaks = c("other","target"), labels = c("Other Cells","Predicted Essential"), name = "") +
    ylab("") +
    xlab(paste0(driver_type,"-",plot_driver)) +
    #geom_text(data = plot_data %>% filter(group==loi,cell_ID== coi), label=coi, colour = "red", angle=90, hjust = 1, nudge_x = 0.6) +
    #geom_text_repel(data = plot_data %>% filter(group==loi, cell_ID!=coi), aes(label=label), box.padding = unit(1,"cm"), position = pos, max.overlaps = Inf) +
    coord_flip() +
    ggtitle(plot_target) +
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5, size = 8), axis.text.y = element_blank(), axis.title.y = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) +
    guides(fill="none", colour="none")
  
  pp2 <- ggplot(plot_data %>% filter(!is.na(drug_ID)), aes(x=drug_ID,y=sensitivity_value)) +
    geom_violinhalf(data=plot_data %>% filter(!is.na(drug_ID),group=="other"),aes(fill=group), alpha = 0.75) +
    geom_violinhalf(data=plot_data %>% filter(!is.na(drug_ID),group=="target"),aes(fill=group), flip = T, alpha = 0.75) +
    geom_point(data = plot_data %>% filter(!is.na(drug_ID),group=="target"), aes(fill=group), colour="black", shape=21, alpha=1, size=2, position = pos, show.legend = T) +
    
    #scale_x_discrete(limits=rev) +
    scale_fill_manual(values=c("#8CB2E5","red"), breaks = c("other","target"), labels = c("All Cells","Predicted Target Cells"), name = "") +
    ylab("") +
    xlab("") +
    ggtitle(plot_drug) +
    #geom_text(data = plot_data %>% filter(group==loi,cell_ID== coi), label=coi, colour = "red", angle=90, hjust = 1, nudge_x = 0.6) +
    #geom_text_repel(data = plot_data %>% filter(group==loi, cell_ID!=coi), aes(label=label), box.padding = unit(1,"cm"), position = pos, max.overlaps = Inf) +
    coord_flip() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 8), axis.text.y = element_blank(), axis.title.y = element_text(size=8), plot.margin = unit(c(0,0,0,0), "cm")) 
    
  
  if(show_x_axis){
    pp1 <- pp1 + ylab("Gene Effect")
    pp2 <- pp2 + ylab("Drug Sensitivity")
  }
  
  if(!show_guides){
    #p1 <- p1 + guides(fill="none", colour="none")
    pp2 <- pp2 + guides(fill="none", colour="none")
  }
  
  
  
  pp1 | pp2 + plot_layout(nrow=1, guides = "collect", widths = 1) & theme(legend.position = "bottom")
  
}

p2.1 <- plot_indiv_predictions2(
  plot_driver = "HRAS",
  plot_target = "HRAS",
  plot_drug = "Kobe2602",
  max_rank = Inf
)

p2.2 <- plot_indiv_predictions2(
  plot_driver = "SIRT2",
  plot_target = "SIRT2",
  plot_drug = "Panobinostat",
  max_rank = Inf
)

p2.3 <- plot_indiv_predictions2(
  plot_driver = "MAPK1",
  plot_target = "MTOR",
  plot_drug = "OSI-027",
  max_rank = Inf,
  show_x_axis = T,
  show_guides = T
)


p2.1/p2.2/p2.3 + plot_layout(axis_titles = "collect", widths = 1)

ggsave("plots/benchmark/best_predictions.svg", device = svglite, width = 20, height = 8, units = "cm")


# which driver predictions frequently led to strong essential gene / drug predictions?

#columns
#cancer type
#driver gene
#frequency per algorithm
#average rank per algorithm
#essential gene
#drug
#average UIg
#average UId
