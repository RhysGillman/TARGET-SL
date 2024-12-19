#SL_pipeline_proof_of_concept.R



# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(ggpubr, quietly = T))
suppressPackageStartupMessages (library(svglite, quietly = T))
suppressPackageStartupMessages (library(effectsize, quietly = T))
suppressPackageStartupMessages (library(rstatix, quietly = T))
suppressPackageStartupMessages (library(patchwork, quietly = T))

# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="consensus", 
              help="algorithms to include in comparison separated by spaces, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Threads")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


network_choice <- opt$network
algorithms <- opt$algorithms
threads <- opt$threads

SL_th <- 0.5
max_SL <- 5
driver_algorithm_top_n <- 1
drug_source="GDSC2"

if(threads>1){
  #registerDoParallel(cores=threads)
  cl <- makeCluster(threads, outfile = "log/proof_of_concept.log")
  registerDoParallel(cl)
}

source("scripts/benchmark_functions.R")

###############
# Sample Info #
###############

sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv"), show_col_types = F)
samples <- sample_info$cell_ID %>% unique() %>% sort()


#############################
# Gene List
#############################

genes <- fread(paste0("benchmark_data/network_",network_choice,"/counts.csv"), select = "gene_ID") %>% pull(gene_ID)


#####################################
# Read In Predicted Essential Genes #
#####################################

message(paste0("Reading in essential gene predictions"))
all_essential_genes <- fread(paste0("results/benchmark/network_",network_choice,"/all_essential_genes.csv")) %>%
  # keep only the top N driver predictions
  filter(sample_ID %in% samples, toupper(algorithm)==toupper(algorithms)) %>%
  group_by(sample_ID) %>%
  filter(driver_rank%in%head(sort(unique(driver_rank)),driver_algorithm_top_n)) %>%
  dplyr::select(algorithm,sample_ID, gene_ID=target)
message(paste0("Reading in rare essential gene predictions"))
rare_essential_genes <- fread(paste0("results/benchmark/network_",network_choice,"/rare_essential_genes.csv")) %>%
  # keep only the top N driver predictions
  filter(sample_ID %in% samples, toupper(algorithm)==toupper(algorithms)) %>%
  group_by(sample_ID) %>%
  filter(driver_rank%in%head(sort(unique(driver_rank)),driver_algorithm_top_n)) %>%
  dplyr::select(algorithm,sample_ID, gene_ID=target)


###########################
# Read In Predicted Drugs #
###########################

# Drug-Gene Interactions

drug_targets <- fread(paste0("benchmark_data/network_",network_choice,"/drug_targets.csv")) %>%
  dplyr::select(drug_ID,gene_ID,source) %>%
  unique()


message(paste0("Reading in drug predictions"))
all_drug_predictions <- fread(paste0("results/benchmark/network_",network_choice,"/all_drug_predictions.csv")) %>%
  filter(sample_ID %in% samples, toupper(algorithm)==toupper(algorithms)) %>%
  left_join(drug_targets %>% dplyr::select(drug_ID,source) %>% unique(), by = "drug_ID") %>%
  filter(source==drug_source) %>%
  # keep only the top N predictions that *are available in the right database*
  group_by(sample_ID) %>%
  filter(driver_rank%in%head(sort(unique(driver_rank)),driver_algorithm_top_n)) %>%
  dplyr::select(algorithm,sample_ID, drug_ID)

message(paste0("Reading in rare drug predictions"))
rare_drug_predictions <- fread(paste0("results/benchmark/network_",network_choice,"/rare_drug_predictions.csv")) %>%
  filter(sample_ID %in% samples, toupper(algorithm)==toupper(algorithms)) %>%
  left_join(drug_targets %>% dplyr::select(drug_ID,source) %>% unique(), by = "drug_ID") %>%
  filter(source==drug_source) %>%
  # keep only the top N predictions that *are available in the right database*
  group_by(sample_ID) %>%
  filter(driver_rank%in%head(sort(unique(driver_rank)),driver_algorithm_top_n)) %>%
  dplyr::select(algorithm,sample_ID, drug_ID)


#################################################
# Predict Essentiality Using CGC Tier 1 Drivers #
#################################################


    ##########################################
    # Read In Synthetic Lethality Predictions
    ##########################################
    
    SL <- fread(paste0("benchmark_data/network_",network_choice,"/SL_partners_max_",max_SL,".csv")) %>%
      filter(score>=SL_th) %>%
      dplyr::select(gene_ID=gene1,partner=gene2,SL_rank=rank)
    
    
    #############################
    # Read In LOF/GOF Predictions
    #############################
    
    LOF_GOF <- fread("benchmark_data/LOF_GOF_annotations.csv")
    if("cell_ID"%in%colnames(LOF_GOF)){
      LOF_GOF <- LOF_GOF %>% dplyr::rename(sample_ID=cell_ID)
    }

    ###################################
    # Predict Essential Genes and Drugs
    ###################################
    
    ## CGC
    
    CGC <- read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv") %>%
      filter(Tier==1) %>%
      dplyr::select(gene_ID = `Gene Symbol`) %>%
      filter(gene_ID %in% genes) %>%
      pull(gene_ID) %>%
      unique()
    
    
    ## Altered Genes
    
    mutation <- fread(paste0("benchmark_data/network_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
      filter(gene_ID %in% genes) %>%
      arrange(gene_ID) %>%
      pivot_longer(-gene_ID, names_to = "sample_ID", values_to = "mutated") %>%
      filter(mutated!=0) %>%
      dplyr::select(-mutated)
    
    cnv <- fread(paste0("benchmark_data/network_",network_choice,"/cnv.csv"), select = c("gene_ID",samples)) %>%
      filter(gene_ID %in% genes) %>%
      arrange(gene_ID) %>%
      pivot_longer(-gene_ID, names_to = "sample_ID", values_to = "mutated") %>%
      filter(mutated!=0) %>%
      dplyr::select(-mutated)
    
    # Combine mutation and CNV data to get list of altered genes
    altered_genes <- rbind(mutation,cnv) %>% 
      unique()
    
    CGC_drivers <- altered_genes %>%
      filter(gene_ID %in% CGC) %>%
      mutate(algorithm="CGC_Tier_1") %>%
      dplyr::rename(driver="gene_ID")
    
    CGC_essential_genes <- CGC_drivers %>%
      left_join(LOF_GOF, by = c("driver"="gene_ID","sample_ID"), relationship = "many-to-many") %>%
      # assume that missing annotations are LOF
      mutate(annotation=ifelse(is.na(annotation),"LOF",annotation)) %>%
      # If mutation is GOF, add the gene itself to the target list
      mutate(target = ifelse(annotation %in% c("GOF"), driver, NA)) %>%
      # Add SL partners
      left_join(SL, by = c("driver"="gene_ID"), relationship = "many-to-many") %>%
      # Remove SL partners for GOF mutations
      mutate(partner = ifelse(annotation == "GOF", NA, partner),
             SL_rank = ifelse(annotation == "GOF", NA, SL_rank)) %>%
      unique() %>%
      # Make the SL partner the target for LOF and both
      mutate(target = ifelse(annotation %in% c("LOF", "both"), partner, target))
    
    # Add the original gene_ID as an additional target for genes with both GOF and LOF mutations
    CGC_essential_genes <- CGC_essential_genes %>%
      rbind(CGC_essential_genes %>% 
              filter(annotation == "both") %>% 
              mutate(target = driver, SL_rank = 0, partner = NA) %>%
              unique()) %>%
      group_by(sample_ID, algorithm) %>%
      # For each cell, after arranging by driver rank, then SL_rank, remove duplicated target genes
      filter(!duplicated(target)) %>%
      # Finally remove missing targets (LOF mutations without an SL partner)
      filter(!is.na(target)) %>%
      ungroup() %>%
      arrange(sample_ID) %>%
      dplyr::select(algorithm,sample_ID,gene_ID=target)
    
    CGC_drugs <- CGC_essential_genes %>%
      # Add drugs for targets
      left_join(drug_targets, by = "gene_ID", relationship="many-to-many") %>%
      filter(!is.na(drug_ID)) %>%
      filter(source=="GDSC2")

essentiality <- fread("benchmark_data/gene_effect_z_scores.csv") %>%
  dplyr::rename(sample_ID=cell_ID) %>%
  dplyr::select(sample_ID,gene_ID,gene_effect,uniqueness=weighted_average)

drug_sensitivity <- fread("benchmark_data/drug_sensitivity.csv") %>%
  dplyr::rename(sample_ID=cell_ID) %>%
  filter(sensitivity_source==drug_source) %>%
  dplyr::select(sample_ID,drug_ID,sensitivity_value,uniqueness=weighted_average)

proof_of_concept_plot <- function(type,measure,predictions,ylims,ypos,seed){
  set.seed(seed)
  
  predictions <- predictions %>% mutate(algorithm="PDPA") %>% ungroup()
  
  if(type=="gene"){
    plot_data <- bind_rows(predictions,CGC_essential_genes) %>%
      right_join(essentiality, by = c("sample_ID","gene_ID")) %>%
      left_join(altered_genes %>% mutate(altered=T), by= c("sample_ID", "gene_ID")) %>%
      mutate(altered=ifelse(is.na(altered),F,altered)) %>%
      mutate(algorithm=ifelse(is.na(algorithm)&!altered, "control",algorithm)) %>%
      filter(!is.na(algorithm)) %>%
      dplyr::select(algorithm,sample_ID,prediction=gene_ID,raw=gene_effect,uniqueness)
  }else if(type=="drug"){
    plot_data <- bind_rows(predictions,CGC_drugs) %>%
      right_join(drug_sensitivity, by = c("sample_ID","drug_ID")) %>%
      mutate(algorithm=ifelse(is.na(algorithm), "control",algorithm)) %>%
      dplyr::select(algorithm,sample_ID,prediction=drug_ID,raw=sensitivity_value,uniqueness)
  }
  if(measure=="raw"){
    plot_data <- plot_data %>%
      dplyr::select(algorithm,sample_ID,prediction,measure=raw)
  }else if(measure=="uniqueness"){
    plot_data <- plot_data %>%
      dplyr::select(algorithm,sample_ID,prediction,measure=uniqueness)
  }
  
  plot_data <- plot_data %>%
    mutate(algorithm=factor(algorithm, levels=c("control","CGC_Tier_1","PDPA")))
  
  # get the max sample size of comparison
  n_predictions <- plot_data %>%
    filter(algorithm!="control") %>%
    group_by(algorithm) %>%
    summarise(count=n())
  max_predictions <- max(n_predictions$count)
  
  # make equally sized control comparison
  
  # Keep all rows for predicted essential genes
  keep_rows <- which(plot_data$algorithm!="control")
  # Also keep a random sample of 10000 other rows
  keep_rows <- append(keep_rows, sample(which(!1:nrow(plot_data) %in% keep_rows), max_predictions))
  
  plot_data <- plot_data[keep_rows,]
  
  pairwise_stats <- compare_means(measure~algorithm,data=plot_data,p.adjust.method = "BH",method = "wilcox.test") %>%
    add_significance(p.col = "p.adj", output.col = "p.adj.signif") %>%
    dplyr::rename(group1_tmp=group1, group2_tmp=group2) %>%
    rowwise() %>%
    mutate(group1=min(group1_tmp,group2_tmp), group2=max(group1_tmp,group2_tmp)) %>%
    dplyr::select(-ends_with("tmp"))
  
  effect_size <- rstatix::cohens_d(measure~algorithm,data=plot_data) %>%
    dplyr::rename(group1_tmp=group1, group2_tmp=group2) %>%
    rowwise() %>%
    mutate(group1=min(group1_tmp,group2_tmp), group2=max(group1_tmp,group2_tmp)) %>%
    dplyr::select(-ends_with("tmp"))
  
  stats <- left_join(pairwise_stats, effect_size, by = c(".y.","group1","group2")) %>%
    ungroup() %>%
    dplyr::select(`.y.`, group1, group2, p.adj.signif,effsize) %>%
    mutate(effsize=round(effsize,2))
  
  ggplot(plot_data, aes(x=algorithm, y=measure)) +
    geom_boxplot(aes(fill=algorithm),outlier.shape = NA, notch = T,lwd=0.3) +
    stat_pvalue_manual(stats,label="{p.adj.signif}, D = {effsize}", 
                       y.position = ypos, color = "red", bracket.size = 0.5, tip.length = 0.01, size = 2.5
    ) +
    scale_fill_manual(labels=c("control"="Not Predicted (Background Control)",
                               "CGC_Tier_1"="TARGET-SL + CGC Tier 1",
                               "PDPA"="TARGET-SL + PDPA"),
                      values=c("#A0B1BA","#FF1F5B","#009ADE")) +
    guides(fill=guide_legend("Predictions")) +
    ylab(ifelse(measure=="raw", ifelse(type=="gene", "Raw Gene Effect", "Drug lnIC50"), ifelse(type=="gene", "Gene Uniqueness Index", "Drug Uniqueness Index"))) +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    coord_cartesian(ylim=ylims) +
    theme(text = element_text(size = 10))

}   

p1 <- proof_of_concept_plot(type="gene",measure = "raw", predictions = rare_essential_genes, ylims = c(-3,2), ypos = c(1.5,1.9,1.5), seed=999)
p2 <- proof_of_concept_plot(type="gene",measure = "uniqueness", predictions = rare_essential_genes, ylims = c(-3,3), ypos = c(2.7,3,2.7), seed=999)    
p3 <- proof_of_concept_plot(type="drug",measure = "raw", predictions = rare_drug_predictions, ylims = c(-5,10.7), ypos = c(9.5,10.5,9.5), seed=999)  
p4 <- proof_of_concept_plot(type="drug",measure = "uniqueness", predictions = rare_drug_predictions, ylims = c(-3,3), ypos = c(2.7,3,2.7), seed=999) 

p1 + p2 + p3 + p4 + plot_layout(guides = "collect", nrow = 1) & theme(legend.position = "bottom")

ggsave("plots/benchmark/proof_of_concept.png", width = 20, height = 8)
ggsave("plots/benchmark/proof_of_concept.svg", width = 20, height = 8, device = svglite)
