#predict_drug_sensitivity.r


# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(patchwork, quietly = T))
suppressPackageStartupMessages (library(purrr, quietly = T))
suppressPackageStartupMessages (library(readxl, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(coin, quietly = T))




# Handling input arguments
option_list = list(
  make_option(c("-m", "--mode"), type="character", default="benchmark", 
              help="mode (benchmark or predict)", metavar ="Mode"),
  make_option(c("-s", "--sampleinfo"), type="character", default=NULL, 
              help="path to sample info file", metavar ="Sample Info"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="ALL", 
              help="algorithms to include in comparison separated by spaces, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Threads"),
  make_option(c("-c", "--cancertype"), type="character", default="all", 
              help="Cancer types to include in the analysis separated by semicolons, or 'ALL' (Default)", metavar ="Cancer Type")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

run_mode <- opt$mode
network_choice <- opt$network
algorithms <- opt$algorithms
algorithms <- str_split(algorithms,";") %>% unlist()
threads <- opt$threads
cancer <- opt$cancertype

sample_info_path <- opt$sampleinfo

if(run_mode=="predict"){
  sample_info_path <- opt$sampleinfo
}


if(threads>1){
  cl <- makeCluster(threads, outfile = "log/predict_essential_genes.log")
  registerDoParallel(cl)
}


#############################
# Sample Info
#############################

if(run_mode=="predict"){
  sample_info <- fread(sample_info_path) %>% filter(get_results)
  if(toupper(cancer) != "ALL"){
    sample_info <- sample_info %>% filter(cancer_type == cancer)
  }
  samples <- sample_info$sample_ID %>% unique() %>% sort()
}else if(run_mode=="benchmark"){
  sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv"))
  if(toupper(cancer) != "ALL"){
    sample_info <- sample_info %>% filter(cancer_type == cancer)
  }
  samples <- sample_info$cell_ID %>% unique() %>% sort()
}

###################################
# Read In Predicted Essential Genes
###################################

message(paste0("Reading in essential gene predictions"))
all_essential_genes <- fread(paste0("results/",run_mode,"/network_",network_choice,"/all_essential_genes.csv"))
message(paste0("Reading in rare essential gene predictions"))
rare_essential_genes <- fread(paste0("results/",run_mode,"/network_",network_choice,"/rare_essential_genes.csv"))

###################################
# Read In RandomDrug Predictions
###################################

#if(run_mode=="benchmark"){
#  message(paste0("Reading in randomDrug predictions"))
#  randomDrug <- foreach(result_dir=list.dirs(paste0("results/benchmark/network_",network_choice,"/randomDrug"), recursive = F), .combine = "rbind") %do% {
#      indiv_result <- foreach(result_file=list.files(result_dir),.combine = "rbind") %do% {
#                                fread(paste0(result_dir,"/",result_file))
#    }
#  } %>% dplyr::rename(drug_rank=rank)
#}


###################################
# Drug-Gene Interactions
###################################
if(run_mode=="predict"){
  drug_targets <- fread("data/inhibitory_drug_targets.csv") %>%
    dplyr::select(drug_ID,gene_ID) %>%
    unique() %>%
    # Get the number of gene targets for each drug
    group_by(drug_ID) %>%
    mutate(n_targets=length(gene_ID)) %>%
    ungroup()
}else if(run_mode=="benchmark"){
  drug_targets <- fread(paste0("benchmark_data/network_",network_choice,"/drug_targets.csv")) %>%
    dplyr::select(drug_ID,gene_ID) %>%
    unique() %>%
    # Get the number of gene targets for each drug
    group_by(drug_ID) %>%
    mutate(n_targets=length(gene_ID)) %>%
    ungroup()
}


####################################
# All Predictions
####################################

all_predictions <- all_essential_genes %>%
  dplyr::rename(gene_rank=final_rank) %>%
  # Add drugs for targets
  left_join(drug_targets, by = c("target"="gene_ID"), relationship="many-to-many") %>%
  filter(!is.na(drug_ID)) %>%
  # Get the best gene rank for each drug in case it occurs multiple times
  group_by(algorithm,sample_ID,drug_ID) %>%
  arrange(gene_rank) %>%
  mutate(min_gene_rank=min(gene_rank)) %>%
  ungroup() %>%
  group_by(algorithm,sample_ID) %>%
  # prioritise drugs that target fewer genes
  arrange(min_gene_rank,n_targets) %>%
  #filter(!duplicated(drug_ID)) %>%
  mutate(drug_rank=row_number()) %>%
  ungroup() %>%
  arrange(algorithm,cancer_type,sample_ID,drug_rank)

#if(run_mode=="benchmark"){
#  all_predictions <- all_predictions %>%
#    bind_rows(randomDrug)
#}

write_csv(all_predictions, paste0("results/",run_mode,"/network_",network_choice,"/all_drug_predictions.csv"))

rare_predictions <- rare_essential_genes %>%
  dplyr::rename(gene_rank=final_rank) %>%
  # Add drugs for targets
  left_join(drug_targets, by = c("target"="gene_ID"), relationship="many-to-many") %>%
  filter(!is.na(drug_ID)) %>%
  # Get the best gene rank for each drug in case it occurs multiple times
  group_by(algorithm,sample_ID,drug_ID) %>%
  arrange(gene_rank) %>%
  mutate(min_gene_rank=min(gene_rank)) %>%
  ungroup() %>%
  group_by(algorithm,sample_ID) %>%
  # prioritise drugs that target fewer genes
  arrange(min_gene_rank,n_targets) %>%
  #filter(!duplicated(drug_ID)) %>%
  mutate(drug_rank=row_number()) %>%
  ungroup() %>%
  arrange(algorithm,cancer_type,sample_ID,drug_rank)

#if(run_mode=="benchmark"){
#  all_predictions <- rare_predictions %>%
#    bind_rows(randomDrug)
#}


write_csv(rare_predictions, paste0("results/",run_mode,"/network_",network_choice,"/rare_drug_predictions.csv"))
