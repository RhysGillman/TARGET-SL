#get_consensus_drivers.r

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-m", "--mode"), type="character", default="benchmark", 
              help="mode (benchmark or predict)", metavar ="Mode"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="OncoImpact;PersonaDrive;sysSVM2;DawnRank", 
              help="algorithms to include in consensus separated by semicolons, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-M", "--method"), type="character", default="custom", 
              help="Rank aggregation method to be used", metavar ="Rank Aggregation Method"),
  make_option(c("-c", "--cancertype"), type="character", default="all", 
              help="Cancer types to analyse separated by semicolons, or 'ALL' (Default)", metavar ="Cancer Type"),
  make_option(c("-s", "--sampleinfo"), type="character", default=NULL, 
              help="path to sample info file", metavar ="Sample Info")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

run_mode <- opt$mode
network_choice <- opt$network
algorithms <- opt$algorithms
RAmethod <- opt$method
cancer <- opt$cancertype
if(run_mode=="predict"){
  sample_info_path <- opt$sampleinfo
}

algorithms <- str_split(algorithms,";") %>% unlist()
cancer <- str_split(cancer,";") %>% unlist()


#############################
# Read In Results
#############################

source("scripts/read_driver_results.R")

aggregated_results <- read_driver_results(run_mode,algorithms,cancer)

#############################
# Rank Aggregation
#############################

if(RAmethod=="custom"){
  exp <- 2
  topn <- 20
  
  custom_transformation <- function(x){
    result <- (topn + 1 - x)^exp
    return(result)
  }
  
  
  # Take the top predictions only
  filtered_results <- aggregated_results %>% filter(rank <= topn)
  
  
  custom_scores <- filtered_results %>%
    #dplyr::select(-lineage) %>%
    pivot_wider(names_from = "algorithm", values_from = "rank")
  custom_scores <- custom_scores %>%
    cbind(
      
      custom_scores %>% 
        dplyr::select(-c(cancer_type,sample_ID,driver)) %>%
        mutate(across(where(is.numeric), custom_transformation)) %>%
        setNames(paste0("score_",names(.)))
      
    ) %>%
    
    mutate(final_score = rowSums(across(starts_with("score_")), na.rm=T))
  
  
  
  
  consensus_drivers <- custom_scores %>%
    dplyr::select(cancer_type,sample_ID, driver, final_score) %>%
    group_by(sample_ID) %>%
    arrange(desc(final_score)) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    dplyr::select(-final_score) %>%
    mutate(algorithm = "consensus")
  
  suppressWarnings(dir.create(paste0("results/",run_mode,"/network_",network_choice,"/consensus")))
  
  write_csv(consensus_drivers,paste0("results/",run_mode,"/network_",network_choice,"/consensus/consensus_drivers.csv"))
  
  
  
}


