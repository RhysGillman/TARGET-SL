#predict_essential_genes.r


# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
#suppressPackageStartupMessages (library(ggpubr, quietly = T))
suppressPackageStartupMessages (library(coin, quietly = T))
#suppressPackageStartupMessages (library(ggpointdensity, quietly = T))
#suppressPackageStartupMessages (library(svglite, quietly = T))
suppressPackageStartupMessages (library(ggh4x, quietly = T))

# Handling input arguments
option_list = list(
  make_option(c("-m", "--mode"), type="character", default="benchmark", 
              help="mode (benchmark or predict)", metavar ="Mode"),
  make_option(c("-s", "--sampleinfo"), type="character", default=NULL, 
              help="path to sample info file", metavar ="Sample Info"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="ALL", 
              help="algorithms to include in comparison separated by semicolons, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
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
cancer <- str_split(cancer,";") %>% unlist()


if(run_mode=="predict"){
  sample_info_path <- opt$sampleinfo
}
# SL-partner score filtering threshold (default 0.5)
SL_th <- 0.5
# pvalue filtering threshold for mann-whitney test for rare essential genes
pval <- 0.05
# max number of SL partners to consider for each gene
max_SL <- 5


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

#############################
# Read In Results
#############################

source("scripts/read_driver_results.R")

aggregated_results <- read_driver_results(run_mode,algorithms,cancer)

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

####################################
# All Predicted Essential Genes    #
####################################


all_essential <- aggregated_results %>%
  dplyr::rename(driver_rank=rank) %>%
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
all_essential <- all_essential %>%
  rbind(all_essential %>% 
          filter(annotation == "both") %>% 
          mutate(target = driver, SL_rank = 0, partner = NA) %>%
          unique()) %>%
  arrange(driver_rank,SL_rank) %>%
  #filter(!is.na(target)) %>%
  group_by(sample_ID, algorithm) %>%
  # For each cell, after arranging by driver rank, then SL_rank, remove duplicated target genes
  filter(!duplicated(target)) %>%
  # Finally remove missing targets (LOF mutations without an SL partner)
  filter(!is.na(target)) %>%
  mutate(final_rank = row_number()) %>%
  ungroup() %>%
  arrange(cancer_type,algorithm, sample_ID, final_rank)


write_csv(all_essential, paste0("results/",run_mode,"/network_",network_choice,"/all_essential_genes.csv"))


#####################################
# Rare Predicted Essential Genes    #
#####################################

  
all_algs <- unique(all_essential$algorithm)

message("Running Mann-Whitney U tests to identify rare essential genes, check status in 'log/predict_essential_genes.log'")

# Loop through each algorithm
rare_test_essential <- foreach(alg=all_algs, .combine = "rbind", .packages = c("tidyverse","coin", "foreach")) %dopar% {
  
  # Get all predicted essential genes for the algorithm
  tmp1 <- all_essential %>% filter(algorithm==alg)
  
  # Loop through each sample
  sample_result <- foreach(sample=unique(tmp1$sample_ID), .combine = "rbind") %do% {
    
    print(paste0("Mann-Whitney U Test to find rare predicted essential genes. Algorithm: ", 
                 "(",which(all_algs == alg),"/",length(all_algs),"). Sample: "
                 , "(",which(unique(tmp1$sample_ID) == sample),"/",length(unique(tmp1$sample_ID)),")"))
    
    # Get the predicted essential genes for an individual sample
    sample_predictions <- tmp1 %>% filter(sample_ID==sample) 
    
    # Get the lineage of the sample
    lin <- sample_predictions %>% pull(cancer_type) %>% head(1)
    
    sample_predictions <- sample_predictions %>% pull(target)
    
    # Get the maximum number of predicted essential genes by that algorithm in that lineage
    max_pred <- tmp1 %>% filter(cancer_type==lin) %>% group_by(sample_ID) %>% summarise(count=n()) %>% pull(count) %>% max()
    
    # Get ranks of each predicted essential gene from all samples in the lineage
    lineage_ranks <- tmp1 %>% 
      filter(cancer_type==lin) %>%
      dplyr::select(sample_ID,target,final_rank) %>%
      pivot_wider(names_from = "target", values_from = "final_rank", values_fill = max_pred+1)
    
    
    gene_result <- foreach(gene=sample_predictions, .combine = "rbind") %do% {
      
      # Get the rank of a predicted essential gene in the sample compared with others in the lineage
      sample_rank <- lineage_ranks[c("sample_ID",gene)] %>%
        mutate(sample_ID=ifelse(sample_ID!=sample, "other", "sample")) %>%
        mutate(sample_ID=factor(sample_ID, levels = c("sample","other"))) %>%
        dplyr::rename(rank=2) %>%
        arrange(sample_ID)
      
      
      stat <- wilcox_test(rank~sample_ID, data=sample_rank, conf.level=0.95,ties.method="mid-ranks", alternative = "less", paired=FALSE)
      
      data.frame(
        algorithm=alg,
        cancer_type=lin,
        sample_ID=sample,
        target=gene,
        final_rank=sample_rank[1,"rank"],
        pvalue=pvalue(stat),
        effect_size=as.numeric(stat@statistic@linearstatistic),
        expectation=as.numeric(expectation(stat))
      )
      
      
    }
    
  }
  
  
  sample_result %>% mutate(padj=p.adjust(sample_result$pvalue, method = "BH"))
  
  
}


rare_essential <- all_essential %>%
  left_join(rare_test_essential, by = c("algorithm","cancer_type", "sample_ID","target","final_rank"="rank")) %>%
  filter(padj < pval) %>%
  group_by(algorithm,cancer_type,sample_ID) %>%
  arrange(final_rank) %>%
  mutate(final_rank=row_number()) %>%
  ungroup() %>%
  arrange(algorithm,cancer_type,sample_ID,final_rank)

write_csv(rare_essential, paste0("results/",run_mode,"/network_",network_choice,"/rare_essential_genes.csv"))
