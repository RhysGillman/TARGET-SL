# format_sysSVM2_results.R

# This code reformats the output from PersonaDrive to be consistent with the other algorithms

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
#suppressPackageStartupMessages (library(biomaRt, quietly = T))

# Handling input arguments
option_list = list(
  make_option(c("-m", "--mode"), type="character", default="benchmark", 
              help="mode (benchmark or predict)", metavar ="Mode"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-c", "--cancertype"), type="character", default="Liver", 
              help="cancer type to analyse", metavar ="Cell Type")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

run_mode <- opt$mode
network_choice <- opt$network
cancer <- opt$cancertype


# Prepare entrez ID mapping file

if(run_mode=="predict"){
  entrez_id_mapping <- fread("data/entrez_ids.csv")
}else if(run_mode=="benchmark"){
  entrez_id_mapping <- read_csv(paste0("benchmark_data/network_", network_choice,"/entrez_ids.csv"))
}


# Get sample list from results files
samples <- list.files(paste0("results/",run_mode,"/network_",network_choice,"/sysSVM2/",cancer,"/")) 
samples <- gsub(".csv","", samples)

suppressWarnings(rm(sysSVM2_results))

for(s in samples){
  sample_drivers <- read_csv(paste0("results/",run_mode,"/network_",network_choice,"/sysSVM2/",cancer,"/",s,".csv"), col_types = cols()) %>%
    dplyr::select(sample_ID=sample,entrez_ID=entrez,rank=sample_rank)
  if(exists("sysSVM2_results")){
    sysSVM2_results <- rbind(sysSVM2_results, sample_drivers)
  }else{
    sysSVM2_results <- sample_drivers
  }
}

sysSVM2_results <- sysSVM2_results %>%
  inner_join(entrez_id_mapping, by = c("entrez_ID"="entrezgene_id"), multiple="all") %>%
  dplyr::select(sample_ID,gene_ID=hgnc_symbol,rank) %>%
  group_by(sample_ID) %>%
  filter(!duplicated(gene_ID)) %>%
  arrange(rank) %>%
  mutate(rank=row_number()) %>%
  ungroup() %>%
  arrange(sample_ID, rank) %>%
  mutate(cancer_type=cancer) %>%
  dplyr::select(cancer_type, sample_ID, driver=gene_ID, rank)

write_csv(sysSVM2_results, paste0("results/",run_mode,"/network_",network_choice,"/sysSVM2/",cancer, ".csv"))
