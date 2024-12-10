#format_combined_de_novo_methods_results.R
# This code reformats the output from the combined de novo methods to be consistent with the other algorithms

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-m", "--mode"), type="character", default="benchmark", 
              help="mode (benchmark or predict)", metavar ="Mode"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-c", "--cancertype"), type="character", default="Liver", 
              help="Cancer type to analyse", metavar ="Cancer Type")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

run_mode <- opt$mode
network_choice <- opt$network
cancer_type <- opt$cancertype

if(run_mode=="benchmark"){

#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv")) %>% filter(lineage==cancer_type)
samples <- sample_info$cell_ID %>% sort()


#############################
# Genetically Altered Genes
#############################
# Mutation

mutation <- fread(paste0("benchmark_data/network_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID")

# CNV

cnv <- fread(paste0("benchmark_data/network_",network_choice,"/cnv.csv"), select = c("gene_ID",samples)) %>%
  column_to_rownames("gene_ID")

# Get a list of altered genes for each sample

mutated <- which(mutation==1, arr.ind = TRUE) %>% 
  as.data.frame() %>%
  mutate(sample = colnames(mutation)[col]) %>%
  mutate(gene_ID = rownames(mutation)[row]) %>%
  dplyr::select(sample, gene_ID)

copy_altered <- which(cnv!=0, arr.ind = TRUE) %>% 
  as.data.frame() %>%
  mutate(sample = colnames(cnv)[col]) %>%
  mutate(gene_ID = rownames(cnv)[row]) %>%
  dplyr::select(sample, gene_ID)

altered_genes <- rbind(mutated,copy_altered)
rownames(altered_genes) <- NULL
altered_genes <- unique(altered_genes) %>%
  mutate(is_altered = TRUE)

#############################
# Results
#############################

# Combine results into a single dataframe, add whether the driver nodes are mutated or not

files <- list.files(paste0("results/",run_mode,"/network_",network_choice,"/combined_de_novo_methods/",cancer_type)) %>% stringr::str_subset("result.csv$")

suppressWarnings(rm(de_novo_networks_results))
for(file in files){
  method <- strsplit(file, "_") %>%
    unlist() %>%
    head(2)
  
  indiv_result <- read_csv(paste0("results/",run_mode,"/network_",network_choice,"/combined_de_novo_methods/",cancer_type,"/",file), show_col_types = FALSE) %>% 
    column_to_rownames("Row")
  
  driver_nodes <- which(indiv_result==1, arr.ind = TRUE) %>% 
    as.data.frame() %>%
    mutate(sample = colnames(indiv_result)[col]) %>%
    mutate(gene_ID = rownames(indiv_result)[row]) %>%
    dplyr::select(sample, gene_ID)
  
  rownames(driver_nodes) <- NULL
  
  driver_nodes <- driver_nodes %>%
    mutate(network_method = method[1], control_method = method[2]) %>%
    unique() %>%
    left_join(altered_genes, by = c("sample", "gene_ID")) %>%
    mutate(is_altered = ifelse(is.na(is_altered), FALSE, is_altered))
  
  if(!exists("de_novo_networks_results")){
    de_novo_networks_results <- driver_nodes
  }else{
    de_novo_networks_results <- rbind(de_novo_networks_results, driver_nodes)
  }
  
}

# Collect the in and out degrees for ranking

methods <- list.files(paste0("results/",run_mode,"/network_",network_choice,"/combined_de_novo_methods/",cancer_type)) %>% 
  stringr::str_subset("deg.csv$") %>% gsub(pattern = "_.+$",replacement = "") %>% unique()


suppressWarnings(rm(de_novo_networks_degrees))
for(method in methods){
  
  indiv_out <- read_csv(paste0("results/",run_mode,"/network_",network_choice,"/combined_de_novo_methods/",cancer_type,"/",method, "_out_deg.csv"),show_col_types = FALSE)
  indiv_out[is.na(indiv_out)] <- 0
  indiv_out <- indiv_out %>% pivot_longer(!gene_ID, names_to = "sample", values_to = "out_degree")
  
  indiv_in <- read_csv(paste0("results/",run_mode,"/network_",network_choice,"/combined_de_novo_methods/",cancer_type,"/",method, "_in_deg.csv"),show_col_types = FALSE)
  indiv_in[is.na(indiv_in)] <- 0
  indiv_in <- indiv_in %>% pivot_longer(!gene_ID, names_to = "sample", values_to = "in_degree")
  indiv <- full_join(indiv_out,indiv_in, by = c("gene_ID", "sample")) %>%
    mutate(network_method = method)
  
  if(!exists("de_novo_networks_degrees")){
    de_novo_networks_degrees <- indiv
  }else{
    de_novo_networks_degrees <- rbind(de_novo_networks_degrees, indiv)
  }
  
}

de_novo_networks_results <- de_novo_networks_results %>% 
  left_join(de_novo_networks_degrees, by = c("gene_ID", "sample", "network_method")) %>%
  filter(is_altered) %>%
  arrange(desc(out_degree)) %>%
  group_by(sample, network_method, control_method) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  mutate(algorithm = paste0(network_method, "_", control_method)) %>%
  mutate(cancer_type) %>%
  dplyr::select(cancer_type, sample_ID = sample, driver = gene_ID, rank, algorithm) %>%
  arrange(algorithm, sample_ID, rank)


write_csv(de_novo_networks_results, paste0("results/",run_mode,"/network_",network_choice,"/combined_de_novo_methods/",cancer_type,".csv"))
}