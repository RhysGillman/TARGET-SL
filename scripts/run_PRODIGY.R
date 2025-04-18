# run_PRODIGY.R
# This master scripts uses the scripts within scripts/DawnRank to run DawnRank on the validation data

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(MASS, quietly = T))
suppressPackageStartupMessages (library(DESeq2, quietly = T))
suppressPackageStartupMessages (library(igraph, quietly = T))
suppressPackageStartupMessages (library(graphite, quietly = T))
suppressPackageStartupMessages (library(ff, quietly = T))
suppressPackageStartupMessages (library(plyr, quietly = T))
suppressPackageStartupMessages (library(biomaRt, quietly = T))
suppressPackageStartupMessages (library(PCSF, quietly = T))
suppressPackageStartupMessages (library(mixtools, quietly = T))
suppressPackageStartupMessages (library(ggplot2, quietly = T))
suppressPackageStartupMessages (library(cowplot, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-m", "--mode"), type="character", default="benchmark", 
              help="mode (benchmark or predict)", metavar ="Mode"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-c", "--cancertype"), type="character", default="Liver", 
              help="cancer type to analyse", metavar ="Cell Type"),
  make_option(c("-t", "--threads"), type="character", default=4, 
              help="Number of threads to use", metavar ="Threads")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

run_mode <- opt$mode
network_choice <- opt$network
cancer_type <- opt$cancertype
threads <- opt$threads

if(run_mode=="benchmark"){

#############################
# Functions
#############################

# Scripts come from https://github.com/Shamir-Lab/PRODIGY

for(s in list.files("scripts/PRODIGY/")){
  source(paste0("scripts/PRODIGY/",s))
}



#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv")) %>% filter(lineage==cancer_type)
samples <- sample_info$cell_ID %>% sort()

#############################
# Prepare Input Data
#############################

rna <- fread(paste0("benchmark_data/network_",network_choice,"/counts.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID") %>% 
  as.matrix()

mutation <-  fread(paste0("benchmark_data/network_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID") %>% 
  as.matrix()

if(network_choice=="own"){
  
  network <- fread("data/own_networks/PRODIGY_personadrive.csv") %>%
    as.matrix()
  
}else{
  network <- fread(paste0("benchmark_data/network_",network_choice,"/network_undirected.csv"))
  colnames(network) <- c("source", "destination", "score")
  network <- network %>%
    mutate(score = score/1000) %>%
    as.matrix()
}
genes <- Reduce(intersect, list(
  rownames(rna),
  rownames(mutation),
  append(network[,"source"], network[,"destination"])
))

rna <- rna[genes, samples]
mutation <- mutation[genes, samples]


# Pathways

if(!file.exists("data/PRODIGY_pathways.rds")){
  pathways <- get_pathway_list_from_graphite(source = "reactome",minimal_number_of_nodes = 10,num_of_cores = threads)
  write_rds(pathways, "data/PRODIGY_pathways.rds")
}else{
  pathways <- read_rds("data/PRODIGY_pathways.rds")
}

#############################
# Running PRODIGY
#############################

res_df <- data.frame(cancer_type = character(), sample = character(), driver = character(), rank = factor())

for(sample in samples){
  
  message(paste0("Beginning PRODIGY analysis for sample: ", sample, " (", which(samples==sample), " of ", length(samples), ")"))
  #Set all other samples as "normal" samples
  sample_origins <- ifelse(samples == sample, "tumor", "normal")
  #Getting DEGs in sample
  message("Calulcating DEGs...")
  
  DEGs <- get_DEGs(expression_matrix = rna,
                   samples = sample,
                   sample_origins = sample_origins,
                   beta = 2,
                   gamma = 0.05)
  diff_genes <- DEGs[[sample]]
  
  #Getting mutated genes in samples
  message("Retrieving Mutated Genes...")
  sample_mutations <- names(mutation[mutation[,sample] == 1,sample])
  #Running PRODIGY
  
  message("Running PRODIGY")
  
  sink("log/PRODIGY_log.txt")
  res <- PRODIGY(
    mutated_genes = sample_mutations,
    expression_matrix = rna,
    network = network,
    sample = sample,
    diff_genes = diff_genes,
    alpha = 0.05,
    pathway_list = pathways,
    num_of_cores = threads,
    sample_origins = sample_origins,
    write_results = F,
    results_folder = "./",
    beta = 2,
    gamma = 0.05,
    delta = 0.05
  )
  
  
  res <- analyze_PRODIGY_results(res)
  
  sink()
  
  if(is_empty(res)){
    res <- c(NA)
  }else{
    res <- res[[1]]
  }
  
  if(is_empty(res)){
    res <- c(NA)
  }
  
  
  temp_res_df <- data.frame(cancer_type, sample_ID = sample, driver = res, rank = 1:length(res))
  
  res_df <- rbind(res_df, temp_res_df)
}

#############################
# Save Result
#############################
if(!dir.exists(paste0("results/",run_mode,"/network_",network_choice,"/PRODIGY"))){
  dir.create(paste0("results/",run_mode,"/network_",network_choice,"/PRODIGY"))
}

write_csv(res_df, paste0("results/",run_mode,"/network_",network_choice,"/PRODIGY/",cancer_type,".csv"))

}