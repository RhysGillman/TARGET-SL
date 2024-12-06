# badDriver.R

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-c", "--cancertype"), type="character", default="all", 
              help="Cancer types to include in the analysis separated by semicolons, or 'ALL' (Default)", metavar ="Cancer Type"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Rank Aggregation Method"),
  make_option(c("-s", "--simnumber"), type="integer", default=10, 
              help="the number of simulations to run", metavar ="Simulation Number"),
  make_option(c("-S", "--seed"), type="integer", default=9999, 
              help="Set the seed for reproducible randomised predictions", metavar ="Seed")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
cancer <- opt$cancertype
cancer <- str_split(cancer,";") %>% unlist()
threads <- opt$threads
n_sims <- opt$simnumber
custom_seed <- opt$seed


if(threads>1){
  #registerDoParallel(cores=threads)
  cl <- makeCluster(threads, outfile = "log/random_predictions.log")
  registerDoParallel(cl)
}

#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv"))
if(toupper(cancer)!="ALL"){
  sample_info <- sample_info %>% filter(lineage %in% cancer)
}
cancer_types <- sample_info$lineage %>% unique()


for(c in cancer_types){
  
  samples <- sample_info %>% filter(lineage==c) %>% pull(cell_ID)
  
  #############################
  # Create Directories
  #############################
  
  if(!dir.exists(paste0("results/benchmark/network_",network_choice,"/randomDriver/",c))){
    dir.create(paste0("results/benchmark/network_",network_choice,"/randomDriver/",c), recursive = T)
  }
  
  if(!dir.exists(paste0("results/benchmark/network_",network_choice,"/randomDrug/",c))){
    dir.create(paste0("results/benchmark/network_",network_choice,"/randomDrug/",c), recursive = T)
  }
    
  #############################
  # Mutation Data
  #############################
  
  # Mutations
  
  mutation <- fread(paste0("benchmark_data/network_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
    arrange(gene_ID) %>%
    column_to_rownames("gene_ID")
  
  cnv <- fread(paste0("benchmark_data/network_",network_choice,"/cnv.csv"), select = c("gene_ID",samples)) %>%
    column_to_rownames("gene_ID")
  
  mutation <- mutation[sort(rownames(mutation)),sort(colnames(mutation))]
  cnv <- cnv[sort(rownames(cnv)),sort(colnames(cnv))]
  
  
  if(!all(colnames(mutation)==colnames(cnv)) & all(rownames(mutation) == rownames(cnv))){
    stop("Error: colnames/rownames of cnv and mutation files do not match")
  }
  
  
  altered_genes <- cnv != 0 | mutation != 0
  altered_genes[] <- as.integer(altered_genes)
  
  
  #############################
  # randomDriver prediction
  #############################
  
  
  for(i in 1:n_sims){
    
    random_drivers <- foreach(sample=samples, .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {

      sample_altered_genes <- rownames(altered_genes)[which(altered_genes[,sample]==1)] %>% 
        # randomise order
        sample(replace = F)
    
      
      data.frame(cancer_type = c,sample_ID = sample, driver = sample_altered_genes) %>% mutate(rank = row_number(), algorithm=paste0("randomDriver_",i))
      
      
    }
    write_csv(random_drivers, paste0("results/benchmark/network_",network_choice,"/randomDriver/",c,"/randomDriver_",i,".csv"))
  }
  
  #############################
  ## Get Available Drugs
  #############################
  
  
  all_drugs <- read_csv("benchmark_data/drug_sensitivity.csv") %>%
    dplyr::select(cell_ID,drug_ID) %>%
    unique() %>%
    filter(cell_ID %in% samples) %>%
    group_by(cell_ID) %>%
    summarise(drugs = list(drug_ID)) %>%
    deframe()
  
  #############################
  ## BadDrug Prediction
  #############################
  
  
  for(i in 1:n_sims){
    
    bad_drugs <- foreach(sample=samples[samples %in% names(all_drugs)], .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
      
      sample_drugs <- all_drugs[[sample]] %>%
        sample(replace=F)
      
      data.frame(cancer_type = c,sample_ID = sample, drug_ID = sample_drugs) %>% mutate(rank = row_number(), algorithm=paste0("randomDrug_",i))
      
      
    }
    
    
    write_csv(bad_drugs, paste0("results/benchmark/network_",network_choice,"/randomDrug/",c,"/randomDrug_",i,".csv"))
    
  }
  
  
}
