# run_DawnRank.R
# This master scripts uses the scripts within scripts/DawnRank to run DawnRank on the validation data

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(snowfall, quietly = T))
suppressPackageStartupMessages (library(maxstat, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-m", "--mode"), type="character", default="benchmark", 
              help="mode (benchmark or predict)", metavar ="Mode"),
  make_option(c("-s", "--sampleinfo"), type="character", default=NULL, 
              help="path to sample info file", metavar ="Sample Info"),
  make_option(c("-r", "--rna"), type="character", default=NULL, 
              help="path to rna file", metavar ="RNA Data"),
  make_option(c("-d", "--dna"), type="character", default=NULL, 
              help="path to mutations file", metavar ="Mutation Data"),
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

if(run_mode=="predict"){
  sample_info_path <- opt$sampleinfo
  rna_path <- opt$rna
  dna_path <- opt$dna
}



#############################
# Functions
#############################

# Scripts come from https://github.com/MartinFXP/DawnRank

for(s in list.files("scripts/DawnRank/")){
  source(paste0("scripts/DawnRank/",s))
}

#############################
# Network
#############################

# read in chosen GIN

if(network_choice == "own"){
  
  load("data/own_networks/dawnrank.rda")
  network <- originalUnweighted
  rm(originalUnweighted)
  
  
}else{
  
  # Need to convert 2-column dataframe of interactions into a complete matrix
  network <-  fread(paste0("benchmark_data/network_",network_choice,"/network_directed.csv")) %>%
    mutate(confidence = ifelse(confidence>0, 1, 0)) %>%
    pivot_wider(names_from = 2, values_from = 3) %>%
    complete(protein_1 = names(.)[-1]) %>%
    data.table::transpose(make.names = "protein_1", keep.names = "protein_2") %>%
    complete(protein_2 = names(.)[-1]) %>%
    data.table::transpose(make.names = "protein_2", keep.names = "protein_1") %>%
    column_to_rownames("protein_1") %>%
    mutate_all(~replace_na(., 0)) %>%
    as.matrix()
}


if(run_mode=="predict"){

    #############################
    # Sample Info
    #############################

    sample_info <- fread(sample_info_path) %>% filter(cancer_type==cancer, get_results | !is.na(normal_sample_ID))
    samples <- sample_info %>% pull(sample_ID) %>% unique() %>% sort()
    results_samples <- sample_info %>% filter(get_results) %>% pull(sample_ID) %>% unique() %>% sort()

    normal_samples <- sample_info %>% pull(normal_sample_ID) %>% na.omit()
    tumour_samples <- sample_info %>% filter(get_results) %>% pull(tumour_sample_ID) %>% na.omit()

    #############################
    # Prepare Input Data
    #############################

    tumour_rna <- fread(rna_path, select = c("gene_ID",tumour_samples)) %>% 
      arrange(gene_ID) %>%
      column_to_rownames("gene_ID") %>% 
      as.matrix()

    tumour_cols <- data.frame(tumour_sample_ID= colnames(tumour_rna)) %>% 
      left_join(sample_info %>% dplyr::select(sample_ID, tumour_sample_ID), by="tumour_sample_ID") %>%
      unique() %>%
      pull(sample_ID)
    colnames(tumour_rna) <- tumour_cols

    normal_rna <- fread(rna_path, select = c("gene_ID",normal_samples)) %>% 
      arrange(gene_ID) %>%
      column_to_rownames("gene_ID") %>% 
      as.matrix()

    normal_cols <- data.frame(normal_sample_ID= colnames(normal_rna)) %>% 
      left_join(sample_info %>% dplyr::select(sample_ID, normal_sample_ID), by="normal_sample_ID") %>%
      unique() %>%
      pull(sample_ID)
    colnames(normal_rna) <- normal_cols

    mutation <-  fread(dna_path, select = c("gene_ID",samples)) %>% 
      arrange(gene_ID) %>%
      column_to_rownames("gene_ID") %>% 
      as.matrix()


    genes <- Reduce(intersect, list(
      rownames(normal_rna),
      rownames(tumour_rna),
      rownames(mutation),
      rownames(network),
      colnames(network)
    ))

    # Make sure that all columns and rows display data in the same order
    normal_rna <- normal_rna[genes,]
    tumour_rna <- tumour_rna[genes,]
    mutation <- mutation[genes, results_samples]
    network <- network[genes,genes]

    #############################
    # Running DawnRank
    #############################

    damping <- dawnDamping(network, 3)
    dawnMatrix <- DawnMatrix(network) 

    res_df <- data.frame(Gene = character(), Patient = character(), Rank = numeric(), PercentRank = numeric())

    for(sample in results_samples){
      
      message(paste0("Beginning DawnRank analysis for sample: ", sample, " (", which(results_samples==sample), " of ", length(results_samples), ")"))
      
      normalizedDawn <- DawnNormalize(tumorMat = tumour_rna, normalMat = normal_rna)
      
      dawn <- Dawn(dawnMatrix, normalizedDawn[,sample], mutation[,sample], damping, maxit = 100, patientTag = sample, epsilon = 1e-04)
      
      #mutated_genes <- dawn$summaryOutput
      
      indiv_res <- dawn$mutatedRanks %>% arrange(desc(PercentRank))
      
      res_df <- rbind(res_df, indiv_res)
    }

    res_df <- res_df %>%
      dplyr::select(sample_ID = Patient, gene_ID = Gene, PercentRank) %>%
      # Converting PercentRank into ordinal gene ranks
      arrange(desc(PercentRank)) %>%
      group_by(sample_ID) %>%
      mutate(rank = row_number(), cancer_type = cancer) %>%
      dplyr::select(-PercentRank) %>%
      arrange(sample_ID) %>%
      dplyr::select(cancer_type, sample_ID, driver = gene_ID, rank)


}else if(run_mode == "benchmark"){

    #############################
    # Sample Info
    #############################

    sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv")) %>% filter(lineage==cancer) %>%
      dplyr::rename(cancer_type=lineage)
    samples <- sample_info$cell_ID %>% sort()
    results_samples <- samples

    #############################
    # Prepare Input Data
    #############################

    rna <- fread(paste0("benchmark_data/network_",network_choice,"/tpm.csv"), select = c("gene_ID",samples)) %>% 
      arrange(gene_ID) %>%
      column_to_rownames("gene_ID") %>% 
      as.matrix()

    mutation <-  fread(paste0("benchmark_data/network_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
      arrange(gene_ID) %>%
      column_to_rownames("gene_ID") %>% 
      as.matrix()

    
    genes <- Reduce(intersect, list(
      rownames(rna),
      rownames(mutation),
      rownames(network),
      colnames(network)
    ))

    # Make sure that all columns and rows display data in the same order
    rna <- rna[genes, samples]
    mutation <- mutation[genes, samples]
    network <- network[genes,genes]

    #############################
    # Running DawnRank
    #############################

    damping <- dawnDamping(network, 3)
    dawnMatrix <- DawnMatrix(network) 

    res_df <- data.frame(Gene = character(), Patient = character(), Rank = numeric(), PercentRank = numeric())

    for(sample in samples){
      
      message(paste0("Beginning DawnRank analysis for sample: ", sample, " (", which(samples==sample), " of ", length(samples), ")"))
      
      tumour_rna <- rna[,sample] %>% as.matrix()
      colnames(tumour_rna) <- sample
      normal_rna <- rna[,samples[samples != sample]]
      
      normalizedDawn <- DawnNormalize(tumorMat = tumour_rna, normalMat = normal_rna)
      colnames(normalizedDawn) <- sample
      
      dawn <- Dawn(dawnMatrix, normalizedDawn[,sample], mutation[,sample], damping, maxit = 100, patientTag = sample, epsilon = 1e-04)
      
      
      indiv_res <- dawn$mutatedRanks %>% arrange(desc(PercentRank))
      
      res_df <- rbind(res_df, indiv_res)
    }

    res_df <- res_df %>%
      dplyr::select(sample_ID = Patient, gene_ID = Gene, PercentRank) %>%
      # Converting PercentRank into ordinal gene ranks
      arrange(desc(PercentRank)) %>%
      group_by(sample_ID) %>%
      mutate(rank = row_number(), cancer_type = cancer) %>%
      dplyr::select(-PercentRank) %>%
      arrange(sample_ID) %>%
      dplyr::select(cancer_type, sample_ID, driver = gene_ID, rank)

}


#############################
# Save Result
#############################
if(!dir.exists(paste0("results/",run_mode,"/network_",network_choice,"/DawnRank"))){
  dir.create(paste0("results/",run_mode,"/network_",network_choice,"/DawnRank"))
}

write_csv(res_df, paste0("results/",run_mode,"/network_",network_choice,"/DawnRank/",cancer,".csv"))






