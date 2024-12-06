#prepare_OncoImpact_data.R
# This code prepares the input data to run OncoImpact and stores it in tmp/

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-m", "--mode"), type="character", default="benchmark", 
              help="mode (benchmark or predict)", metavar ="Mode"),
  make_option(c("-s", "--sampleinfo"), type="character", default=NULL, 
              help="path to sample info file", metavar ="Sample Info"),
  make_option(c("-r", "--rna"), type="character", default=NULL, 
              help="path to rna file", metavar ="Sample Info"),
  make_option(c("-d", "--dna"), type="character", default=NULL, 
              help="path to mutations file", metavar ="Network"),
  make_option(c("-v", "--cnv"), type="character", default=NULL, 
              help="path to mutations file", metavar ="Network"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-c", "--cancertype"), type="character", default="Liver", 
              help="cancer type to analyse", metavar ="Cell Type"),
  make_option(c("-t", "--threads"), type="character", default=4, 
              help="Number of threads to use", metavar ="Threads"),
  make_option(c("-w", "--workdir"), type="character", default=NULL, 
              help="working directory", metavar ="Working Directory"),
  make_option(c("-T", "--tmp"), type="character", default=NULL, 
              help="temporary directory with symbolic link permission", metavar ="Tmp")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

run_mode <- opt$mode
network_choice <- opt$network
cancer <- opt$cancertype
threads <- opt$threads
WORK_DIR <- opt$workdir
TMP_DIR <- opt$tmp

if(run_mode=="predict"){
  sample_info_path <- opt$sampleinfo
  rna_path <- opt$rna
  dna_path <- opt$dna
  cnv_path <- opt$cnv
}

#############################
# Network
#############################

if(network_choice=="own"){
  network <- fread("data/own_networks/OncoImpact.txt", select = c("Gene1","Gene2"))
}else{
  network <- fread(paste0("benchmark_data/network_",network_choice,"/network_undirected.csv")) %>%
    dplyr::select(1,2)
}


if(run_mode=="predict"){

  #############################
  # Sample Info
  #############################
  
  # check which samples have all required info available
  all_data_avail <- Reduce(intersect, list(
    gsub("_[nt]$","",fread(rna_path, header = F, nrows = 1)[1,]) %>% unlist(),
    fread(dna_path, header = F, nrows = 1)[1,] %>% unlist(),
    fread(cnv_path, header = F, nrows = 1)[1,] %>% unlist()
  )) %>% unlist()
  
  
  sample_info <- fread(sample_info_path) %>% filter(cancer_type==cancer, sample_ID %in% all_data_avail)
  
  rna_samples <- sample_info %>% filter(get_results | !is.na(normal_sample_ID)) %>% pull(sample_ID) %>% unique() %>% sort()
  
  results_samples <- sample_info %>% filter(get_results) %>% pull(sample_ID) %>% unique() %>% sort()
  
  if(length(results_samples) < 10){
    set.seed(999)
    extra_samples <- sample(sample_info %>% filter(!get_results) %>%pull(sample_ID) %>% unique(), 10-length(results_samples))
    results_samples <- append(results_samples, extra_samples) %>% sort()
  }
  
  normal_samples <- sample_info %>% pull(normal_sample_ID) %>% na.omit()
  tumour_samples <- sample_info %>% filter(sample_ID %in% results_samples) %>% pull(tumour_sample_ID) %>% na.omit() %>% sort()
  
  #############################
  # Prepare Input Data
  #############################
  
  # RNA
  
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
  
  #RNA data needs to be converted to log2 fold change data
  #To do this, compare each sample vs the mean of all others
  for(sample in results_samples){
    sample_rna <- tumour_rna[,sample]
    mean_rna <- apply(normal_rna,1,mean)
    diff_rna <- sample_rna - mean_rna %>% as.data.frame()
    colnames(diff_rna) <- sample
    if(sample == results_samples[1]){
      l2fc_rna <- diff_rna 
    }else{
      l2fc_rna <- cbind(l2fc_rna, diff_rna)
    }
  }
  rm(tumour_rna, normal_rna, diff_rna, mean_rna)
  
  # Mutations
  
  mutation <-  fread(dna_path, select = c("gene_ID",results_samples)) %>% 
    arrange(gene_ID) %>%
    column_to_rownames("gene_ID")
  
  # CNV
  
  cnv <- fread(cnv_path, select = c("gene_ID",results_samples)) %>%
    column_to_rownames("gene_ID")
  
}else if(run_mode == "benchmark"){
  
  #############################
  # Sample Info
  #############################
  
  sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv")) %>% filter(lineage==cancer)
  samples <- sample_info$cell_ID %>% sort()
  results_samples <- samples
  
  #############################
  # Prepare Input Data
  #############################
  
  # RNA
  
  rna <- fread(paste0("benchmark_data/network_",network_choice,"/tpm.csv"), select = c("gene_ID",samples)) %>% 
    arrange(gene_ID) %>%
    column_to_rownames("gene_ID") 
  
  #RNA data needs to be converted to log2 fold change data
  #To do this, compare each sample vs the mean of all others
  for(sample in samples){
    tumour_rna <- rna[,sample] %>% as.matrix()
    normal_rna <- rna[,samples[samples != sample]]
    mean_rna <- apply(normal_rna,1,mean)
    diff_rna <- tumour_rna - mean_rna %>% as.data.frame()
    colnames(diff_rna) <- sample
    if(sample == samples[1]){
      l2fc_rna <- diff_rna 
    }else{
      l2fc_rna <- cbind(l2fc_rna, diff_rna)
    }
  }
  rm(rna, tumour_rna, normal_rna, diff_rna, mean_rna)
  
  # Mutations
  
  mutation <- fread(paste0("benchmark_data/network_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
    arrange(gene_ID) %>%
    column_to_rownames("gene_ID")
  
  # CNV
  
  cnv <- fread(paste0("benchmark_data/network_",network_choice,"/cnv.csv"), select = c("gene_ID",samples)) %>%
    column_to_rownames("gene_ID")
  
}
 

#############################
# Filtering / Formatting
#############################

genes <- Reduce(intersect, list(
  rownames(l2fc_rna),
  rownames(mutation),
  rownames(cnv),
  append(network %>% pull(1), network %>% pull(2))
)) %>% sort

l2fc_rna <- l2fc_rna[genes, results_samples]
mutation <- mutation[genes, results_samples]
cnv <- cnv[genes, results_samples] 
  
#############################
# Create Temporary Files
#############################

write_tsv(l2fc_rna %>% rownames_to_column("Genes"), paste0("tmp/tmp_",cancer,"_OncoImpact_EXP.txt"))
write_tsv(mutation %>% rownames_to_column("Genes"), paste0("tmp/tmp_",cancer,"_OncoImpact_SNP.txt"))
write_tsv(network, paste0("tmp/tmp_",cancer,"_OncoImpact_network.txt") ,col_names = FALSE)
write_tsv(cnv %>% rownames_to_column("Genes"), paste0("tmp/tmp_",cancer,"_OncoImpact_cnv.txt"))



if(!dir.exists(paste0("results/",run_mode,"/network_",network_choice,"/OncoImpact"))){
  dir.create(paste0("results/",run_mode,"/network_",network_choice,"/OncoImpact/",cancer), recursive = T)
}



# Create Config File

perl_dir <- gsub("([A-Za-z]):", "/\\L\\1",WORK_DIR, perl = T)

if(!is.null(TMP_DIR)){
  out_dir <- gsub("([A-Za-z]):", "/\\L\\1",TMP_DIR, perl = T)  
}else{
  out_dir <- perl_dir
}

cfg_lines <- c(
  paste0("outDir=",out_dir,"/results/",run_mode,"/network_",network_choice,"/OncoImpact/",cancer),
  paste0("scriptDir=",perl_dir,"/scripts/OncoImpact"),
  paste0("numThreads=",threads),
  paste0("cnv=",perl_dir,"/tmp/tmp_",cancer,"_OncoImpact_cnv.txt"),
  paste0("exp=",perl_dir,"/tmp/tmp_",cancer,"_OncoImpact_EXP.txt"),
  paste0("snp=",perl_dir,"/tmp/tmp_",cancer,"_OncoImpact_SNP.txt"),
  "dataType=RNA_SEQ",
  "testMode=0",
  paste0("network=",perl_dir,"/tmp/tmp_",cancer,"_OncoImpact_network.txt")
  )

write_lines(cfg_lines, paste0("tmp/tmp_",cancer,"_OncoImpact_config.cfg"))
