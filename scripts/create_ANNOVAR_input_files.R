#create_ANNOVAR_input_files.R
# This code converts the CCLE variant calls into the input format required by ANNOVAR

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
  make_option(c("-d", "--dna"), type="character", default=NULL, 
              help="path to mutations file", metavar ="Network"),
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
  dna_path <- opt$dna
  
  if(!dir.exists(paste0("data/ANNOVAR_input/",cancer))){
    dir.create(paste0("data/ANNOVAR_input/",cancer), recursive = T)
  }
  
}

if(run_mode=="predict"){
  #############################
  # Sample Info
  #############################
  
  sample_info <- fread(sample_info_path) %>% filter(cancer_type==cancer, get_results)
  samples <- sample_info$sample_ID %>% unique() %>% sort()
  
  #############################
  # Create ANNOVAR Input Files
  #############################
  
  # Annovar format Chromosome ("chr" prefix is optional), Start, End, Reference Allele, Alternative Allele
  # tab separated
  # https://annovar.openbioinformatics.org/en/latest/user-guide/input/
  
  mutations <- fread(dna_path)
  
  for(sample in samples){
    annovar_input <- mutations %>%
      filter(patient == sample) %>%
      dplyr::select(chromosome,start_position,end_position,reference_allele,alt_allele)
    write_tsv(annovar_input, col_names = F, paste0("data/ANNOVAR_input/",network_choice,"/",cancer,"/",sample,".avinput"))
    
  }
}else if(run_mode=="benchmark"){
  
  #############################
  # Sample Info
  #############################
  
  sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv")) %>% filter(lineage==cancer)
  samples <- sample_info$cell_ID %>% sort()
  
  #############################
  # Create ANNOVAR Input Files
  #############################
  
  # Annovar format Chromosome ("chr" prefix is optional), Start, End, Reference Allele, Alternative Allele
  # tab separated
  # https://annovar.openbioinformatics.org/en/latest/user-guide/input/
  
  mutations <- fread(paste0("benchmark_data/network_",network_choice,"/mutations_MAF.csv"), sep = ",")
  
  for(sample in samples){
    annovar_input <- mutations %>%
      filter(cell_ID == sample) %>%
      mutate(End = ifelse(nchar(Ref) > 2 & Alt=="-", Pos + nchar(Ref) - 1, Pos)) %>%
      dplyr::select(Chrom,Start=Pos,End,Ref,Alt)
    write_tsv(annovar_input, col_names = F, paste0("benchmark_data/network_", network_choice, "/ANNOVAR_input/",cell_type,"/",sample,".avinput"))
    
  }
  
}