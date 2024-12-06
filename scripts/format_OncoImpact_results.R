#format_OncoImpact_results.R
# This code reformats the output from OncoImpact to be consistent with the other algorithms

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
              help="cancer type to analyse", metavar ="CancerType"),
  make_option(c("-T", "--tmp"), type="character", default=NULL, 
              help="temporary directory with symbolic link permission", metavar ="Tmp")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

run_mode <- opt$mode
network_choice <- opt$network
cancer <- opt$cancertype
TMP_DIR <- opt$tmp

if(TMP_DIR==""){
  TMP_DIR=NULL
}

if(!is.null(TMP_DIR)){
  out_dir <- TMP_DIR 
}else{
  out_dir <- "./"
}


suppressWarnings(rm(OncoImpact_results))

for(result_file in list.files(paste0(out_dir,"/results/",run_mode,"/network_",network_choice,"/OncoImpact/",cancer,"/sample_driver_list/"), pattern = ".txt$")){
  sample_ID <- gsub("\\.txt", "", result_file)
  # Read in and fix first line
  result1 <- fread(paste0(out_dir,"/results/",run_mode,"/network_",network_choice,"/OncoImpact/",cancer,"/sample_driver_list/", result_file), 
                   nrows = 1, 
                   header = F)
  
  if(ncol(result1)>8){
    result1 <- result1 %>%
      dplyr::select(gene_ID=8, impact=10) %>% 
      mutate(gene_ID=gsub("PAN_CANCER","",gene_ID))  %>%
      mutate(sample_ID) %>%
      relocate(sample_ID)
    result2 <- fread(paste0(out_dir,"/results/",run_mode,"/network_",network_choice,"/OncoImpact/",cancer,"/sample_driver_list/", result_file),
                     select = c(1,3),fill = T,header = T) %>%
      dplyr::select(gene_ID = `#GENE`, impact = SAMPLE_IMPACT) %>%
      mutate(sample_ID) %>%
      relocate(sample_ID)
    result <- rbind(result1,result2)
  }else{
    next
  }
  
  if(!exists("OncoImpact_results")){
    OncoImpact_results <- result
  }else{
    OncoImpact_results <- rbind(OncoImpact_results, result)
  }
}

OncoImpact_results <- OncoImpact_results %>%
  arrange(desc(impact)) %>%
  group_by(sample_ID) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(-impact) %>%
  arrange(sample_ID) %>%
  dplyr::mutate(cancer_type = cancer) %>%
  dplyr::select(cancer_type, sample_ID, driver = gene_ID, rank)



write_csv(OncoImpact_results, paste0("results/",run_mode,"/network_",network_choice,"/OncoImpact/",cancer,".csv"))
