#format_PersonaDrive_results.R
# This code reformats the output from PersonaDrive to be consistent with the other algorithms

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
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-c", "--cancertype"), type="character", default="Liver", 
              help="cancer type to analyse", metavar ="Cell Type")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
cancer <- opt$cancertype
run_mode <- opt$mode

if(run_mode=="predict"){
  sample_info_path <- opt$sampleinfo
  sample_info <- fread(sample_info_path) %>% filter(cancer_type == cancer, get_results)
  samples <- sample_info$patient %>% unique()
  
}else{
  sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv")) %>% filter(lineage==cancer)
  samples <- sample_info$cell_ID %>% sort()
}





suppressWarnings(rm(PersonaDrive_results))

results_length <- length(read_lines(paste0("results/",run_mode,"/network_",network_choice,"/PersonaDrive/",cancer, "/PersonaDrive.txt")))

for(i in 1:results_length){
  #message(i)
  sample_drivers <- fread(paste0("results/",run_mode,"/network_",network_choice,"/PersonaDrive/",cancer, "/PersonaDrive.txt"), header = FALSE, skip = (i-1), nrows = 1)
  indiv_sample <- sample_drivers[1,1] %>% as.character()
  if(!indiv_sample %in% gsub("_","-",samples)) {
    next
  }else{
    
    colnames(sample_drivers) <- c("sample_ID", 1:(ncol(sample_drivers)-1))
    sample_drivers <- sample_drivers %>%
      gather("rank", "driver",-sample_ID) %>%
      na.omit()
    if(exists("PersonaDrive_results")){
      PersonaDrive_results <- rbind(PersonaDrive_results, sample_drivers)
    }else{
      
      PersonaDrive_results <- sample_drivers
    }
    
  }
  
}

PersonaDrive_results <- PersonaDrive_results %>%
  mutate(cancer_type=cancer) %>%
  dplyr::select(cancer_type, sample_ID, driver, rank) %>%
  mutate(sample_ID=gsub("-","_",sample_ID))

write_csv(PersonaDrive_results, paste0("results/",run_mode,"/network_",network_choice,"/PersonaDrive/",cancer, ".csv"))


