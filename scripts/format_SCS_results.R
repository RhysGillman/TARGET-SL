#format_SCS_results.R
# This code reformats the output from SCS to be consistent with the other algorithms

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
              help="cancer type to analyse", metavar ="cancer Type")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

run_mode <- opt$mode
network_choice <- opt$network
cancer_type <- opt$cancertype

if(run_mode=="benchmark"){

suppressWarnings(rm(SCS_results))
results_names <- read_csv(paste0("results/",run_mode,"/network_",network_choice,"/SCS/",cancer_type,"/sample_names.txt")) %>% deframe()
for(result in names(results_names)){
  path <- paste0("results/",run_mode,"/network_",network_choice,"/SCS/",cancer_type,"/result_sample_",result,".csv")
  if(file.exists(path)){
    tmp <- read_csv(path, skip = 1, col_names = c("gene_ID", "module", "impact"))
    
    # The below is required for some broken outputs that do not quote the module vector
    
    if(ncol(tmp) > 3){
      tmp <- tmp %>%
        unite("module_fixed", 2:(ncol(tmp)-1), sep = ",", na.rm=T)
      colnames(tmp) <- c("gene_ID", "module", "impact")
    }
    
    tmp <- tmp %>%
      arrange(desc(impact)) %>%
      mutate(rank = row_number(), sample_ID=results_names[result])
    if(!exists("SCS_results")){
      SCS_results <- tmp
    }else{
      SCS_results <- rbind(SCS_results, tmp)
    }
  }
  
}

SCS_results <- SCS_results %>%
  dplyr::mutate(cancer_type) %>%
  arrange(sample_ID) %>%
  dplyr::select(cancer_type, sample_ID, driver = gene_ID, rank)

write_csv(SCS_results, paste0("results/",run_mode,"/network_",network_choice,"/SCS/",cancer_type,".csv"))

}