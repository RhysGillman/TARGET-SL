#prepare_LOFGOF_annotations.R

# set working directory
#setwd("/home/workspace/files/TARGET_SL/")

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(parallel, quietly = T))

# Handling input arguments
option_list = list(
  #Sample info
  #make_option(c("-s", "--sampleinfo"), type="character", default="data/sample_info.csv", 
  #            help="path to sample info file", metavar ="Sample Info"),
  make_option(c("-l", "--LOFGOFfilepath"), type="character", default="data/LOFGOF/logofunc-predictions/LoGoFuncVotingEnsemble_preds_final.csv", 
              help="Path to LoGoFunc Prediction file", metavar ="LOFGOF File Path"),
  make_option(c("-m", "--mutationfilepath"), type="character", default="data/CCLE/OmicsSomaticMutations.csv", 
              help="Path to mutations file", metavar ="Mutations File Path"),
  make_option(c("-o", "--outputfilepath"), type="character", default="data/LOFGOF/itan1_annotated_mutations.csv", 
              help="Path to output file", metavar ="Output File Path"),
  make_option(c("-t", "--threads"), type="numeric", default=detectCores()-1, 
              help="Number of threads to use", metavar ="Threads"),
  make_option(c("-T", "--TempDir"), type="character", default="tmp", 
              help="Temporary directory to save tmp files", metavar ="TempDir"),
  make_option(c("-S", "--saveMemory"), type="logical", default=T, 
              help="Save memory by performing join in smaller bins and saving to temporary files", metavar ="saveMemory"),
  make_option(c("-L", "--Lines"), type="numeric", default=82468699, 
              help="Number of lines in the LOFGOF file", metavar ="Lines")
  
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

LOFGOFfilepath <- opt$LOFGOFfilepath
mutationfilepath <- opt$mutationfilepath
outputfilepath <- opt$outputfilepath
threads <- opt$threads
TempDir <- opt$TempDir
saveMemory <- opt$saveMemory
LOF_GOF_lines <- opt$Lines
#sample_info_path <- opt$sampleinfo

if(threads>1){
  cl <- makeCluster(threads, outfile = "log/map_genomic_LOF_GOFs.log")
  registerDoParallel(cl)
}


#############################
# Sample Info
#############################

#sample_info <- fread(sample_info_path) %>% filter(get_results)
#patients <- sample_info %>% pull(patient) %>% unique() %>% sort()
#mutations <- fread(mutationfilepath) %>% filter(patient %in% patients)


if(saveMemory){
  
  if(!dir.exists(TempDir)){
    dir.create(TempDir)
  }
  
  bins <- 20
  bin_size <- round(LOF_GOF_lines/bins)
  
  foreach(bin=1:bins, .packages = "tidyverse") %dopar% {
    #for(bin in 1:bins) {
    
    skip <- (bin-1)*bin_size
    
    if(bin!=bins){
      max_lines <- bin_size
    }else{
      max_lines <- Inf
    }
    
    
    
    annotations <- data.table::fread(LOFGOFfilepath, sep = "\t", header = T, skip = skip, nrows = max_lines, nThread = threads,
                                     col.names =  c("chromosome", "start_position", "reference_allele", "alt_allele", "LOF_GOF_id", "LOF_GOF_prediction", "score_neutral", "score_GOF", "score_LOF"),
                                     colClasses =  c("character","numeric","character","character","character","character","numeric","numeric","numeric")) %>%
      mutate(chromosome = ifelse(str_detect(chromosome,"chr"), chromosome,paste0("chr",chromosome)))
    
    annotated_mutations <- inner_join(mutations, annotations, by = c("chromosome","start_position","reference_allele", "alt_allele"))
    
    write_csv(annotated_mutations, paste0(TempDir,"/annotations_bin_",bin,".csv"))
    
    
  }
  
  
  annotated_mutations <- foreach(file=list.files(TempDir, pattern = "annotations_bin"), .combine = "rbind") %do% {
    fread(paste0(TempDir,"/",file))
  }
  
  final_annotated_mutations <- mutations %>%
    left_join(annotated_mutations)
  
  write_csv(final_annotated_mutations,outputfilepath)
  
  system(paste0("rm ", TempDir,"/annotations_bin_*"))
  
  
}else{
  
  LOF_GOF_predictions <- data.table::fread(LOFGOFfilepath, sep = "\t",
                                           col.names =  c("chromosome", "start_position", "reference_allele", "alt_allele", "LOF_GOF_id", "LOF_GOF_prediction", "score_neutral", "score_GOF", "score_LOF"),
                                           colClasses =  c("character","numeric","character","character","character","character","numeric","numeric","numeric")) %>%
    mutate(chromosome = ifelse(str_detect(chromosome,"chr"), chromosome,paste0("chr",chromosome)))
  
  
  annotated_mutations <- fread(mutationfilepath) %>%
    left_join(LOF_GOF_predictions, 
              by = c("chromosome","start_position","reference_allele", "alt_allele"))
  write_csv(annotated_mutations,outputfilepath)
  
}










# 82468699 lines