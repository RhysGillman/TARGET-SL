#prepare_PersonaDrive_data.R
# This data prepares the input data to run PersonaDrive and stores it in tmp/

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

# Using code from https://github.com/shahcompbio/drivernet/blob/master/R/getPatientOutlierMatrix.R as the authors of PersonaDrive said they did.

# Personadrive uses a threshold of 0.5 instead of 2

getPatientOutlierMatrix <- function(patExpMatrix, th=2)
  {
  # Get sd of each column
  expSd   <- apply(patExpMatrix, 2, sd)
  # Only keep columns with sd > 0
  id <- expSd > 0
  expSd <- expSd[id]
  patExpMatrix <- patExpMatrix[, id]
  # Get the mean of each column
  expMean <- apply(patExpMatrix, 2, mean)
	
  # dnorma returns the probability density function for a normal distribution given the mean, sd
  # I believe this is the y axis value on a normal distribution
  # Overall, this part converts the values to their probability density value given the sd and mean of the column
  # Thus, higher value means closer to mean
  num <- dnorm(x=t(patExpMatrix), mean=expMean, sd=expSd, log=T)
  num <- t(num)
  # Note, adding a vector to a matrix adds each value to the entire row
  # Thus, if columns are genes, getting the mean + 0.5sd of each gene and getting the threshold probability density value corresponding to this
  numSd <- dnorm( x=t(expMean+th*expSd), mean=expMean, sd=expSd, log=T )
  
  # Creating a matrix that repeats the threshold values per gene for every sample
  y <- rep(numSd, each=dim(patExpMatrix)[1])
  y <- matrix(y, nrow=dim(patExpMatrix)[1], ncol=dim(patExpMatrix)[2])
	
  # TRUE values are outliers because they have probability densitites lower than the threshold
  patOutMatrix <- num <= y
  
  id <- colSums(patOutMatrix) > 0
  patOutMatrix <- patOutMatrix[, id]
  
  return(patOutMatrix)
}


#############################
# Network
#############################


#The default network used by PersonaDrive is weighted and seems to be non-directed, though the authors don't state this.
#It is in the form on a long table (3 columns protein 1, protein 2, weight.) It is a weighted network with weight values 
#ranging from 0.701 to 0.999. Every row in the table is unique, however every single edge is duplicated in reverse 
#(gene1-gene2 is also gene2-gene1), and thus every single gene is present in both columns of the table.

if(network_choice=="own"){
  network <- fread("data/own_networks/PRODIGY_personadrive.csv", col.names = c("genes1","genes2","Score"))
}else{
  network <- fread(paste0("benchmark_data/network_",network_choice,"/network_undirected.csv"))
  colnames(network) <- c("protein_1", "protein_2", "confidence")
  reverse_edges <- network %>% dplyr::select(protein_1 = protein_2, protein_2 = protein_1, confidence)
  network <- rbind(network, reverse_edges) %>% unique()
  rm(reverse_edges)
  network <- network %>% mutate(confidence = confidence/1000)
  colnames(network) <- c("genes1","genes2","Score")
}


if(run_mode=="predict"){

  #############################
  # Sample Info
  #############################
  
  sample_info <- read_csv(sample_info_path) %>% filter(cancer_type==cancer, !is.na(tumour_sample_ID))
  samples <- sample_info$sample_ID %>% unique() %>% sort()
  rna_samples <- sample_info %>% column_to_rownames("sample_ID")
  rna_samples <- rna_samples[samples,]
  rna_samples <- rna_samples$tumour_sample_ID
  
  #############################
  # Prepare Input Data
  #############################
  
  # RNA
  
  rna <- fread(rna_path) %>%
    arrange(gene_ID) %>%
    data.table::transpose(make.names = "gene_ID", keep.names = "tumour_sample_ID") %>%
    inner_join(sample_info %>% dplyr::select(tumour_sample_ID,sample_ID), by = "sample_ID") %>%
    dplyr::select(-tumour_sample_ID) %>%
    column_to_rownames("sample_ID")
  
  # Mutations
  
  mutation <- fread(dna_path) %>% 
    arrange(gene_ID) %>%
    column_to_rownames("gene_ID")
  
  # Filtering Patients
  
  avail_samples <- intersect(rownames(rna), colnames(mutation))
  avail_samples <- intersect(avail_samples, samples) %>% sort()
    
  
  rna <- rna[avail_samples,]
  mutation <- mutation[,avail_samples]
    
  
  rna_outliers <- getPatientOutlierMatrix(rna, th = 0.5)
  
}else if(run_mode == "benchmark"){
  
  #############################
  # Sample Info
  #############################
  
  sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv")) %>% filter(lineage==cancer)
  samples <- sample_info$cell_ID %>% sort()
  
  #############################
  # Prepare Input Data
  #############################
  
  # RNA
  
  rna <- fread(paste0("benchmark_data/network_",network_choice,"/tpm.csv"), select = c("gene_ID",samples)) %>% 
    arrange(gene_ID) %>%
    data.table::transpose(make.names = "gene_ID", keep.names = "sample_ID") %>%
    column_to_rownames("sample_ID")
  
  rna_outliers <- getPatientOutlierMatrix(rna, th = 0.5)
  
  # Mutations
  
  mutation <- fread(paste0("benchmark_data/network_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
    arrange(gene_ID) %>%
    column_to_rownames("gene_ID")
  
}


# Formatting

genes <- Reduce(intersect, list(
  colnames(rna_outliers),
  rownames(mutation),
  append(network %>% pull(1), network %>% pull(2))
)) %>% sort()

rna_outliers <- rna_outliers[,genes]
mutation <- mutation[genes,]

#underscores are problematic for PersonaDrive
fixed_sample_IDs <- gsub("_","-",rownames(rna_outliers))

if(all(rownames(rna_outliers)==colnames(mutation))){
  rownames(rna_outliers) <- fixed_sample_IDs
  colnames(mutation) <- fixed_sample_IDs
}else{
  warning("Problem with sample or patient IDs, please fix!")
}


#############################
# Create Temporary Files
#############################


write.table(rna_outliers, "tmp/tmp_PersonaDrive_DEGs.csv", sep = "," )
write_tsv(network, "tmp/tmp_PersonaDrive_network.tsv")
write.csv(mutation, "tmp/tmp_PersonaDrive_MUT.csv", quote = FALSE)


if(!dir.exists(paste0("results/",run_mode,"/network_",network_choice,"/PersonaDrive/",cancer))){
  dir.create(paste0("results/",run_mode,"/network_",network_choice,"/PersonaDrive/",cancer), recursive = T)
}



