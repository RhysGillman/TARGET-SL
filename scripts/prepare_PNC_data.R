# prepare_PNC_data.R
# This code prepares the input data to run PNC and stores it in tmp/

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
              help="cell type to analyse", metavar ="Cell Type")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

run_mode <- opt$mode
network_choice <- opt$network
cancer_type <- opt$cancertype

if(run_mode=="benchmark"){
  


#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv")) %>% filter(lineage==cancer_type)
samples <- sample_info$cell_ID %>% sort()

#############################
# Prepare Input Data
#############################

# RNA
# PNC requires paired tumour and normal expression data
# Because this is not possible with the cell line data, we will make pseudonormal samples based on the expression of all 'other' cells iteratively



rna <- fread(paste0("benchmark_data/network_",network_choice,"/tpm.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID") 

suppressWarnings(rm(pseudonormal_rna))
for(sample in samples){
  # Collect expression data for all samples other than the current sample
  other_rna <- rna[,samples[samples != sample]]
  # Collect mean expression for other samples as the pseudo-normal expression for each cell
  tmp_pseudonormal_rna <- enframe(rowMeans(other_rna), name = "gene_ID", value = sample) %>% column_to_rownames("gene_ID")
  
  
  if(!exists("pseudonormal_rna")){
    pseudonormal_rna <- tmp_pseudonormal_rna
  }else{
    pseudonormal_rna <- cbind(pseudonormal_rna, tmp_pseudonormal_rna)
  }
  
}

rm(tmp_pseudonormal_rna, other_rna)

# Formatting

genes <- rownames(rna) %>% unique() %>% sort()

rna <- rna[genes, samples]
pseudonormal_rna <- pseudonormal_rna[genes, samples]



#############################
# Create Temporary Files
#############################

write_tsv(rna %>% rownames_to_column("Genes"), "tmp/tmp_PNC_tumour_expression.txt")
write_tsv(pseudonormal_rna %>% rownames_to_column("Genes"), "tmp/tmp_PNC_pseudonormal_expression.txt")



if(!dir.exists(paste0("results/",run_mode,"/network_",network_choice,"/PNC"))){
  dir.create(paste0("results/",run_mode,"/network_",network_choice,"/PNC"))
}
}
