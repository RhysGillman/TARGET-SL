# run_sysSVM2.r.r

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-m", "--mode"), type="character", default="benchmark", 
              help="mode (benchmark or predict)", metavar ="Mode"),
  make_option(c("-a", "--annovarDir"), type="character", default="../annovar/", 
              help="the directory where annovar is installed", metavar ="Annovar Directory"),
  make_option(c("-s", "--sampleinfo"), type="character", default=NULL, 
              help="path to sample info file", metavar ="Sample Info"),
  make_option(c("-v", "--cnv"), type="character", default=NULL, 
              help="path to binary cnv file", metavar ="CNV"),
  make_option(c("-V", "--cnr"), type="character", default=NULL, 
              help="path to cony number ratio file", metavar ="CNR"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-c", "--cancertype"), type="character", default="Liver", 
              help="cancer type to analyse", metavar ="cancer Type")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

run_mode <- opt$mode
network_choice <- opt$network
cancer <- opt$cancertype
annovar_dir <- opt$annovarDir

if(run_mode=="predict"){
  sample_info_path <- opt$sampleinfo
  cnv_path <- opt$cnv
  cnr_path <- opt$cnr
  annovar_input_dir <- "data/ANNOVAR_input/"
}else if(run_mode=="benchmark"){
  cnv_path <- paste0("benchmark_data/network_",network_choice,"/cnv.csv")
  cnr_path <- paste0("benchmark_data/network_",network_choice,"/copy_numbers.csv")
  annovar_input_dir <- paste0("benchmark_data/network_",network_choice,"/ANNOVAR_input/")
}


#############################
# Functions
#############################

source("scripts/sysSVM2/annotation_functions_modified.R")
source("scripts/sysSVM2/train_predict_functions.R")


#############################
# ANNOVAR Setup
#############################


current_dir <- getwd()
setwd(annovar_dir)
if(!file.exists(paste0(annovar_dir,"/humandb/hg38_refGene.txt"))){
    message("hg38_refGene database missing")
    system("perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/")
}
if(!file.exists(paste0(annovar_dir,"/humandb/hg38_dbnsfp35a.txt"))){
    message("hg38_dbnsfp35a database missing")
    system("perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp35a humandb/")
}
if(!file.exists(paste0(annovar_dir,"/humandb/hg38_dbscsnv11.txt"))){
    message("hg38_dbscsnv11 database missing")
    system("perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbscsnv11 humandb/")
}
message("All required ANNOVAR databases available")
setwd(current_dir)


if(run_mode=="predict"){


  #############################
  # Sample Info
  #############################
  
  sample_info <- fread(sample_info_path) %>% filter(cancer_type==cancer, get_results)
  samples <- sample_info$sample_ID %>% unique() %>% sort()
  
  
  #############################
  # Entrez IDs
  #############################

  
  if(!file.exists("data/entrez_ids.csv")){
    suppressPackageStartupMessages (library(biomaRt, quietly = T))
    
    genes <- data.table::fread(cnv_path, select = c("gene_ID")) %>% pull(gene_ID) %>% unique()
    mart <- biomaRt::useMart("ensembl")
    mart <- useDataset("hsapiens_gene_ensembl",mart=mart)
    attr <- listAttributes(mart)
    entrez_id_mapping <- getBM(attributes=c('hgnc_symbol','entrezgene_id'), mart = mart, filters = 'hgnc_symbol', values = genes) %>%
      unique()
    write_csv(entrez_id_mapping, "data/entrez_ids.csv")
    detach(package:biomaRt)
    
  }else{
    entrez_id_mapping <- read_csv("data/entrez_ids.csv")
  }
  
  
  
}else if(run_mode=="benchmark"){
  
  #############################
  # Sample Info
  #############################
  
  sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv")) %>% filter(lineage==cancer)
  samples <- sample_info$cell_ID %>% sort()
  
  
  #############################
  # Entrez IDs
  #############################
  
  
  if(!file.exists(paste0("benchmark_data/network_", network_choice,"/entrez_ids.csv"))){
    suppressPackageStartupMessages (library(biomaRt, quietly = T))
    
    genes <- data.table::fread(paste0("benchmark_data/network_",network_choice,"/cnv.csv"), select = c("gene_ID")) %>% pull(gene_ID) %>% unique()
    mart <- biomaRt::useMart("ensembl")
    mart <- useDataset("hsapiens_gene_ensembl",mart=mart)
    attr <- listAttributes(mart)
    entrez_id_mapping <- getBM(attributes=c('hgnc_symbol','entrezgene_id'), mart = mart, filters = 'hgnc_symbol', values = genes)
    write_csv(entrez_id_mapping, paste0("benchmark_data/network_", network_choice,"/entrez_ids.csv"))
    detach(package:biomaRt)
    
  }else{
    entrez_id_mapping <- read_csv(paste0("benchmark_data/network_", network_choice,"/entrez_ids.csv"))
  }
  
  
  
}

################################
# Prepare Directories
################################

if(!dir.exists("tmp/annovar_tmp")){
  dir.create("tmp/annovar_tmp")
}


if(!dir.exists(paste0("results/",run_mode,"/network_",network_choice,"/sysSVM2"))){
  dir.create(paste0("results/",run_mode,"/network_",network_choice,"/sysSVM2/",cancer),recursive = T)
}

#########################################
# Get Pre-Trained PANCANCER sysSVM Model
#########################################

PANCAN_trained_sysSVM <- readRDS("scripts/sysSVM2/trained_models/PANCAN_trained_sysSVM.rds")

for(s in samples){

  #############################
  # Prepare Input Data
  #############################
  
  # Prepare input mutation file for ANNOVAR annotation
  avinput <- paste0(annovar_input_dir,cancer,"/",s,".avinput")
  
  # Prepare input CNV file for bedtools annotation
  #segment_cnv <- read_csv(paste0("validation_data/CCLE_",network_choice,"/cnv_segment.csv"), col_names = c("sample","chromosome","start","end","segment_mean"), col_types = "ccnnn") %>%
  #  filter(sample==s)
  
  #############################
  # Run ANNOVAR annotation
  #############################
  
  ssms_annotated <- annotate_ssms(
    avinput=avinput, 
    sample = s, 
    annovar_dir, 
    genome_version = "hg38", 
    gene_aliases_entrez = "scripts/sysSVM2/annotation_reference_files/gene_aliases_entrez.tsv", 
    hotspots = "scripts/sysSVM2/annotation_reference_files/tcga_pancancer_hotspots_oncodriveclust.tsv",
    temp_dir = "tmp/annovar_tmp/"
  )
  
  #############################
  # Run Bedtools annotation
  #############################

  #cnvs_annotated <- annotate_cnvs(
  #    cnv_segments=segment_cnv,            # Table with the following columns: sample; chromosome; start; end; and copy_number or segment_mean
  #    ploidy = NULL,           # Table with the following columns: sample; ploidy. Leave null if unavailable (assumes diploidy)
  #    ploidy_threshold = 2,    # Threshold for determining gene amplifications: CN >= ploidy_threshold*ploidy
  #    gene_coords = "scripts/sysSVM2/annotation_reference_files/gene_coords_hg38.tsv",             # gene_coords_hg19.tsv or gene_coords_hg38.tsv from the sysSVM2 GitHub repository
  #    bedtools_bin_dir = bedtools_dir, # Directory where the bedtools binary executable is located, if not in $PATH
  #    temp_dir = "tmp/bedtools_tmp"     # Directory for temporary files to be created
  #)
  
  # No longer running bedtools due to inconsistencies between algorithms with calculating CNV. Therefore, instead using the already prepared CNV information
  
  
  if(run_mode=="predict"){
    cnv_change <- fread(cnv_path, select = c("gene_ID",s)) %>%
      dplyr::rename(cnv_change=2) %>%
      mutate(sample = s, 
             CNVGain = cnv_change > 0, 
             CNVLoss = cnv_change < 0)
    cnv_number <- fread(cnr_path, select = c("gene_ID",s)) %>%
      dplyr::rename(cnr=2) %>%
      mutate(sample = s, Copy_number=round(2^cnr*2))
    cnv <- inner_join(cnv_change, cnv_number, by = c("sample", "gene_ID"), ) %>%
      inner_join(entrez_id_mapping, by = c("gene_ID"="hgnc_symbol"), multiple="all") %>%
      # Converting CN ratio to copy number
      #mutate(Copy_number = ifelse(cnr >= 1 & cnr < 2, 2,
      #                            ifelse(cnr >= 3 & cnr < 4, 3, round(2*cnr,0)))) %>%
      # Remove any duplicate entrez ids by prioritising the version with a cnv change
      group_by(sample,entrezgene_id) %>%
      arrange(cnv_change == 0) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(sample, entrez=entrezgene_id, Copy_number, CNVGain, CNVLoss) %>%
      mutate(Copy_number = as.integer(Copy_number),
             CNVGain = as.integer(CNVGain),
             CNVLoss = as.integer(CNVLoss)
      )
  }else if(run_mode=="benchmark"){
    cnv_change <- fread(cnv_path, select = c("gene_ID",s)) %>%
      dplyr::rename(cnv_change=2) %>%
      mutate(sample = s, 
             CNVGain = cnv_change > 0, 
             CNVLoss = cnv_change < 0)
    cnv_number <- fread(cnr_path, select = c("gene_ID",s)) %>%
      dplyr::rename(Copy_number=2) %>%
      mutate(sample = s)
    cnv <- inner_join(cnv_change, cnv_number, by = c("sample", "gene_ID"), ) %>%
      inner_join(entrez_id_mapping, by = c("gene_ID"="hgnc_symbol"), multiple="all") %>%
      # Converting CN ratio to copy number
      mutate(Copy_number = ifelse(Copy_number >= 1 & Copy_number < 2, 2,
                                  ifelse(Copy_number >= 3 & Copy_number < 4, 3, round(2*Copy_number,0)))) %>%
      # Remove any duplicate entrez ids by prioritising the version with a cnv change
      group_by(sample,entrezgene_id) %>%
      arrange(cnv_change == 0) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(sample, entrez=entrezgene_id, Copy_number, CNVGain, CNVLoss) %>%
      mutate(Copy_number = as.integer(Copy_number),
             CNVGain = as.integer(CNVGain),
             CNVLoss = as.integer(CNVLoss)
      )
  }
  


  
  
  #############################
  # Combine Annotations
  #############################
  
  #molecular_data <- make_totalTable(
  #  ssms_annotated, cnvs_annotated, 
  #  canonical_drivers = "scripts/sysSVM2/example_data/canonical_drivers.rds"
  #)
  
  
  canonical_drivers <- readRDS("scripts/sysSVM2/example_data/canonical_drivers.rds")
  canonical_drivers <- canonical_drivers %>% select(entrez, geneType)
  
  # Summarise mutations at the level of sample-gene pairs
  sample_gene_muts = ssms_annotated %>%
    group_by(sample, entrez) %>%
    summarise(
      no_ALL_muts = n(), 
      no_TRUNC_muts = sum(TRUNC_mut), 
      no_NTDam_muts = sum(NTDam_mut), 
      no_GOF_muts = sum(GOF_mut)
    ) %>%
    ungroup

  
  # Join
  totalTable = full_join(
    sample_gene_muts,
    cnv,
    by = c("sample", "entrez")
  ) %>%
    replace_na(list(
      no_ALL_muts = 0, no_TRUNC_muts = 0, no_NTDam_muts = 0, no_GOF_muts = 0,
      Copy_number = 2, CNVGain = 0, CNVLoss = 0
    ))
  
  
  # Subset for damaged genes (including GoF in oncogenes and LoF in TSGs)
  totalTable = totalTable %>% 
    subset(no_TRUNC_muts + no_NTDam_muts + no_GOF_muts > 0 | Copy_number == 0 | CNVGain == 1) %>%
    left_join(canonical_drivers, by = "entrez") %>%
    subset(
      is.na(geneType) |
        (geneType == "oncogene" & no_NTDam_muts + no_GOF_muts >= 1 | CNVGain == 1) |
        (geneType == "tumourSuppressor" & no_TRUNC_muts + no_NTDam_muts >= 1 | Copy_number == 0) |
        (geneType == "TP53" & no_TRUNC_muts + no_NTDam_muts + no_GOF_muts >= 1 | Copy_number == 0)
    )
  
  
  # Output
  totalTable = totalTable %>% dplyr::select(-geneType)
  
  
  
  
  
  #############################
  # Run predictions Annotations
  #############################
  
  predictions <- predict_sysSVM2(
    trained_sysSVM = PANCAN_trained_sysSVM, 
    molecular_data = totalTable, 
    systemsLevel_data = "scripts/sysSVM2/example_data/systemsLevel_features_allGenes.tsv"
  )
  
  
  write_csv(predictions, paste0("results/",run_mode,"/network_",network_choice,"/sysSVM2/",cancer,"/",s,".csv"))

}