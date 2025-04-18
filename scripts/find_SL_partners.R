#find_SL_partners.R
# This code finds SL partners for CCLE mutations

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(reshape2, quietly = T))
suppressPackageStartupMessages (library(readxl, quietly = T))
suppressPackageStartupMessages (library(biomaRt, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-m", "--max"), type="numeric", default=5, 
              help="Maximum number of SL-partners to find for each gene", metavar ="max"),
  make_option(c("-r", "--rna"), type="character", default="data/merged_vst_rna.csv", 
              help="path to rna file", metavar ="Sample Info"),
  make_option(c("-n", "--network"), type="character", default="own", 
              help="network being used", metavar ="Network")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

max_SL <- opt$max
rna_path <- opt$rna
network_choice <- opt$network


######################
# Info on sources
######################

# SynLethDB
# Combination of experimental screens (RNAi, CRISPR, drugs), in silico screens, text mining and scouring other databases
# Scores have confidence values based on quality of evidence
# v2 released 2022

# BioGRID
# Open access repository of CRISPR screens
# Includes multiple levels of synthetic lethality, also dosage lethality
# Database has been pre-filtered and aggregated, code in Bioinformatics/BioGRID_SL/filtering.rmd
# Release 4.4.224, 25/7/23

# Berrena_et_al_gMCSs
# Predicted minimal cut sets based on metabolic and functional genetic pathways
# Filtered to only size=2 MCSs and given confidence scores based on PPV of each model
# Last updated 05/23


# Slorth
# Random-forest classifier based on PPI network features and GO terms

# SLOAD
# Random forest predictions
# Trained on experimentally validated SLs features (mutation, CN, methylation, mRNA)
# Classified against candidate SLs from DAISY and MiSL
# Cancer specific

# CGIdb
# Based on mutual exclusivity of mutations in gene-pairs
# Based on data from TCGA, CCLE and GDSC
# Last updated 2018-12-18

# SiLi
# Predictions in liver cancer only based on Stuart's rank aggregation of (1) function similarity (2) differential expression (3) correlation (4) survival hazard ratio

######################
# Genes to consider
######################

genes <- fread(rna_path,sep = ",", select = c("gene_ID")) %>%
  pull(gene_ID)



######################
# Prepare gMCSs
######################

if(!file.exists("data/SL/Berrena_et_al_gMCSs/gMCS_SL_pairs.csv")){
  
  # To score the gMCS SL-predictions, predictions were scored between 0 and 1 based on the average positive predictive value of the model in the original paper's analysis
  
  
  # Get names of various models used for predictions
  
  names <- c("Human1",
             "Human1-O1", "Human1-O2", "Human1-O3",
             "Human1-D1", "Human1-D2", "Human1-D3",
             "Human1-T1", "Human1-T2", "Human1-T3")
  
  # Only using the data capped at gMCS length = 5
  
  length <- 5
  
  # First using the papers evaluation results. Using PPV to rank the models to give confidence values
  # Read in all results and rbind
  
  HartResultList <- lapply(names, function(x){do.call(rbind,readRDS(paste0("data/SL/Berrena_et_al_gMCSs/Hart_AllThresholdsResults_",  x, "_length_", length, ".RDS")))})
  AllThresholdAnalysis <- do.call(rbind,HartResultList) %>% mutate(source="Hart")
  DepMapResultList <- lapply(names, function(x){do.call(rbind,readRDS(paste0("data/SL/Berrena_et_al_gMCSs/DepMap_AllThresholdsResults_",  x, "_length_", length, ".RDS")))})
  AllThresholdAnalysis <- rbind(AllThresholdAnalysis,
                                do.call(rbind,DepMapResultList) %>% mutate(source="DepMap"))
  
  Melted_gMCS_All_Statistics_th <- melt(data = AllThresholdAnalysis,
                                        id.vars = c("source","Order", "cellLine"),
                                        measure.vars = c("True Positives", "False Positives",
                                                         "False Negatives", "True Negatives",
                                                         "Accuracy", "Sensitivity", "Specificity",
                                                         "Positive Predicted Value", "Matthew's Cor. Coef."))
  # Get PPVs
  
  PPV <- Melted_gMCS_All_Statistics_th %>%
    filter(variable=="Positive Predicted Value") %>%
    group_by(Order) %>%
    summarise(mean_PPV=mean(value))
  
  maxPPV <- max(PPV$mean_PPV)
  minPPV <- min(PPV$mean_PPV)
  
  # Scale models from 0-1 based on average PPV across all cell lines
  
  ranking <- PPV %>%
    mutate(ranking=(mean_PPV-minPPV)/(maxPPV-minPPV)) %>%
    dplyr::select(source=Order, ranking)
  
  # Read in predicted gMCSs for all models
  
  gMCS_list <- lapply(names, function(x){read_xlsx("data/SL/Berrena_et_al_gMCSs/gMCSs_list.xlsx", sheet = x, col_types = "text") %>% as.data.frame() %>% mutate(source = x)})
  gMCS_SL_pairs <- do.call(rbind,gMCS_list) %>%
    # Only keep length = 2 gMCSs for SL pairs
    filter(rowSums(!is.na(.[1:5]))==2) %>%
    # Add confidence scores
    left_join(ranking, by = "source") %>%
    dplyr::select(gene1=V1,gene2=V2,confidence=ranking) %>%
    group_by(gene1,gene2) %>%
    summarise(confidence=max(confidence)) %>%
    separate_longer_delim(gene1,"_") %>%
    separate_longer_delim(gene1,"_")
  
  # Using BioMart to convert ensembl IDs to HGNC IDs
  
  #mart <- useMart("ENSEMBL_MART_ENSEMBL")
  #mart <- useDataset('hsapiens_gene_ensembl', mart)
  #gene_ids <- getBM(
  #  mart = mart,
  #  attributes = c('ensembl_gene_id','hgnc_symbol'),
  #  filter = 'ensembl_gene_id',
  #  values = unique(append(gMCS_SL_pairs$gene1,gMCS_SL_pairs$gene2)),
  #  uniqueRows = TRUE)
  
  gene_ids <- fread("data/gene_id_mapping.tsv") %>% dplyr::select(ensembl_gene_id, hgnc_symbol)
  
  gMCS_SL_pairs <- gMCS_SL_pairs %>%
    left_join(gene_ids, by = c("gene1"="ensembl_gene_id")) %>%
    rename(gene1_hgnc=hgnc_symbol) %>%
    mutate(gene1_hgnc=ifelse(!str_detect(gene1,"^ENSG"),gene1,gene1_hgnc)) %>%
    left_join(gene_ids, by = c("gene2"="ensembl_gene_id")) %>%
    rename(gene2_hgnc=hgnc_symbol) %>%
    mutate(gene2_hgnc=ifelse(!str_detect(gene2,"^ENSG"),gene2,gene2_hgnc)) %>%
    dplyr::select(gene1=gene1_hgnc,gene2=gene2_hgnc,confidence) %>%
    group_by(gene1,gene2) %>%
    summarise(confidence=max(confidence))
  
  write_csv(gMCS_SL_pairs,"data/SL/Berrena_et_al_gMCSs/gMCS_SL_pairs.csv")
  
  suppressWarnings(rm(gMCS_list,AllThresholdAnalysis,DepMapResultList,HartResultList,gene_ids,mart,Melted_gMCS_All_Statistics_th,PPV,ranking,length,maxPPV,minPPV,names))
  
}else{
  gMCS_SL_pairs <- read_csv("data/SL/Berrena_et_al_gMCSs/gMCS_SL_pairs.csv")
}

gMCS_SL_pairs <- gMCS_SL_pairs %>% filter(gene1 %in% genes & gene2 %in% genes)

######################
# Prepare bioGRID
######################

if(!file.exists("data/SL/BioGRID/SSL_filtered.csv")){
  
  # Read in entire BioGRID database, keeping only relevant columns
  BioGRID <- data.table::fread("data/SL/BioGRID/BIOGRID-ALL-4.4.224.tab3.txt", sep = "\t", select = c(
    "Official Symbol Interactor A","Official Symbol Interactor B","Synonyms Interactor A","Synonyms Interactor B","Experimental System","Throughput", "Organism Name Interactor A", "Organism Name Interactor B"
  )) %>%
    # Filter to human data
    filter(`Organism Name Interactor A` == "Homo sapiens" & `Organism Name Interactor B` == "Homo sapiens") %>%
    # Filter to data related to synthetic lethality and synthetic doasge lethality
    filter(`Experimental System` %in% c("Synthetic Lethality","Synthetic Growth Defect","Dosage Lethality","Dosage Growth Defect","Negative Genetic")) %>%
    # Removing incompatible gene IDs and expanding synonyms
    mutate(`Synonyms Interactor A`=paste0(`Synonyms Interactor A`,"|-")) %>%
    separate_longer_delim(cols = "Synonyms Interactor A", delim = "|") %>%
    filter(!str_detect(`Synonyms Interactor A`,"[a-z]| |L00|.+\\-|S00|\\+|\\(|\\.")) %>%
    mutate(`Synonyms Interactor B`=paste0(`Synonyms Interactor B`,"|-")) %>%
    separate_longer_delim(cols = "Synonyms Interactor B", delim = "|") %>%
    filter(!str_detect(`Synonyms Interactor B`,"[a-z]| |L00|.+\\-|S00|\\+|\\(|\\.")) %>%
    mutate(`Official Symbol Interactor A`=ifelse(`Synonyms Interactor A`!="-",`Synonyms Interactor A`,`Official Symbol Interactor A`)) %>%
    mutate(`Official Symbol Interactor B`=ifelse(`Synonyms Interactor B`!="-",`Synonyms Interactor B`,`Official Symbol Interactor B`)) %>%
    filter(!str_detect(`Official Symbol Interactor A`,"[a-z]| |L00|.+\\-|S00|\\+|\\(|\\.")) %>%
    filter(!str_detect(`Official Symbol Interactor B`,"[a-z]| |L00|.+\\-|S00|\\+|\\(|\\.")) %>%
    # Initially scoring based on throughput, high throughput = 0, low throughput = 1
    mutate(score=ifelse(str_detect(Throughput,"Low Throughput"), 1, 0)) %>%
    dplyr::select(geneA=`Official Symbol Interactor A`,geneB=`Official Symbol Interactor B`, type=`Experimental System`,score) %>%
    group_by(geneA,geneB,type) %>%
    summarise(score=max(score)) %>%
    # Additionally scoring based on strength of relationship
    # Best (+1)= Synthetic Lethality, Dosage Lethality
    # Middle (+0.5)= Synthetic Growth Defect, Dosage Growth Defect
    # Worst (+0)= Negative Genetic
    mutate(score=ifelse(type %in% c("Synthetic Lethality","Dosage Lethality"), score+1,score)) %>%
    mutate(score=ifelse(type %in% c("Synthetic Growth Defect","Dosage Growth Defect"), score+0.5,score)) %>%
    mutate(type=ifelse(type %in% c("Synthetic Lethality","Synthetic Growth Defect","Negative Genetic"), "SSL",type)) %>%
    mutate(type=ifelse(type %in% c("Dosage Lethality","Dosage Growth Defect"), "SDL",type)) %>%
    group_by(geneA,geneB,type) %>%
    summarise(score=max(score)) %>%
    mutate(score=(score)/2) %>%
    filter(geneA!=geneB) %>%
    dplyr::select(gene1=geneA,gene2=geneB,confidence=score,type)
  
  BioGRID_SSL <- BioGRID %>% filter(type=="SSL") %>% dplyr::select(-type)
  BioGRID_SDL <- BioGRID %>% filter(type=="SDL") %>% dplyr::select(-type)
  
  write_csv(BioGRID_SSL, "data/SL/BioGRID/SSL_filtered.csv")
  write_csv(BioGRID_SDL, "data/SL/BioGRID/SDL_filtered.csv")
  
  suppressWarnings(rm(BioGRID,BioGRID_SDL))
  
  
}else{
  BioGRID_SSL <- read_csv("data/SL/BioGRID/SSL_filtered.csv")
}

BioGRID_SSL <- BioGRID_SSL %>% filter(gene1 %in% genes & gene2 %in% genes)


######################
# Prepare CGIdb
######################


# This version is downloaded from the website http://www.medsysbio.org/CGIdb/data/all_SL/   v0.1
CGIdb_old <- read_csv("data/SL/CGIdb/All_synthetic_lethality_gene_pairs.csv") %>%
  dplyr::select(geneid1,symbol1,geneid2,symbol2) %>%
  unique()

# This version was supplied by the authors and includes scores
CGIdb <- read_csv("data/SL/CGIdb/CGIdb_all_prediction_table.csv") %>%
  filter(type=="SL") %>%
  dplyr::select(gene_id_1,gene_id_2,score) %>%
  unique() %>%
  inner_join(CGIdb_old, by=c("gene_id_1"="geneid1","gene_id_2"="geneid2"))


max_score <- max(CGIdb$score)
min_score <- min(CGIdb$score)

CGIdb <- CGIdb %>%
  mutate(confidence = (score-min_score)/(max_score-min_score)) %>%
  dplyr::select(gene1=symbol1,gene2=symbol2,confidence)

CGIdb <- CGIdb %>% filter(gene1 %in% genes & gene2 %in% genes)

suppressWarnings(rm(max_score,min_score))


######################
# Prepare SiLi
######################

SiLi <- read_csv("data/SL/SiLi/tableS3.csv") %>%
  # RAS = aggregated ranking scores, higher is better
  separate(col = "TSG-DT pairs", into = c("gene1","gene2"), sep = "-")

max_RAS <- max(SiLi$RAS)
min_RAS <- min(SiLi$RAS)

SiLi <- SiLi %>%
  mutate(confidence = (RAS-min_RAS)/(max_RAS-min_RAS)) %>%
  dplyr::select(gene1,gene2,confidence) %>%
  unique()

SiLi <- SiLi %>% filter(gene1 %in% genes & gene2 %in% genes)

suppressWarnings(rm(max_RAS,min_RAS))

######################
# Prepare SLOAD
######################

SLOAD <- lapply(list.files("data/SL/SLOAD/", pattern = "RESULT.TXT$"), function(x){fread(paste0("data/SL/SLOAD/",x), sep = "\t") %>% as.data.frame() %>% mutate(cancer = gsub("_RESULT.TXT","",x))})
SLOAD <- do.call(rbind,SLOAD) %>%
  group_by(GENEA,GENEB,verified) %>%
  summarise(n_cancers=n())

max_occurences <- max(SLOAD$n_cancers)
min_occurences <- min(SLOAD$n_cancers)

SLOAD <- SLOAD %>%
  mutate(occurence_score = (n_cancers-min_occurences)/(max_occurences-min_occurences)) %>%
  mutate(total_score = occurence_score+verified)

max_score <- max(SLOAD$total_score)
min_score <- min(SLOAD$total_score)

SLOAD <- SLOAD %>%
  mutate(confidence = (total_score-min_score)/(max_score-min_score)) %>%
  dplyr::select(gene1=GENEA,gene2=GENEB,confidence) %>%
  unique()

SLOAD <- SLOAD %>% filter(gene1 %in% genes & gene2 %in% genes)

rm(max_occurences,min_occurences,max_score,min_score)

######################
# Prepare Slorth
######################

slorth <- fread("data/SL/Slorth/h.sapiens_ssl_predictions.csv", sep = "\t", col.names = c("gene1","gene2","ENSG1","ENSG2","species","source","n_sources","score")) %>%
  filter(source=="Slorth") %>%
  dplyr::select(gene1,gene2,score)

max_score <- max(slorth$score)
min_score <- min(slorth$score)

slorth <- slorth %>%
  mutate(confidence = (score-min_score)/(max_score-min_score)) %>%
  dplyr::select(gene1,gene2,confidence) %>%
  unique()

slorth <- slorth %>% filter(gene1 %in% genes & gene2 %in% genes)

rm(max_score,min_score)

######################
# Prepare synlethdb
######################

synlethdb <- fread("data/SL/synlethdbv2/Human_SL.csv", col.names = c("gene1","id1","gene2","id2","cell_line","pubmed_id","source","score"))

max_score <- max(synlethdb$score)
min_score <- min(synlethdb$score)

synlethdb <- synlethdb %>%
  mutate(confidence = (score-min_score)/(max_score-min_score)) %>%
  dplyr::select(gene1,gene2,confidence) %>%
  unique()

synlethdb <- synlethdb %>% filter(gene1 %in% genes & gene2 %in% genes)

rm(max_score,min_score)



######################
# SL Pair Aggregation
######################

# Combined all sources together

SL <- rbind(
  BioGRID_SSL %>% mutate(source="BioGRID", confidence = 2*(confidence + 1)),   # Increasing scores for experimentally validated databases
  CGIdb %>% mutate(source="CGIdb"),
  gMCS_SL_pairs %>% mutate(source="gMCS"),
  SiLi %>% mutate(source="SiLi"),
  SLOAD %>% mutate(source="SLOAD"),
  slorth %>% mutate(source="slorth"),
  synlethdb %>% mutate(source="synlethdb", confidence = 2*(confidence + 1))   # Increasing scores for experimentally validated databases
)

# Copy all pairs in reverse

SL <- rbind(
  SL,
  SL %>% dplyr::select(gene1=gene2,gene2=gene1,confidence,source)
) %>%
  group_by(gene1,gene2,source) %>%
  summarise(confidence=max(confidence)) %>%
  ungroup() 

# Combine scores

SL <- SL %>%
  pivot_wider(names_from = "source", values_from = "confidence") %>%
  mutate(final_score = rowSums(!is.na(dplyr::select(., -c(gene1, gene2)))) + rowSums(dplyr::select(., -c(gene1, gene2)), na.rm=T)) %>%
  group_by(gene1) %>%
  arrange(desc(final_score)) %>%
  mutate(rank=row_number()) %>%
  filter(rank <= max_SL) %>%
  ungroup() %>%
  arrange(gene1, rank)

max_score <- max(SL$final_score)
min_score <- min(SL$final_score)

SL <- SL %>% mutate(score = (final_score-min_score)/(max_score-min_score)) 

write_csv(SL,paste0("data/SL_partners_max_",max_SL,"_all_evidence.csv"))

write_csv(SL %>% dplyr::select(gene1,gene2,rank,score),paste0("data/SL_partners_max_",max_SL,".csv"))

