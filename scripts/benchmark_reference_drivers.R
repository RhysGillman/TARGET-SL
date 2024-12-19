# evaluate_reference_drivers.r


# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(readxl, quietly = T))
suppressPackageStartupMessages (library(svglite, quietly = T))
suppressPackageStartupMessages (library(patchwork, quietly = T))
suppressPackageStartupMessages (library(ggpubr, quietly = T))
suppressPackageStartupMessages (library(rstatix, quietly = T))



# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default=c("DawnRank",
                                                                   "OncoImpact",
                                                                   "PNC",
                                                                   "PRODIGY",
                                                                   "PersonaDrive",
                                                                   "SCS",
                                                                   "sysSVM2",
                                                                   "PhenoDriverR",
                                                                   "CSN_NCUA",
                                                                   "consensus"),
              help="algorithms to include in comparison separated by semicolons, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Threads"),
  make_option(c("-c", "--cancertype"), type="character", default="all", 
              help="Cancer types to include in the analysis separated by semicolons, or 'ALL' (Default)", metavar ="Cancer Type"),
  make_option(c("-N", "--n_predictions"), type="integer", default=10, 
              help="Number of Predictions to Benchmark", metavar ="N Predictions")
);



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


network_choice <- opt$network
algorithms <- opt$algorithms
algorithms <- str_split(algorithms,";") %>% unlist()
threads <- opt$threads
cancer <- opt$cancertype
cancer <- str_split(cancer,";") %>% unlist()
n_predictions <- opt$n_predictions


if(toupper(algorithms[1])!="ALL"){
  algorithms <- append(algorithms,"randomDriver")
}


# A minimum number of gold standards allowed to each cell. If less are available, the cell is removed from calculations
min_gs <- 10


if(threads>1){
  cl <- makeCluster(threads, outfile = "log/benchmark_reference_drivers.log")
  registerDoParallel(cl)
}



alg_colours <- read_csv("data/algorithm_colours.csv", show_col_types = F) %>% deframe()
if(any(!algorithms %in% names(alg_colours))){
  no_col_algs <- algorithms[!algorithms %in% names(alg_colours)]
  add_cols <- rep("#999999", length(no_col_algs))
  names(add_cols) <- no_col_algs
  alg_colours <- alg_colours %>% append(add_cols)
  
}

source("scripts/benchmark_functions.R")


###############
# Sample Info #
###############

sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv"), show_col_types = F)
if(toupper(cancer) != "ALL"){
  sample_info <- sample_info %>% filter(cancer_type == cancer)
}
samples <- sample_info$cell_ID %>% unique() %>% sort()



#############################
# Gene List
#############################

genes <- fread(paste0("benchmark_data/network_",network_choice,"/counts.csv"), select = "gene_ID") %>% pull(gene_ID)


#############################
# Altered Genes
#############################

# Mutations

mutation <- fread(paste0("benchmark_data/network_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID")

cnv <- fread(paste0("benchmark_data/network_",network_choice,"/cnv.csv"), select = c("gene_ID",samples)) %>%
  column_to_rownames("gene_ID")

mutation <- mutation[sort(rownames(mutation)),sort(colnames(mutation))]
cnv <- cnv[sort(rownames(cnv)),sort(colnames(cnv))]


if(!all(colnames(mutation)==colnames(cnv)) & all(rownames(mutation) == rownames(cnv))){
  stop("Error: colnames/rownames of cnv and mutation files do not match")
}


altered_genes <- cnv != 0 | mutation != 0

altered_genes <- altered_genes %>%
  as.data.frame() %>%
  rownames_to_column("gene_ID") %>%
  pivot_longer(cols = -gene_ID, names_to = "sample_ID", values_to = "altered") %>%
  filter(altered) %>%
  dplyr::select(sample_ID, gene_ID) %>%
  arrange(sample_ID)


#############################
# Read In Results
#############################

source("scripts/read_driver_results.R")

aggregated_results <- read_driver_results("benchmark",algorithms,cancer)

if(toupper(cancer) != "ALL"){
  aggregated_results <- aggregated_results %>% filter(cancer_type %in% cancer)
}

#############################
# Keywords
#############################

if(!file.exists("data/cancer_reference_genes/keywords.xlsx")){
  CGC_keywords <- read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv", show_col_types = F) %>%
    mutate(keywords = paste(`Tumour Types(Somatic)`,`Tumour Types(Germline)`,`Cancer Syndrome`, sep = " "),
           source = "CGC") %>%
    dplyr::select(source, keywords) %>%
    separate_longer_delim(cols = keywords, delim = " ") %>%
    mutate(keywords = gsub(",|NA|\\.", "", keywords)) %>%
    unique()
  
  
  NCG_keywords <- read_tsv("data/cancer_reference_genes/NCG/NCG_cancerdrivers_annotation_supporting_evidence.tsv", show_col_types = F) %>%
    mutate(keywords = paste(cancer_type,primary_site, sep = " "),
           source = "NCG") %>%
    dplyr::select(source, keywords) %>%
    separate_longer_delim(cols = keywords, delim = " ") %>%
    mutate(keywords = gsub(",|NA|\\.", "", keywords)) %>%
    unique()
  
  
  CancerMine_keywords <- read_tsv("data/cancer_reference_genes/CancerMine/cancermine_collated.tsv", show_col_types = F) %>%
   #Filtering for >1 citation count as done in PersonaDrive
    filter(citation_count >1) %>%
    mutate(keywords = paste(cancer_normalized, sep = " "),
           source = "CancerMine") %>%
    dplyr::select(source, keywords) %>%
    separate_longer_delim(cols = keywords, delim = " ") %>%
    mutate(keywords = gsub(",|NA|\\.", "", keywords)) %>%
    unique()
  
  all_keywords <- rbind(CGC_keywords,NCG_keywords,CancerMine_keywords) %>%
    dplyr::select(keywords) %>%
    unique()
  
  
  write_csv(all_keywords,"data/cancer_reference_genes//raw_keywords.csv")
}else{
  all_keywords <- read_xlsx("data/cancer_reference_genes/keywords.xlsx") 
}

all_keywords <- all_keywords %>%
  pivot_longer(cols = everything(), names_to = "cancer_type", values_to = "keyword", values_drop_na = T)



#############################
# Reference Drivers
#############################

CGC_all <- read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv", show_col_types = F) %>%
  dplyr::select(gene_ID = `Gene Symbol`) %>%
  unique()

CGC_specific <- read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv", show_col_types = F) %>%
  pivot_longer(cols = c(`Tumour Types(Somatic)`,`Tumour Types(Germline)`,`Cancer Syndrome`), names_to = "keyword_type", values_to = "keyword") %>%
  separate_longer_delim(cols = keyword, delim = ",") %>%
  mutate(keyword = gsub("^ ", "", keyword)) %>%
  dplyr::select(gene_ID = `Gene Symbol`, keyword) %>%
  rbind(
    read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv", show_col_types = F) %>%
      mutate(keyword = paste(`Tumour Types(Somatic)`,`Tumour Types(Germline)`,`Cancer Syndrome`, sep = " ")) %>%
      separate_longer_delim(cols = keyword, delim = " ") %>%
      mutate(keyword = gsub(",|NA|\\.", "", keyword)) %>%
      dplyr::select(gene_ID = `Gene Symbol`, keyword)
  ) %>%
  left_join(all_keywords, by = "keyword", multiple = "all", relationship = "many-to-many") %>%
  dplyr::select(gene_ID, cancer_type) %>%
  na.omit() %>%
  unique()


NCG_all <- read_tsv("data/cancer_reference_genes/NCG/NCG_cancerdrivers_annotation_supporting_evidence.tsv", show_col_types = F) %>%
  dplyr::select(gene_ID = symbol) %>%
  unique()

NCG_specific <- read_tsv("data/cancer_reference_genes/NCG/NCG_cancerdrivers_annotation_supporting_evidence.tsv", show_col_types = F) %>%
  pivot_longer(cols = c(cancer_type,primary_site), names_to = "keyword_type", values_to = "keyword") %>%
  separate_longer_delim(cols = keyword, delim = ",") %>%
  mutate(keyword = gsub("^ ", "", keyword)) %>%
  dplyr::select(gene_ID = symbol, keyword) %>%
  rbind(
    read_tsv("data/cancer_reference_genes/NCG/NCG_cancerdrivers_annotation_supporting_evidence.tsv", show_col_types = F) %>%
      mutate(keyword = paste(cancer_type,primary_site, sep = " ")) %>%
      separate_longer_delim(cols = keyword, delim = " ") %>%
      mutate(keyword = gsub(",|NA|\\.", "", keyword)) %>%
      dplyr::select(gene_ID = symbol, keyword)
    ) %>%
  left_join(all_keywords, by = "keyword", multiple = "all",relationship = "many-to-many") %>%
  dplyr::select(gene_ID, cancer_type) %>%
  na.omit() %>%
  unique()


CancerMine_all <- read_tsv("data/cancer_reference_genes/CancerMine/cancermine_collated.tsv", show_col_types = F) %>%
  # Filtering for >1 citation count as done in PersonaDrive
  filter(citation_count >1) %>%
  dplyr::select(gene_ID = gene_normalized) %>%
  unique()

CancerMine_specific <- read_tsv("data/cancer_reference_genes/CancerMine/cancermine_collated.tsv", show_col_types = F) %>%
  # Filtering for >1 citation count as done in PersonaDrive
  filter(citation_count >1) %>%
  dplyr::select(gene_ID = gene_normalized, keyword = cancer_normalized) %>%
  rbind(
    read_tsv("data/cancer_reference_genes/CancerMine/cancermine_collated.tsv", show_col_types = F) %>%
      separate_longer_delim(cancer_normalized, delim = " ") %>%
      mutate(keyword = gsub(",|NA|\\.", "", cancer_normalized)) %>%
      dplyr::select(gene_ID = gene_normalized, keyword)
  ) %>%
  left_join(all_keywords, by = "keyword", multiple = "all",relationship = "many-to-many") %>%
  dplyr::select(gene_ID, cancer_type) %>%
  na.omit() %>%
  unique()

reference_drivers_all <- rbind(
  CGC_all %>% mutate(source="CGC"),
  NCG_all %>% mutate(source = "NCG"),
  CancerMine_all %>% mutate(source = "CancerMine")
) %>%
  filter(gene_ID %in% genes)

reference_drivers_specific <- rbind(
  CGC_specific %>% mutate(source="CGC"),
  NCG_specific %>% mutate(source = "NCG"),
  CancerMine_specific %>% mutate(source = "CancerMine")
) %>%
  filter(gene_ID %in% genes)

rm(CGC_all,CGC_specific,NCG_all,NCG_specific,CancerMine_all,CancerMine_specific,all_keywords)


#############################
# Results
#############################




# average n results per algorithm

n_per_algorithm <- aggregated_results %>%
  #filter(algorithm %in% c(
  #  "DawnRank", "OncoImpact", "SCS", "PersonaDrive", "PhenoDriverR", "PNC", "PRODIGY", "sysSVM2", 
  #  "CSN_DFVS", "CSN_MDS", "CSN_MMS", "CSN_NCUA", "LIONESS_DFVS", "LIONESS_MDS", "LIONESS_MMS", "LIONESS_NCUA", "SPCC_DFVS", "SPCC_MDS", "SPCC_MMS", "SPCC_NCUA", "SSN_DFVS", "SSN_MDS","SSN_MMS", "SSN_NCUA"
  #)) %>%
  group_by(algorithm, sample_ID) %>%
  summarise(n_drivers=n()) %>%
  group_by(algorithm) %>%
  summarise(mean_drivers=mean(n_drivers))

write_csv(n_per_algorithm, paste0("results/benchmark/network_",network_choice,"/average_n_drivers_per_sample.csv"))
  




run_evaluation <- function(ref_set,gold_standard_type){
  
  if(gold_standard_type == "all"){
    gs <- reference_drivers_all %>% filter(source==ref_set) %>% pull(gene_ID) %>% unique()
    samples <- aggregated_results$sample_ID %>% unique()
  } else if(gold_standard_type == "specific"){
    
    if(ref_set == "all"){
      gold_standard <- reference_drivers_specific
    }else{
      gold_standard <- reference_drivers_specific %>% filter(source==ref_set)
    }
    
    samples <- aggregated_results %>% filter(cancer_type %in% gold_standard$cancer_type) %>% pull(sample_ID) %>% unique()
  }
  
    prediction_stats <- foreach(sample=samples, .combine = "rbind", .packages = c("tidyverse","foreach"), .export = c("aggregated_results","altered_genes","min_gs")) %dopar% {
      
      lin <- aggregated_results %>% filter(sample_ID==sample) %>% pull("cancer_type") %>% unique()
      
      alt_genes <- altered_genes %>% filter(sample_ID==sample) %>% pull(gene_ID)
      
      if(gold_standard_type == "specific"){
        gs_specific <- gold_standard %>% filter(cancer_type==lin) %>% pull(gene_ID) %>% unique()
        gs_specific <- gs_specific[which(gs_specific %in% alt_genes)]
      }else{
        gs_specific <- gs[which(gs %in% alt_genes)]
      }
      
      if(length(gs_specific)<min_gs){return(NULL)}else{
      
      
      
      tmp1 <- aggregated_results %>% filter(sample_ID == sample)
      
      foreach(alg=unique(tmp1$algorithm), .combine = "rbind", .packages = c("tidyverse","foreach"), .export = "run_evaluation") %do% {
        
        tmp2 <- tmp1 %>% filter(algorithm == alg)
        if(nrow(tmp2)==0){return(NULL)}else{
        max_drivers <- max(tmp2 %>% pull(rank))
        
        print(paste0("Calculating stats for ", sample, "(",which(samples == sample),"/",length(samples),")", " targets using ", alg))
        
        foreach(n=1:max_drivers, .combine = "rbind", .packages = c("tidyverse"), .export = "run_evaluation") %do% {
          
          predicted <- tmp2 %>% filter(rank <= n) %>% pull(driver)
          correct <- predicted[which(predicted %in% gs_specific)]
          wrong <- predicted[which(!predicted %in% gs_specific)]
          TP <- length(correct)
          FP <- length(wrong)
          precision <- TP/n
          recall <- TP/length(gs_specific)
          F1 <- 2*((precision*recall)/(precision + recall))
          if(is.nan(F1)){
            F1 <- 0
          }
          n_gs <- length(gs_specific)
          
          data.frame(cancer_type=lin,
                     sample_ID = sample, 
                     algorithm = alg, 
                     n = n,
                     correct = paste0(correct, collapse = ";"), 
                     TP = TP,
                     FP = FP,
                     precision = precision,
                     recall = recall,
                     F1 = F1,
                     n_gs = n_gs)
          
        }
        }
        
        
        
      }
      }
      
    }
    return(prediction_stats)
}


if(!file.exists(paste0("cache/stats_reference_drivers.csv"))){
prediction_stats <- run_evaluation("CGC","all")
write_csv(prediction_stats, paste0("results/benchmark/network_",network_choice,"/stats_reference_drivers.csv"))
write_csv(prediction_stats, paste0("cache/stats_reference_drivers.csv"))
}else{
  message("*****************Stats file already in cache. Reading in previous results.*****************")
  message(paste0("*****************To stop this, clear cache*****************"))
  prediction_stats <- read_csv(paste0("cache/stats_reference_drivers.csv"), show_col_types = F)
}



stats_summarised <- summarise_stats(prediction_stats,algorithms,10,10)


plot_without_opacity(summarised_stats=stats_summarised, 
                     plot_measure=c("Mean Precision"="mean_precision"), 
                     N_max=10, 
                     title="Driver Gene Predictions vs CGC")

ggsave("plots/benchmark/reference_drivers_precision.svg",device = svglite,width = 15,height = 8, units = "cm")
ggsave("plots/benchmark/reference_drivers_precision.png",width = 15,height = 8, units = "cm")

top_n_plot_mean_CI(title="All Predictions vs CGC (Top 10)",
                   stats_df = prediction_stats,
                   algs = algorithms,
                   N_max = 10,
                   measures = "precision",
                   comparisons = list(c("Consensus","PersonaDrive"), c("Consensus","SCS")),
                   y_pos = c(0.6,0.7),
                   y_pos_ref = 0.8,
                   ref="randomDriver")

ggsave("plots/benchmark/reference_drivers_precision_top10.svg", device = svglite, width = 19, height = 12, units = "cm")
ggsave("plots/benchmark/reference_drivers_precision_top10.png", width = 19, height = 12, units = "cm")

####################
# Unused Functions #
####################

calculate_AP <- function(precision_scores,topn){
  
  data <- precision_scores %>% dplyr::select(sample_ID, algorithm, n, TP, precision, n_gs) %>%
    filter(n<=topn) %>%
    group_by(sample_ID, algorithm) %>%
    arrange(n) %>%
    #checks whether there has been an increase in the number of TP since the last row
    mutate(relevant_retrieval=ifelse(TP > dplyr::lag(TP, n=1), 1, 0)) %>%
    #add relevance if the first result is a TP
    mutate(relevant_retrieval=ifelse(n==1,ifelse(TP==1,1,0),relevant_retrieval)) %>%
    ungroup() %>%
    #keep only relevant entries
    filter(relevant_retrieval==1) %>%
    group_by(sample_ID,algorithm) %>%
    # calculate average precision (AP)
    summarise(AP=mean(precision), .groups = "drop")
  
  
}

compare_mAP <- function(AP_scores,algs,comparisons,ypos_specific=0.8,y_pos_badDriver=0.8,title){
  
  alg_colours <- read_csv("data/algorithm_colours.csv", show_col_types = F) %>% deframe()
  
  AP_scores_trim <- AP_scores %>%
    mutate(algorithm=ifelse(str_detect(algorithm,"randomDriver"), "randomDriver", algorithm)) %>%
    filter(algorithm %in% algs) %>%
    mutate(algorithm = ifelse(str_detect(algorithm,"consensus"), "Consensus", algorithm))
  
  #ggplot(AP_scores_trim, aes(x=algorithm, y=AP)) + geom_boxplot()

  pairwise_stats <- compare_means(AP~algorithm,data=AP_scores_trim,p.adjust.method = "BH",method = "wilcox.test") %>%
    mutate(custom=ifelse(p.signif=="****", "<0.0001", "")) %>%
    mutate(custom=ifelse(p.signif=="***", "<0.001", custom)) %>%
    mutate(custom=ifelse(p.signif=="**", "<0.01", custom)) %>%
    mutate(custom=ifelse(p.signif=="*", "<0.05", custom)) %>%
    mutate(custom=ifelse(p.signif=="ns", "= ns",custom))
  
  vs_badDriver <- pairwise_stats %>% filter(group1=="randomDriver"|group2=="randomDriver") %>% mutate(xpos=ifelse(group1=="randomDriver",group2,group1))
  specific_comparisons <- pairwise_stats %>% filter(group1 %in% comparisons[[1]] & group2 %in% comparisons[[1]])
  
  plot_data <- AP_scores_trim %>%
    group_by(algorithm) %>%
    summarise(
      n=n(),
      mAP=mean(AP, na.rm = T), 
      s_AP=sd(AP, na.rm=T)
    ) %>%
    ungroup() %>%
    mutate(
      ci95_AP=qt(0.975,df=n-2)*s_AP/sqrt(n),
    ) %>%
    mutate(algorithm = factor(algorithm, levels = names(alg_colours)))
  
  ggplot(plot_data, aes(x = algorithm, y = mAP, colour = algorithm)) +
    #geom_jitter(data=AP_scores_trim, mapping = aes(y=AP,x=algorithm,colour=algorithm), alpha = 0.2) +
    geom_point(position = position_dodge(width=0.75)) +
    geom_errorbar(aes(ymin=mAP-ci95_AP,ymax=mAP+ci95_AP),
                  position = position_dodge(width=0.75)) +
    scale_colour_manual(breaks = names(alg_colours),values = alg_colours) +
    stat_pvalue_manual(specific_comparisons,label="p.signif", 
                       y.position = ypos_specific, color = "black", bracket.size = 0.5, tip.length = 0.01
    ) +
    stat_pvalue_manual(vs_badDriver, x="xpos", y.position = y_pos_badDriver, label = "p.signif") +
    ylab("Mean Average Precision +/- 95%CI") +
    xlab("") +
    theme_bw() +
    guides(colour="none") +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) +
    ggtitle(title) +
    coord_cartesian(ylim = c(0,NA))

}
