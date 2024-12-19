#benchmark_essential_gene_predictions.R


# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(coin, quietly = T))
suppressPackageStartupMessages (library(svglite, quietly = T))
suppressPackageStartupMessages (library(patchwork, quietly = T))
suppressPackageStartupMessages (library(ggpubr, quietly = T))
suppressPackageStartupMessages (library(ggh4x, quietly = T))
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
                                                                   "Consensus",
                                                                   "pandrugs2"),
              help="algorithms to include in comparison separated by semicolons, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Threads"),
  make_option(c("-c", "--cancertype"), type="character", default="all", 
              help="Cancer types to include in the analysis separated by semicolons, or 'ALL' (Default)", metavar ="Cancer Type"),
  make_option(c("-N", "--n_predictions"), type="integer", default=10, 
              help="Number of Predictions to Benchmark", metavar ="N Predictions"),
  make_option(c("-p", "--pandrugs2"), type="character", default="../PanDrugs2/results/", 
              help="Location of PanDrugs2 predictions for all samples", metavar ="PanDrugs2 Predictions")
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
panDrugs2_path <- opt$pandrugs2

if(threads>1){
  cl <- makeCluster(threads, outfile = "log/benchmark_drug_predictions.log")
  registerDoParallel(cl)
}

if(toupper(algorithms[1])!="ALL"){
  algorithms <- append(algorithms,c("randomDriver","randomDrug"))
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

###########################
# Read In Predicted Drugs #
###########################

#randomDrug

randomDrug_predictions <- foreach(dir=list.dirs(paste0("results/benchmark/network_",network_choice,"/randomDrug"), recursive = F),.combine = "rbind") %do% {
  foreach(result_file=list.files(dir), .combine = "rbind") %do% {
            fread(paste0(dir,"/",result_file))
          }
} %>% dplyr::rename(drug_rank=rank)


message(paste0("Reading in drug predictions"))
all_drug_predictions <- fread(paste0("results/benchmark/network_",network_choice,"/all_drug_predictions.csv")) %>%
  bind_rows(randomDrug_predictions) %>%
  filter(sample_ID %in% samples)
message(paste0("Reading in rare drug predictions"))
rare_drug_predictions <- fread(paste0("results/benchmark/network_",network_choice,"/rare_drug_predictions.csv"))%>%
  bind_rows(randomDrug_predictions) %>%
  filter(sample_ID %in% samples)


##################
# PanDrugs2
##################

if("PANDRUGS2"%in%toupper(algorithms)){
  
all_drug_synonyms <- fread("data/drug_ids/all_drug_synonyms.csv")

pandrugs2_results_raw <- foreach(result_path=list.files(panDrugs2_path), .combine = "rbind") %do%{
  indiv_cell <- gsub("_pandrugs2.csv","",result_path)
  indiv_result <- fread(paste0(panDrugs2_path,"/",result_path)) %>% 
    mutate(sample_ID=indiv_cell, algorithm="pandrugs2") %>%
    #get the raw ranking for each drug
    mutate(raw_rank=row_number())
}

pandrugs2_results <- pandrugs2_results_raw %>% 
  # First matching drug names with drug names from screening data
  ## make a column with names for matching with dug screening
  dplyr::rename(possible_IDs=drug_name) %>%
  mutate(original_name=drug_show_name) %>%
  pivot_longer(c(possible_IDs,original_name),names_to = "group", values_to = "possible_IDs") %>%
  dplyr::select(-group) %>% 
  mutate(possible_IDs = toupper(gsub("[[:punct:], ]", "", possible_IDs))) %>%
  ## join it with all available drug synonyms
  left_join(all_drug_synonyms %>% mutate(synonyms=toupper(gsub("[[:punct:], ]", "", synonyms))),
            by = c("possible_IDs"="synonyms")
  ) %>%
  unique() %>%
  # Only keep drugs predicted to be sensitive, not resistant
  filter(drug_response=="SENSITIVITY") %>%
  # Only keep results with a match in drug screening data
  filter(!is.na(original)) %>%
  # Make colnames more clear
  dplyr::rename(pandrugs_drug_ID=drug_show_name, db_match_source=source,db_match_ID=original)

# no longer ranking by g_score and d_score separately, instead just taking the
# order of results from original files which seems to combine d_score and g_score

pandrugs_rankings <- pandrugs2_results  %>%
  left_join(sample_info %>% dplyr::select(cell_ID,cancer_type=lineage), by = c("sample_ID"="cell_ID")) %>%
  dplyr::select(algorithm,cancer_type,sample_ID, drug_ID=db_match_ID, raw_rank) %>%
  # Take lowest rank for multimatches
  group_by(algorithm,cancer_type,sample_ID, drug_ID) %>%
  summarise(raw_rank=min(raw_rank)) %>%
  ungroup() %>%
  #rank by dscore
  arrange(raw_rank) %>%
  group_by(sample_ID) %>%
  mutate(rank=row_number()) %>%
  ungroup() %>%
  mutate(algorithm="pandrugs2") %>%
  mutate(n_targets=NA, min_gene_rank=NA) %>%
  dplyr::select(algorithm,cancer_type,sample_ID, drug_ID,n_targets,min_gene_rank,drug_rank=rank)

pandrugs2_samples <- unique(pandrugs2_results$sample_ID)

all_drug_predictions <- all_drug_predictions %>%
  bind_rows(pandrugs_rankings) %>%
  filter(sample_ID%in%pandrugs2_samples)

rare_drug_predictions <- rare_drug_predictions %>%
  bind_rows(pandrugs_rankings) %>%
  filter(sample_ID%in%pandrugs2_samples)
}


#############################
# Read In Ground Truth Data #
#############################

gold_standard_drug_sensitivity <- read_csv(paste0("benchmark_data/network_",network_choice,"/gold_standard_drug_sensitivity.csv"), show_col_types = F)

all_gold_standards <- gold_standard_drug_sensitivity %>%
  filter(global_sensitive) %>%
  dplyr::select(sample_ID=cell_ID,drug_ID) %>%
  unique() %>%
  group_by(sample_ID) %>%
  summarise(sensitive_drugs = list(drug_ID)) %>%
  deframe()

rare_gold_standards <- gold_standard_drug_sensitivity %>%
  filter(rare_sensitive) %>%
  dplyr::select(sample_ID=cell_ID,drug_ID) %>%
  unique() %>%
  group_by(sample_ID) %>%
  summarise(sensitive_drugs = list(drug_ID)) %>%
  deframe()

########################
# All vs All Benchmark #
########################

if(!file.exists(paste0("cache/stats_drug_sensitivity_all_vs_all.csv"))){
  
  all_vs_all_results <- traditional_benchmark(all_drug_predictions,all_gold_standards,n_predictions,"drug")
  
  write_csv(all_vs_all_results, paste0("results/benchmark/network_",network_choice,"/stats_drug_sensitivity_all_vs_all.csv"))
  write_csv(all_vs_all_results, paste0("cache/stats_drug_sensitivity_all_vs_all.csv"))
}else{
  message("*****************Stats file already in cache. Reading in previous results.*****************")
  message(paste0("*****************To stop this, clear cache*****************"))
  all_vs_all_results <- fread(paste0("cache/stats_drug_sensitivity_all_vs_all.csv"))
}

all_vs_all_results_summarised <- summarise_stats(all_vs_all_results,algorithms,10,10)

plot_without_opacity(summarised_stats=all_vs_all_results_summarised, 
                     plot_measure=c("Mean Precision"="mean_precision"), 
                     N_max=10, 
                     title="All Predicted Sensitive Drugs vs All Ground Truth Sensitive Drugs")

ggsave("plots/benchmark/drug_sensitivity_all_vs_all.svg",device = svglite,width = 15,height = 8, units = "cm")
ggsave("plots/benchmark/drug_sensitivity_all_vs_all.png", width = 15,height = 8, units = "cm")

top_n_plot_mean_CI(title="All Predicted Sensitive Drugs vs All Ground Truth Sensitive Drugs (Top 10)",
                   stats_df = all_vs_all_results,
                   algs = algorithms,
                   N_max = 10,
                   measures = "precision",
                   y_pos_ref = 0.19,
                   ref="randomDrug")

ggsave("plots/benchmark/drug_sensitivity_precision_all_vs_all_top10.svg", device = svglite, width = 19, height = 12, units = "cm")
ggsave("plots/benchmark/drug_sensitivity_precision_all_vs_all_top10.png", width = 19, height = 12, units = "cm")

##########################
# Rare vs Rare Benchmark #
##########################

if(!file.exists(paste0("cache/stats_drug_sensitivity_rare_vs_rare.csv"))){
  
  rare_vs_rare_results <- traditional_benchmark(rare_drug_predictions,rare_gold_standards,n_predictions,"drug")
  
  write_csv(rare_vs_rare_results, paste0("cache/stats_drug_sensitivity_rare_vs_rare.csv"))
  write_csv(rare_vs_rare_results, paste0("results/benchmark/network_",network_choice,"/stats_drug_sensitivity_rare_vs_rare.csv"))
}else{
  message("*****************Stats file already in cache. Reading in previous results.*****************")
  message(paste0("*****************To stop this, clear cache*****************"))
  rare_vs_rare_results <- fread(paste0("cache/stats_drug_sensitivity_rare_vs_rare.csv"))
}

rare_vs_rare_results_summarised <- summarise_stats(rare_vs_rare_results,algorithms,10,10)

plot_without_opacity(summarised_stats=rare_vs_rare_results_summarised, 
                     plot_measure=c("Mean Precision"="mean_precision"), 
                     N_max=10, 
                     title="Rare Predicted Sensitive Drugs vs Rare Ground Truth Sensitive Drugs")

ggsave("plots/benchmark/drug_sensitivity_rare_vs_rare.svg",device = svglite,width = 15,height = 8, units = "cm")
ggsave("plots/benchmark/drug_sensitivity_rare_vs_rare.png", width = 15,height = 8, units = "cm")

top_n_plot_mean_CI(title="Rare Predicted Sensitive Drugs vs Rare Ground Truth Sensitive Drugs (Top 10)",
                   stats_df = rare_vs_rare_results,
                   algs = algorithms,
                   N_max = 10,
                   measures = "precision",
                   y_pos_ref = 0.12,
                   ref="randomDrug")

ggsave("plots/benchmark/drug_sensitivity_precision_rare_vs_rare_top10.svg", device = svglite, width = 15, height = 8, units = "cm")
ggsave("plots/benchmark/drug_sensitivity_precision_rare_vs_rare_top10.svg", width = 15, height = 8, units = "cm")

#####################
# Quantitative Plot #
#####################


# For the quantitative plots, still onyl plot samples with a decent number of rare-gold standards
#quant_plot_samples <- names(rare_gold_standards[lapply(rare_gold_standards, length)>=10])

quant_precily_only <- F
quant_pandrugs_only <- T

if(quant_precily_only){
  quant_plot_samples <- precily_samples
}else if(quant_pandrugs_only){
  quant_plot_samples <- intersect(pandrugs2_samples, names(rare_gold_standards[lapply(rare_gold_standards, length)>=10]))
}else{
  quant_plot_samples <- names(rare_gold_standards[lapply(rare_gold_standards, length)>=10])
}

top_colour = 25

drug_sensitivity <- fread("benchmark_data/drug_sensitivity.csv")

drug_sensitivity_results <- foreach(alg=unique(rare_drug_predictions$algorithm), .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
    
    message(paste0(
      "Getting cumulative means for ", 
      alg,
      " (", 
      which(unique(rare_drug_predictions$algorithm)==alg), 
      "/",
      length(unique(rare_drug_predictions$algorithm)),
      ")"
      
    ))
    
    # Get SL-partner level results for specific algorithm
    alg_results <- rare_drug_predictions %>% 
      # Only keep cells with enouhg gold-standards
      filter(sample_ID %in% quant_plot_samples) %>%
      filter(algorithm==alg) %>% 
      filter(drug_rank <= 100) %>%
      # Just keep the necessary info
      dplyr::select(algorithm,cancer_type,sample_ID,drug_ID,drug_rank)
    # Only keep data with > 10 replicates
    #group_by(lineage,algorithm,rank) %>%
    #filter(n()>=10) %>%
    #ungroup()
    
    # Combine this information with the z-scores
    
    alg_results <- alg_results %>%
      left_join(drug_sensitivity,
                by = c("drug_ID","sample_ID"="cell_ID","cancer_type"="lineage"), 
                relationship = "many-to-many") %>%
      # Remove cell / drug combos that aren't available
      na.omit() %>%
      group_by(algorithm,cancer_type,sample_ID) %>%
      arrange(drug_rank) %>%
      mutate(rank=row_number()) %>%
      ungroup()
    
    
    foreach(rank_n=seq(1,100), .combine = "rbind") %do% {
      # Get cumulative mean values for each rank
      
      rank_n_results <- alg_results %>%
        filter(rank <= rank_n) %>%
        group_by(algorithm,sensitivity_source) %>%
        summarise(mean_local_z = mean(local_z, na.rm = T), 
                  mean_global_z = mean(global_z, na.rm = T), 
                  mean_weighted_z = mean(weighted_average, na.rm = T), 
                  mean_sensitivity = mean(sensitivity_value, na.rm = T),
                  .groups = "drop" ) %>%
        ungroup() %>%
        mutate(final_rank=rank_n)
      
      rank_n_results
      
    }
}

plot_source <- c("GDSC1","GDSC2")

plot_data <- drug_sensitivity_results %>%
  filter(sensitivity_source %in% plot_source) %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"randomDriver"), "randomDriver", algorithm)) %>%
  mutate(algorithm=ifelse(str_detect(algorithm,"randomDrug"), "randomDrug", algorithm)) %>%
  group_by(algorithm,final_rank) %>%
  summarise(mean_local_z=mean(mean_local_z, na.rm = T), 
            mean_global_z=mean(mean_global_z, na.rm = T),
            mean_weighted_z=mean(mean_weighted_z, na.rm = T),
            mean_sensitivity=mean(mean_sensitivity, na.rm = T)
  ) %>%
  ungroup() %>%
  filter(algorithm %in% algorithms) %>%
  mutate(algorithm = ifelse(str_detect(algorithm,"consensus"), "Consensus", algorithm)) %>%
  mutate(algorithm = factor(algorithm, levels = names(alg_colours))) %>%
  mutate(point_scale=   ifelse(
    final_rank <= top_colour, top_colour + 1 - final_rank, (1 - (   (final_rank - min(final_rank)) / (max(final_rank) - min(final_rank))   ))
  )) %>%
  arrange(desc(final_rank))



strip_colours <- alg_colours[unique(sort(unique(plot_data$algorithm)))]
text_colours <- strip_colours=="#000000"
text_colours[text_colours==TRUE] <- "white"
text_colours[text_colours==FALSE] <- "black"

strip <- strip_themed(background_x = elem_list_rect(fill = strip_colours), text_x = elem_list_text(color = text_colours))

ggplot(plot_data, aes(x=mean_weighted_z, y=mean_sensitivity, colour=point_scale)) +
  geom_point(size=0.9) +
  scale_colour_gradient(high = "red", low = "blue", 
                        #breaks = seq(top_colour,0), 
                        #labels = c(seq(1,top_colour),"100 (cap)"),
                        breaks=c(25,20,15,10,5,1,0),
                        labels=c(1,5,10,15,20,25,100),
                        guide = "colourbar") +
  geom_vline(xintercept = 0, colour = "black", alpha = 0.25) +
  #geom_hline(yintercept = 0, colour = "black", alpha = 0.25) +
  guides(colour=guide_colorbar(title="Top N Predictions")) +
  labs(x="Drug Uniqueness Index (Cumulative Average)", y= "GDSC1/2 lnIC50 (Cumulative Average)") +
  theme(panel.background = element_rect(fill="lightgrey")) +
  facet_wrap2(~algorithm, strip = strip) +
  theme_bw()+
  theme(legend.key.height = unit(0.7,"cm"),
        text = element_text(size = 10), 
        legend.text = element_text(size=7), 
        legend.title = element_text(size=8),
        axis.text = element_text(size=5)
  )

ggsave(paste0("plots/benchmark/drug_sensitivity_quantitative.png"),width = 15, height = 10, units = "cm", dpi = 300)
ggsave(paste0("plots/benchmark/drug_sensitivity_quantitative.svg"), device = svglite, width = 15, height = 10, units = "cm")