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
                                                                   "randomDriver",
                                                                   "Consensus"),
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

if(threads>1){
  cl <- makeCluster(threads, outfile = "log/benchmark_essential_genes.log")
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

#####################################
# Read In Predicted Essential Genes #
#####################################

message(paste0("Reading in essential gene predictions"))
all_essential_genes <- fread(paste0("results/benchmark/network_",network_choice,"/all_essential_genes.csv")) %>%
  filter(sample_ID %in% samples)
message(paste0("Reading in rare essential gene predictions"))
rare_essential_genes <- fread(paste0("results/benchmark/network_",network_choice,"/rare_essential_genes.csv"))%>%
  filter(sample_ID %in% samples)

#######################################
# Read In Predicted Ground Truth Data #
#######################################

gold_standard <- fread(paste0("benchmark_data/",network_choice,"/gold_standards.csv")) %>%
  dplyr::select(sample_ID=cell_ID,gene_ID) %>%
  unique() %>%
  group_by(sample_ID) %>%
  summarise(sensitive_genes = list(gene_ID)) %>%
  deframe()

rare_gold_standard <- fread(paste0("benchmark_data/",network_choice,"/rare_gold_standards.csv")) %>%
  dplyr::select(sample_ID=cell_ID,gene_ID) %>%
  unique() %>%
  group_by(sample_ID) %>%
  summarise(sensitive_genes = list(gene_ID)) %>%
  deframe()

gene_effect_z_scores <- fread("benchmark_data/gene_effect_z_scores.csv")


########################
# All vs All Benchmark #
########################

if(!file.exists(paste0("cache/stats_essential_genes_all_vs_all.csv"))){

all_vs_all_results <- traditional_benchmark(all_essential_genes,gold_standard,n_predictions,"gene")

write_csv(all_vs_all_results, paste0("results/benchmark/network_",network_choice,"/stats_essential_genes_all_vs_all.csv"))
write_csv(all_vs_all_results, paste0("cache/stats_essential_genes_all_vs_all.csv"))
}else{
  message("*****************Stats file already in cache. Reading in previous results.*****************")
  message(paste0("*****************To stop this, clear cache*****************"))
  all_vs_all_results <- fread(paste0("cache/stats_essential_genes_all_vs_all.csv"))
}

all_vs_all_results_summarised <- summarise_stats(all_vs_all_results,algorithms,10,10)

plot_without_opacity(summarised_stats=all_vs_all_results_summarised, 
                  plot_measure=c("Mean Precision"="mean_precision"), 
                  N_max=10, 
                  title="All Predicted Essential Genes vs All Ground Truth Essential Genes")

ggsave("plots/benchmark/essential_genes_all_vs_all.svg",device = svglite,width = 15,height = 8, units = "cm")
ggsave("plots/benchmark/essential_genes_all_vs_all.png", width = 15,height = 8, units = "cm")

top_n_plot_mean_CI(title="All Predicted Essential Genes vs All Ground Truth Essential Genes (Top 10)",
                   stats_df = all_vs_all_results,
                   algs = algorithms,
                   N_max = 10,
                   measures = "precision",
                   comparisons = list(c("Consensus","CSN_NCUA"), c("Consensus","OncoImpact")),
                   y_pos = c(0.45,0.55),
                   y_pos_randomDriver = 0.65)

ggsave("plots/benchmark/essential_genes_precision_all_vs_all_top10.svg", device = svglite, width = 15, height = 8, units = "cm")
ggsave("plots/benchmark/essential_genes_precision_all_vs_all_top10.png", width = 15, height = 8, units = "cm")

##########################
# Rare vs Rare Benchmark #
##########################

if(!file.exists(paste0("cache/stats_essential_genes_rare_vs_rare.csv"))){
  
  rare_vs_rare_results <- traditional_benchmark(rare_essential_genes,rare_gold_standard,n_predictions,"gene")
  
  write_csv(rare_vs_rare_results, paste0("cache/stats_essential_genes_rare_vs_rare.csv"))
  write_csv(rare_vs_rare_results, paste0("results/benchmark/network_",network_choice,"/stats_essential_genes_rare_vs_rare.csv"))
}else{
  message("*****************Stats file already in cache. Reading in previous results.*****************")
  message(paste0("*****************To stop this, clear cache*****************"))
  rare_vs_rare_results <- fread(paste0("cache/stats_essential_genes_rare_vs_rare.csv"))
}

rare_vs_rare_results_summarised <- summarise_stats(rare_vs_rare_results,algorithms,10,10)

plot_without_opacity(summarised_stats=rare_vs_rare_results_summarised, 
                  plot_measure=c("Mean Precision"="mean_precision"), 
                  N_max=10, 
                  title="Rare Predicted Essential Genes vs Rare Ground Truth Essential Genes")

ggsave("plots/benchmark/essential_genes_rare_vs_rare.svg",device = svglite,width = 15,height = 8, units = "cm")
ggsave("plots/benchmark/essential_genes_rare_vs_rare.png", width = 15,height = 8, units = "cm")

top_n_plot_mean_CI(title="Rare Predicted Essential Genes vs Rare Ground Truth Essential Genes (Top 10)",
                   stats_df = rare_vs_rare_results,
                   algs = algorithms,
                   N_max = 10,
                   measures = "precision",
                   comparisons = list(c("Consensus","CSN_NCUA"),c("Consensus","DawnRank"),c("Consensus","OncoImpact"),c("Consensus","sysSVM2")),
                   y_pos = c(0.07,0.08,0.09,0.10),
                   y_pos_randomDriver = 0.12)

ggsave("plots/benchmark/essential_genes_precision_rare_vs_rare_top10.svg", device = svglite, width = 15, height = 8, units = "cm")
ggsave("plots/benchmark/essential_genes_precision_rare_vs_rare_top10.svg", width = 15, height = 8, units = "cm")

#####################
# Quantitative Plot #
#####################

top_colour = 25
# For the quantitative plots, still onyl plot samples with a decent number of rare-gold standards
quant_plot_samples <- names(rare_gold_standard[lapply(rare_gold_standard, length)>=10])

gene_effect_z_scores <- fread("benchmark_data/gene_effect_z_scores.csv") %>%
  dplyr::rename(sample_ID=cell_ID, cancer_type=lineage)

rare_vs_rare_quantitative_results <- foreach(alg=unique(rare_essential_genes$algorithm), .combine = "rbind", .export = c("rare_essential_genes"), .packages = c("tidyverse","foreach")) %dopar% {
  
  message(paste0(
    "Getting cumulative means for ", 
    alg,
    " (", 
    which(unique(rare_essential_genes$algorithm)==alg), 
    "/",
    length(unique(rare_essential_genes$algorithm)),
    ")"
    
  ))
  
  # Get SL-partner level results for specific algorithm
  alg_results <- rare_essential_genes %>% 
    # Only keep cells with enough gold-standards
    filter(sample_ID %in% quant_plot_samples) %>%
    filter(algorithm==alg) %>% 
    filter(final_rank <= 100) %>%
    # Just keep the necessary info
    dplyr::select(algorithm,cancer_type,sample_ID,target,final_rank) %>%
    # Only keep data with > 10 replicates
    group_by(algorithm,final_rank) %>%
    filter(n()>=10) %>%
    ungroup()
  
  # Combine this information with the gene gene essentiality and z-scores
  
  alg_results <- alg_results %>%
    left_join(gene_effect_z_scores,
              by = c("target"="gene_ID","sample_ID","cancer_type"), 
              relationship = "many-to-one")
  
  
  foreach(rank_n=seq(1,100), .combine = "rbind") %do% {
    # Get cumulative mean values for each rank
    
    rank_n_results <- alg_results %>%
      filter(final_rank <= rank_n) %>%
      group_by(algorithm) %>%
      summarise(mean_local_z = mean(local_z_score, na.rm = T), 
                mean_global_z = mean(global_z_score, na.rm = T), 
                mean_weighted_z = mean(weighted_average, na.rm = T), 
                mean_gene_effect = mean(gene_effect, na.rm = T),
                .groups = "drop" ) %>%
      ungroup() %>%
      mutate(final_rank=rank_n)
    
    rank_n_results
    
    
  }
  
  
}


plot_data <- rare_vs_rare_quantitative_results %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"randomDriver"), "randomDriver", algorithm)) %>%
  group_by(algorithm,final_rank) %>%
  summarise(mean_local_z=mean(mean_local_z, na.rm = T), 
            mean_global_z=mean(mean_global_z, na.rm = T),
            mean_weighted_z=mean(mean_weighted_z, na.rm = T),
            mean_gene_effect=mean(mean_gene_effect, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(algorithm = ifelse(str_detect(algorithm,"consensus"), "Consensus", algorithm)) %>%
  filter(toupper(algorithm) %in% toupper(algorithms)) %>%
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


ggplot(plot_data, aes(x=mean_weighted_z, y=mean_gene_effect, colour=point_scale)) +
  geom_point(size=0.9) +
  scale_colour_gradient(high = "red", low = "blue", 
                        #breaks = seq(top_colour,0), 
                        #labels = c(seq(1,top_colour),"100 (cap)"),
                        breaks=c(25,20,15,10,5,1,0),
                        labels=c(1,5,10,15,20,25,100),
                        guide = "colourbar") +
  geom_vline(xintercept = 0, colour = "black", alpha = 0.25) +
  geom_hline(yintercept = 0, colour = "black", alpha = 0.25) +
  guides(colour=guide_colourbar(title="Top N Predictions")) +
  labs(x="Gene Uniqueness Index (Cumulative Average)", y= "Gene Effect (Cumulative Average)") +
  theme(panel.background = element_rect(fill="lightgrey")) +
  facet_wrap2(~algorithm, strip = strip) +
  theme_bw() +
  theme(legend.key.height = unit(0.7,"cm"),
        text = element_text(size = 10), 
        legend.text = element_text(size=7), 
        legend.title = element_text(size=5)
  )


ggsave("plots/benchmark/essential_genes_quantitative.svg", device = svglite, width = 15, height = 10, units = "cm")
ggsave("plots/benchmark/essential_genes_quantitative.png", width = 15, height = 10, units = "cm")
