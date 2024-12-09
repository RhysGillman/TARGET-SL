#benchmark_essential_gene_predictions.R


# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(ggpubr, quietly = T))
suppressPackageStartupMessages (library(coin, quietly = T))
suppressPackageStartupMessages (library(svglite, quietly = T))
suppressPackageStartupMessages (library(ggh4x, quietly = T))
suppressPackageStartupMessages (library(patchwork, quietly = T))

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



alg_colours <- read_csv("data/algorithm_colours.csv") %>% deframe()
if(any(!algorithms %in% names(alg_colours))){
  no_col_algs <- algorithms[!algorithms %in% names(alg_colours)]
  add_cols <- rep("#999999", length(no_col_algs))
  names(add_cols) <- no_col_algs
  alg_colours <- alg_colours %>% append(add_cols)
  
}

###############
# Sample Info #
###############

sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv"))
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

gold_standard <- fread("benchmark_data/all_gold_standards.csv") %>%
  dplyr::select(sample_ID=cell_ID,gene_ID) %>%
  unique() %>%
  group_by(sample_ID) %>%
  summarise(sensitive_genes = list(gene_ID)) %>%
  deframe()

rare_gold_standard <- fread("benchmark_data/rare_gold_standards.csv") %>%
  dplyr::select(sample_ID=cell_ID,gene_ID) %>%
  unique() %>%
  group_by(sample_ID) %>%
  summarise(sensitive_genes = list(gene_ID)) %>%
  deframe()

gene_effect_z_scores <- fread("benchmark_data/gene_effect_z_scores.csv")

traditional_benchmark_essential_genes <- function(predictions,gold_standard,n_max){
  
  predictions <- predictions %>% filter(final_rank <= n_max)
  
  tmp_samples <- intersect(predictions$sample_ID, names(gold_standard))
  
  result <- foreach(sample=tmp_samples, .combine = "rbind", .packages = c("tidyverse","foreach")) %do% {
    
    lineage <- predictions %>% filter(sample_ID==sample) %>% pull(cancer_type) %>% head(1)
    gs <- gold_standard[sample] %>% unlist()
    tmp1 <- predictions %>% filter(sample_ID == sample)
    
    foreach(alg=unique(tmp1$algorithm), .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
      
      tmp2 <- tmp1 %>% filter(algorithm == alg)
      if(nrow(tmp2)==0){break}
      n_predictions <- max(tmp2 %>% pull(final_rank))
      stop_at <- min(n_predictions,n_max)
      
      message(paste0("Calculating stats for ", sample, "(",which(tmp_samples == sample),"/",length(tmp_samples),")", " targets using ", alg))
      
      foreach(n=1:stop_at, .combine = "rbind", .packages = c("tidyverse")) %do% {
        
        predicted <- tmp2 %>% filter(final_rank <= n) %>% pull(target)
        correct <- predicted[which(predicted %in% gs)]
        wrong <- predicted[which(!predicted %in% gs)]
        TP <- length(correct)
        FP <- length(wrong)
        precision <- TP/n
        recall <- TP/length(gs)
        F1 <- 2*((precision*recall)/(precision + recall))
        if(is.nan(F1)){
          F1 <- 0
        }
        data.frame(cancer_type=lineage,
                   sample_ID = sample, 
                   algorithm = alg, 
                   n = n,
                   correct = paste0(correct, collapse = ";"), 
                   TP = TP,
                   FP = FP,
                   precision = precision,
                   recall = recall,
                   F1 = F1,
                   length_gs = length(gs))
        
        
      }
      
      
      
    }
    
  }
  
  return(result)
  
}

summarise_stats <- function(stats_df,alg_of_interest,min_gs,min_n) {
  n_samples <- length(unique(stats_df$sample_ID))
  stats_df_summarised <- stats_df %>%
    # Combining all badDriver simulations to one mean
    mutate(algorithm=ifelse(str_detect(algorithm,"randomDriver"), "randomDriver", algorithm)) %>%
    mutate(algorithm = ifelse(str_detect(algorithm,"consensus"), "Consensus", algorithm)) %>%
    filter(algorithm %in% alg_of_interest) %>%
    # Only keeping results for samples with > min_gs gold standards
    filter(length_gs >= min_gs) %>%
    group_by(algorithm,n) %>%
    # Only keep measurements where more than min_n cells are available to calculate mean
    filter(n()>=min_n) %>%
    summarise(
      mean_TP = mean(TP),
      mean_precision = mean(precision),
      mean_recall = mean(recall),
      mean_F1 = median(F1),
      sample_size = n()
    ) %>%
    pivot_longer(cols = -c(algorithm,n,sample_size), names_to = "measure", values_to = "value") %>%
    mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1")),
           algorithm = factor(algorithm, levels = names(alg_colours)),
           sample_size = ifelse(algorithm=="randomDriver",n_samples,sample_size))
  return(stats_df_summarised)
}

plot_with_opacity <- function(summarised_stats, plot_measure, N_max, title=NULL){
  
  plot_data <- summarised_stats %>% filter(measure==plot_measure)
  
  alg_linetype <- plot_data$algorithm %>% 
    unique() %>%
    sort() %>%
    as.character()
  alg_linetype[alg_linetype=="Consensus"] <- "dotted"
  alg_linetype[alg_linetype=="randomDriver"] <- "dashed"
  alg_linetype[!alg_linetype%in%c("dotted","dashed")] <- "solid"
  
  p1 <- ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
    geom_line(data=plot_data %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "Consensus|randomDriver")), mapping=aes(alpha = sample_size), size = 0.8) +
    geom_line(data=plot_data %>% filter(n<=N_max) %>% filter(algorithm=="Consensus"), linetype = "dotted", size = 0.8) +
    geom_line(data=plot_data %>% filter(n<=N_max) %>% filter(algorithm=="randomDriver"), linetype = "dashed", size = 0.8) +
    scale_color_manual(breaks = names(alg_colours),values = alg_colours) +
    ylab(names(plot_measure)) +
    xlab("Number of Predicted Sensitive Genes") +
    guides(colour=guide_legend(title="Algorithm (Colour)", override.aes = list(linetype = alg_linetype)),
           alpha="none",
    ) +
    theme_bw() +
    coord_cartesian(ylim = c(0,NA)) +
    ggtitle(title) +
    theme(text = element_text(size = 10), 
          legend.text = element_text(size=7), 
          legend.title = element_text(size=8,hjust = 0.5),
          legend.key.height = unit(0.3,"cm"),
          legend.key.width = unit(1.5,"cm"),
          plot.title = element_text(hjust = 0.5, size = 8)
    )


  p2 <- ggplot(plot_data %>% filter(!str_detect(algorithm, "Consensus|randomDriver")), aes(colour=sample_size, x=n,y=value)) +
    geom_line(alpha=0) +
    scale_color_gradient(high = "black",low = "white",
                         breaks = c(
                           min(plot_data %>% filter(!str_detect(algorithm, "Consensus|randomDriver")) %>% pull(sample_size)) %>% round(0),
                           max(plot_data %>% filter(!str_detect(algorithm, "Consensus|randomDriver")) %>% pull(sample_size)) %>% round(0)
                         )
    ) +
    guides(color=guide_colorbar(title = "Sample Size Remaining\n(Opacity)")) +
    theme(legend.title = element_text(size=8, hjust=0.5), legend.key.height = unit(0.3,"cm"),
          axis.title = element_blank(), axis.line = element_blank(), axis.text = element_blank(), rect = element_blank(), axis.ticks = element_blank())
  
  p1 + p2 + plot_layout(guides = "collect", widths = c(10,0)) & theme(legend.position = "right")
}

plot_without_opacity <- function(summarised_stats, plot_measure, N_max, title=NULL){
  
  plot_data <- summarised_stats %>% filter(measure==plot_measure)
  
  alg_linetype <- plot_data$algorithm %>% 
    unique() %>%
    sort() %>%
    as.character()
  alg_linetype[alg_linetype=="Consensus"] <- "dotted"
  alg_linetype[alg_linetype=="randomDriver"] <- "dashed"
  alg_linetype[!alg_linetype%in%c("dotted","dashed")] <- "solid"
  
  ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
    geom_line(data=plot_data %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "Consensus|randomDriver")), size = 0.8) +
    geom_line(data=plot_data %>% filter(n<=N_max) %>% filter(algorithm=="Consensus"), linetype = "dotted", size = 0.8) +
    geom_line(data=plot_data %>% filter(n<=N_max) %>% filter(algorithm=="randomDriver"), linetype = "dashed", size = 0.8) +
    scale_color_manual(breaks = names(alg_colours),values = alg_colours) +
    ylab(names(plot_measure)) +
    xlab("Number of Predicted Sensitive Genes") +
    guides(colour=guide_legend(title="Algorithm (Colour)", override.aes = list(linetype = alg_linetype)),
           alpha="none",
    ) +
    theme_bw() +
    coord_cartesian(ylim = c(0,NA)) +
    ggtitle(title) +
    theme(text = element_text(size = 10), 
          legend.text = element_text(size=7), 
          legend.title = element_text(size=8,hjust = 0.5),
          legend.key.height = unit(0.3,"cm"),
          legend.key.width = unit(1.5,"cm"),
          plot.title = element_text(hjust = 0.5, size = 8)
    )
}

########################
# All vs All Benchmark #
########################

if(!file.exists(paste0("results/benchmark/network_",network_choice,"/stats_essential_genes_all_vs_all.csv"))){

all_vs_all_results <- traditional_benchmark_essential_genes(all_essential_genes,gold_standard,n_predictions)

write_csv(all_vs_all_results, paste0("results/benchmark/network_",network_choice,"/stats_essential_genes_all_vs_all.csv"))
}else{
  message("Prediction stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/benchmark/network_",network_choice,"/stats_essential_genes_all_vs_all.csv"))
  all_vs_all_results <- fread(paste0("results/benchmark/network_",network_choice,"/stats_essential_genes_all_vs_all.csv"))
}

all_vs_all_results_summarised <- summarise_stats(all_vs_all_results,algorithms,10,10)

plot_without_opacity(summarised_stats=all_vs_all_results_summarised, 
                  plot_measure=c("Mean Precision"="mean_precision"), 
                  N_max=10, 
                  title="All Predictions vs All Ground Truth Essential Genes")

##########################
# Rare vs Rare Benchmark #
##########################

if(!file.exists(paste0("results/benchmark/network_",network_choice,"/stats_essential_genes_rare_vs_rare.csv"))){
  
  rare_vs_rare_results <- traditional_benchmark_essential_genes(rare_essential_genes,rare_gold_standard,n_predictions)
  
  write_csv(rare_vs_rare_results, paste0("results/benchmark/network_",network_choice,"/stats_essential_genes_rare_vs_rare.csv"))
}else{
  message("Prediction stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/benchmark/network_",network_choice,"/stats_essential_genes_rare_vs_rare.csv"))
  rare_vs_rare_results <- fread(paste0("results/benchmark/network_",network_choice,"/stats_essential_genes_rare_vs_rare.csv"))
}

rare_vs_rare_results_summarised <- summarise_stats(rare_vs_rare_results,algorithms,10,10)

plot_without_opacity(summarised_stats=rare_vs_rare_results_summarised, 
                  plot_measure=c("Mean Precision"="mean_precision"), 
                  N_max=10, 
                  title="Rare Predictions vs Rare Ground Truth Essential Genes")

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
  filter(algorithm %in% algorithms) %>%
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
        legend.title = element_text(size=8)
  )


ggsave("", device = svglite, width = 15, height = 10, units = "cm")