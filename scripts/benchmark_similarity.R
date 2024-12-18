#benchmark_similarity.R

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(MASS, quietly = T))
suppressPackageStartupMessages (library(ggrepel, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(UpSetR, quietly = T))
suppressPackageStartupMessages (library(svglite, quietly = T))

suppressPackageStartupMessages (library(ComplexUpset, quietly = T))
suppressPackageStartupMessages (library(ggplot2, quietly = T))
suppressPackageStartupMessages (library(patchwork, quietly = T))

# Handling input arguments
option_list = list(
  make_option(c("-m", "--mode"), type="character", default="benchmark", 
              help="mode (benchmark or predict)", metavar ="Mode"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default=c("DawnRank",
                                                                   "OncoImpact",
                                                                   "PNC",
                                                                   "PRODIGY",
                                                                   "PersonaDrive",
                                                                   "SCS",
                                                                   "sysSVM2",
                                                                   "PhenoDriverR",
                                                                   "CSN_NCUA"
                                                                   ),
              help="algorithms to include in comparison separated by semicolons (;), or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Rank Aggregation Method"),
  make_option(c("-s", "--sampleinfo"), type="character", default=NULL, 
              help="path to sample info file", metavar ="Sample Info")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

run_mode <- opt$mode
network_choice <- opt$network
algorithms <- opt$algorithms
threads <- opt$threads

algorithms <- str_split(algorithms, ";") %>% unlist()

if(run_mode=="predict"){
  sample_info_path <- opt$sampleinfo
}



if(threads>1){

  cl <- makeCluster(threads, outfile = "log/evaluate_similarity.log")
  registerDoParallel(cl)
}

run_divergence_plot <- F

#############################
# Sample Info
#############################

if(run_mode=="predict"){
  sample_info <- fread(sample_info_path)
  samples <- sample_info %>% pull(sample_ID) %>% unique() %>% sort()
}else if(run_mode=="benchmark"){
  sample_info <- read_csv(paste0("benchmark_data/network_",network_choice,"/sample_info.csv"))
  samples <- sample_info$cell_ID %>% unique() %>% sort()
}

#############################
# Read In Results
#############################

source("scripts/read_driver_results.R")

aggregated_results <- read_driver_results(run_mode,algorithms,"all")

##################
# Upset Plots    #
##################


create_custom_upset_plot <- function(n_intersects=10,upset_topn=10,upset_algorithms){
  suppressWarnings(rm(upset_alg_combinations))
  for(i in seq(2,length(upset_algorithms))){
    indiv_comb <- as.list(as.data.frame(combn(upset_algorithms,i)))
    if(!exists("upset_alg_combinations")){
      upset_alg_combinations <- indiv_comb
    }else{
      upset_alg_combinations <- append(upset_alg_combinations,indiv_comb)
    }
  }
  
  unique_sizes <- aggregated_results %>%
    dplyr::select(sample_ID,driver,algorithm,rank) %>%
    # Take top n predicted drivers for the selection of algorithms
    filter(rank <=upset_topn, algorithm %in% upset_algorithms) %>%
    # Ensure there are no duplicated drivers within one algorithm
    dplyr::select(-rank) %>%
    unique() %>%
    group_by(sample_ID,driver) %>%
    # only keep unduplicated lines
    summarise(repeats=n(), algorithm = paste(algorithm, collapse = "&")) %>%
    ungroup() %>%
    filter(repeats==1) %>%
    group_by(algorithm) %>%
    summarise(size=n())
  
  suppressWarnings(rm(intersect_sizes))
  intersect_sizes <- foreach(comb=upset_alg_combinations,.combine = "rbind") %do% {
    
    data <- aggregated_results %>%
      dplyr::select(sample_ID,driver,algorithm,rank) %>%
      # Take top n predicted drivers for the selection of algorithms
      filter(rank <=upset_topn, algorithm %in% unlist(comb)) %>%
      # Ensure there are no duplicated drivers within one algorithm
      dplyr::select(-rank) %>%
      unique() %>%
      # Get only the drivers that are repeated in every algorithm
      group_by(sample_ID,driver) %>%
      summarise(repeats=n(), .groups = "drop") %>%
      #ungroup() %>%
      filter(repeats==length(unlist(comb))) %>%
      # Sum up the number of repeated drivers
      group_by(sample_ID) %>%
      summarise(intersect_size_per_cell=n(), .groups = "drop") %>%
      #ungroup() %>%
      summarise(size=sum(intersect_size_per_cell), .groups = "drop") %>%
      mutate(set=paste(unlist(comb), collapse = "&"))
    
  }
  
  
  top_n_intersect <- intersect_sizes %>% 
    arrange(desc(size)) %>%
    head(n_intersects) %>%
    pull(set)
  
  top_n_intersect <- foreach(int=top_n_intersect) %do% {
    str_split(int,"&") %>% unlist
  }
  
  
  
  plot_intersects <- unique_sizes %>%
    arrange(desc(size)) %>%
    pull(algorithm) %>%
    append(top_n_intersect)
  
  plot_intersects <- foreach(int=plot_intersects) %do% {
    str_split(int,"&") %>% unlist
  }
  
  
  
  set_order <- aggregated_results %>%
    filter(rank <= upset_topn, algorithm %in% upset_algorithms) %>%
    dplyr::select(sample_ID,algorithm,driver) %>%
    unique() %>%
    group_by(algorithm) %>%
    summarise(size=n(), .groups = "drop") %>%
    arrange(desc(size)) %>%
    pull(algorithm)
  
  set_colours <- read_csv("data/algorithm_colours.csv") %>%
    deframe()
  set_colours <- set_colours[set_order]
  set_colours[is.na(set_colours)] <- "grey"
  names(set_colours) <- set_order
  
  
  upset_data <- aggregated_results %>%
    filter(rank<=upset_topn, algorithm %in% upset_algorithms) %>%
    mutate(find=paste0(sample_ID,"_",driver)) %>%
    dplyr::select(find,algorithm) %>%
    table() %>%
    as.data.frame() %>%
    mutate(Freq=as.logical(Freq)) %>%
    pivot_wider(names_from = algorithm, values_from = Freq)
  
  
  
  intersect_mode <- "exclusive_intersection"
  
  
  intersect_ylim = round(max(append(unique_sizes$size,intersect_sizes$size))*1.1,-1)
  
  
  p1 <- ComplexUpset::upset(data = upset_data,
                            min_size = 6,
                            intersect = set_order, 
                            sort_sets = F,
                            intersections = (unique_sizes %>% arrange(desc(size)) %>% pull(algorithm)),
                            sort_intersections= F,
                            mode= intersect_mode,
                            queries = foreach(col=names(set_colours)) %do% {
                              upset_query(set=col, fill=set_colours[col])
                            },
                            #n_intersections = 10,
                            # Customising set sizes plots
                            set_sizes = (
                              upset_set_size() +
                                geom_text(aes(label=after_stat(count)), hjust=0, stat='count') +
                                #expand_limits(y=1500*upset_topn)
                                expand_limits(y=intersect_ylim)
                            ),
                            # Customising Intersect sizes plot
                            
                            base_annotations=list(
                              'Intersection size'=intersection_size(text_mapping=aes(label=paste0(round(
                                !!get_size_mode(intersect_mode)/!!get_size_mode('inclusive_union') * 100
                              ), '%')),
                              mode=intersect_mode) + ylim(0,intersect_ylim) + ggtitle("Exclusive Predictions") + theme(plot.title = element_text(hjust = 0.5))
                            ),
                            # customise matrix
                            matrix=(
                              intersection_matrix(geom=geom_point(shape='circle filled', size=4))
                            ),
                            themes = upset_modify_themes(
                              list(
                                'intersections_matrix'=theme(axis.title.x = element_blank()#, 
                                                             #axis.text.y = element_text(colour=set_colours)
                                )
                              ))
  ) + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
  
  intersect_mode <- "inclusive_intersection"
  
  p2 <- ComplexUpset::upset(data = upset_data, 
                            intersect = set_order, 
                            sort_sets = F,
                            intersections = top_n_intersect,
                            sort_intersections= F,
                            mode= intersect_mode,
                            queries = foreach(col=names(set_colours)) %do% {
                              upset_query(set=col, fill=set_colours[col])
                            },
                            #set_sizes = F,
                            set_sizes = (
                              upset_set_size(position="right") +
                                geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count') +
                                expand_limits(y=1500*upset_topn)
                            ),
                            # Customising Intersect sizes plot
                            base_annotations=list(
                              'Intersect Size'=intersection_size(text_mapping=aes(label=paste0(round(
                                !!get_size_mode(intersect_mode)/!!get_size_mode('inclusive_union') * 100
                              ), '%')),
                              mode=intersect_mode,) + ylab("") + theme(axis.text = element_blank()) + ylim(0,intersect_ylim) + ggtitle("Inclusive Intersections") + theme(plot.title = element_text(hjust = 0.5))
                            ),
                            matrix=(
                              intersection_matrix(geom=geom_point(shape='circle filled', size=4))
                            ),
                            themes = upset_modify_themes(
                              list(
                                'intersections_matrix'=theme(axis.text.y = element_blank(), axis.title.x = element_blank())
                              )
                            )#,
                            #wrap = T
  )+ theme(plot.margin = margin(0, 0, 0, 0, "pt"))
  
  (p1 + plot_layout(widths = c(2.5,4.5)) | p2 + plot_layout(widths = c(5,2.5))) + plot_layout(widths = c(8,7))
}

create_custom_upset_plot(n_intersects = 10,upset_topn = 1,upset_algorithms = algorithms)
ggsave(paste0("plots/",run_mode,"/upset_plot_top_1.png"), width = 30, height = 15, units = "cm")
ggsave(paste0("plots/",run_mode,"/upset_plot_top_1.svg"),device = svglite, width = 30, height = 15, units = "cm")
create_custom_upset_plot(n_intersects = 10,upset_topn = 10,upset_algorithms = algorithms)
ggsave(paste0("plots/",run_mode,"/upset_plot_top_10.png"), width = 30, height = 15, units = "cm")
ggsave(paste0("plots/",run_mode,"/upset_plot_top_10.svg"),device = svglite, width = 30, height = 15, units = "cm")


if(run_divergence_plot){
  


#############################
# Divergence Plots          #
#############################

# Limit all results to those available for every algorithm
# Skip this step for upset plots

samples <- aggregated_results %>%
  dplyr::select(sample_ID,algorithm) %>%
  unique() %>%
  group_by(algorithm) %>%
  summarise(samples=list(sample_ID)) %>%
  deframe()

samples <- Reduce(intersect, samples)


aggregated_results <- aggregated_results %>% filter(sample_ID %in% samples)

top_ns <- c(1,2,3,4,5,6,7,8,9,10)
alg_pairs <- as.list(as.data.frame(combn(unique(aggregated_results$algorithm), 2)))

suppressWarnings(rm(combined_cos_sim))


if(!file.exists("cache/full_results_similarity.csv")){
  
  
  combined_sim <- foreach(alg_pair=alg_pairs, .combine = "rbind",  .packages = c("tidyverse","foreach"), .export = c("aggregated_results", "top_ns")) %dopar% {
    
    alg1 <- alg_pair[1]
    alg2 <- alg_pair[2]
    
    tmp_results_1 <- aggregated_results %>% filter(algorithm %in% c(alg1, alg2))
    
    comp_cells <- intersect(tmp_results_1 %>% filter(algorithm==alg1) %>% pull(sample_ID),
                            tmp_results_1 %>% filter(algorithm==alg2) %>% pull(sample_ID)
    )
    
    message(paste("Calculating cosine and jaccard similarity for", alg1, "vs", alg2, sep = " "))
    
    
    foreach(cell=comp_cells, .combine = "rbind", .packages = c("tidyverse","foreach"), .export = c("tmp_results_1","cell","alg1","alg2")) %do% {
      
      tmp_results_2 <- tmp_results_1 %>% filter(sample_ID == cell)
      
      foreach(n=top_ns, .combine = "rbind", .packages = c("tidyverse","foreach"), .export = c("tmp_results_2","cell","alg1","alg2")) %do% {
        
        # Calculate the cosine similarity of results from each cell
        
        cos_sim_matrix <- tmp_results_2 %>% filter(rank <= n) %>% 
          ungroup() %>% 
          dplyr::select(driver, algorithm) %>%
          table() %>%
          as.data.frame() %>%
          pivot_wider(names_from = algorithm, values_from = Freq)
        
        cos_sim <- (
          # This calculates the sum of the products of the two vectors
          as.numeric(t(pull(cos_sim_matrix,2))%*%pull(cos_sim_matrix,3))
          # Divide
        ) /
          # This multiplies the square root of the sum squares of each vector
          (sqrt(sum(pull(cos_sim_matrix,2)^2))*sqrt(sum(pull(cos_sim_matrix,3)^2)))
        
        
        # Calculate the jaccard similarity of results from each cell
        
        jac_sim <- (
          as.numeric(length(intersect(tmp_results_2 %>% filter(rank <= n) %>% filter(algorithm==alg1) %>% pull(driver),
                                      tmp_results_2 %>% filter(rank <= n) %>% filter(algorithm==alg2) %>% pull(driver)))) / 
            as.numeric(length(union(tmp_results_2 %>% filter(rank <= n) %>% filter(algorithm==alg1) %>% pull(driver),
                                    tmp_results_2 %>% filter(rank <= n) %>% filter(algorithm==alg2) %>% pull(driver))))
        )
        
        
        
        data.frame(sample_ID = cell,
                   algorithm_1 = alg1, 
                   algorithm_2 = alg2,
                   top_n = n,
                   cos_sim = cos_sim,
                   jac_sim = jac_sim
        )
      }
    }
    
    
  }
  
  write_csv(combined_sim, "cache/full_results_similarity.csv")
}else{
  message("*****************Similarity stats file already in cache. Reading in previous results.*****************")
  message(paste0("*****************To stop this, clear cache*****************"))
  combined_sim <- read_csv("cache/full_results_similarity.csv")
}


#############################
# Cosine Similarity
#############################

summarised_cos_sim <- combined_sim %>%
  group_by(algorithm_1,algorithm_2,top_n) %>%
  summarise(mean_cos_sim = mean(cos_sim, na.rm = T)) %>%
  mutate(mean_cos_dist = 1-mean_cos_sim)

if(toupper(algorithms[1]) != "ALL"){
    summarised_cos_sim <- summarised_cos_sim %>% filter(algorithm_1 %in% algorithms & algorithm_2 %in% algorithms)
}


suppressWarnings(rm(cos_similarity_plot))
for(n in top_ns){
  dist_matrix <- summarised_cos_sim %>%
    filter(top_n==n) %>%
    dplyr::select(algorithm_1,algorithm_2,mean_cos_dist)
  dist_matrix <- dist_matrix %>%
    # Need to duplicate values to make full symmetric matrix
    rbind(dist_matrix %>% dplyr::select(algorithm_1 = algorithm_2, algorithm_2 = algorithm_1, mean_cos_dist)) %>%
    unique() %>%
    # Make distance matrix
    pivot_wider(names_from = algorithm_2, values_from = mean_cos_dist, values_fill = 0) %>%
    column_to_rownames("algorithm_1") %>%
    as.matrix()
  dist_matrix <- dist_matrix[sort(unique(append(summarised_cos_sim$algorithm_1,summarised_cos_sim$algorithm_2))),sort(unique(append(summarised_cos_sim$algorithm_1,summarised_cos_sim$algorithm_2)))]
  
  #mds <- isoMDS(dist_matrix, k = 1, 
  #              tol = 0.5, 
  #              p = 2,
  #              maxit = 100)
  mds <- cmdscale(dist_matrix, k=1)
  
  tmp_similarity_plot <- mds %>%
    as.data.frame() %>%
    dplyr::rename(pos=V1) %>%
    rownames_to_column("algorithm") %>%
    mutate(n_drivers = n)
  
  if(!exists("cos_similarity_plot")){
    cos_similarity_plot <- tmp_similarity_plot
  }else{
    cos_similarity_plot <- rbind(cos_similarity_plot, tmp_similarity_plot)
    
  }
  
}

# Invert values if order of algorithms becomes reversed
for(n in top_ns){
  
  if(n!=top_ns[1]){
    current_order <- cos_similarity_plot %>% 
      filter(n_drivers==n) %>%
      arrange(desc(pos)) %>%
      pull(algorithm)
    
    current_match <- sum(current_order==prev_order)
    reverse_match <- sum(current_order==rev(prev_order))
    
    #print(paste0(current_match," vs ",reverse_match," reverse = ", reverse_match>current_match))
    
    if(reverse_match>current_match){
      cos_similarity_plot <- cos_similarity_plot %>%
        mutate(pos=ifelse(n_drivers==n,pos/-1,pos))
    }
    
    
  }
  
  prev_order <- cos_similarity_plot %>% 
    filter(n_drivers==n) %>%
    arrange(desc(pos)) %>%
    pull(algorithm)
  
}



cos_similarity_plot <- cos_similarity_plot %>%
  mutate(algorithm=factor(algorithm, levels = names(alg_colours)))

alg_colours <- read_csv("data/algorithm_colours.csv") %>% deframe()
# A list of algorithms to show coloured. The Rest will be grey
coloured_algs <- c("CSN_NCUA","PersonaDrive","PRODIGY","OncoImpact","sysSVM2","PhenoDriverR","DawnRank","SCS","PNC")

alg_colours[which(!names(alg_colours) %in% coloured_algs)] <- "grey"
plot_ns <- c(1,2,3,4,5,6,7,8,9,10)

ggplot(cos_similarity_plot %>% 
         filter(!str_detect(algorithm,"consensus")) %>% 
         mutate(start_label = if_else(n_drivers == min(n_drivers), as.character(algorithm), NA_character_),
                end_label = if_else(n_drivers == max(n_drivers), as.character(algorithm), NA_character_)) %>%
         filter(n_drivers %in% plot_ns), 
       aes(x=as.factor(n_drivers),
           y=pos, 
           colour=algorithm, 
           group=algorithm
       )) +
  geom_point(size=1.5) +
  geom_line(size=1,alpha=0.75) +
  scale_colour_manual(values = alg_colours) +
  ylab("Cosine Similarity (Dimensionality Reduced)") +
  xlab("Top N Drivers") +
  geom_label_repel(aes(label=start_label),
                   fill = "black",
                   size=3,
                   nudge_x = -3,
                   na.rm = T, direction = "y", arrow = arrow(type = "open", angle = 10, length = unit(0.3,"cm"))
  ) +
  #scale_fill_manual(values = "grey") +
  #geom_label_repel(aes(label=end_label),
  #                 nudge_x = 5,
  #                 na.rm = T) +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_line(colour="light grey"), 
        panel.grid.minor = element_blank(),
        #legend.position = "none",
        #panel.background = element_rect(fill="#1B2631")
  ) +
  guides(colour="none") +
  geom_vline(xintercept = 1, colour = "black", linetype = "dashed") +
  theme(text = element_text(size=10))

ggsave(paste0("plots/",run_mode,"/divergence_plot.png"))

}







