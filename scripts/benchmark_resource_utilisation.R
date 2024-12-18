#benchmark_resource_utilisation.R


# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
#suppressPackageStartupMessages (library(ggpubr, quietly = T))
suppressPackageStartupMessages (library(svglite, quietly = T))
suppressPackageStartupMessages (library(scales, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="ALL", 
              help="algorithms to include in comparison separated by spaces, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-c", "--cancertype"), type="character", default="all", 
              help="Cancer types to include in the analysis separated by semicolons, or 'ALL' (Default)", metavar ="Cancer Type")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
algorithms <- opt$algorithms
celltypes <- opt$cancertype

cell_counts <- fread(paste0("benchmark_data/network_",network_choice,"/summary_info/cell_counts.csv"))

if(topper(celltypes) != "ALL"){
  celltypes <- str_split_1(celltypes, pattern = ";")
}else{
  celltypes <- cell_counts$lineage
}

full_alg_list <- c(
  "DawnRank", "OncoImpact", "PNC", "SCS", "PhenoDriverR", "PRODIGY", "PersonaDrive", "combined_de_novo_methods", "sysSVM2"
)


if(algorithms != "ALL"){
  celltypes <- str_split_1(celltypes, pattern = ";")
} else{
  algorithms <- full_alg_list
}


# Metadata


if(celltypes[1] != "ALL"){
  cell_counts <- cell_counts %>% filter(lineage %in% celltypes)
}

# Stats

stats <- foreach(lin=celltypes, .combine = "rbind") %do% {
  
  foreach(alg=algorithms, .combine = "rbind") %do% {
    
    if(file.exists(paste0("log/",alg,"_",network_choice,"_",lin,"_stats.txt"))){
      indiv_time <- fread(paste0("log/",alg,"_",network_choice,"_",lin,"_stats.txt"), col.names = c("runtime_sec","peak_VM_KiB")) %>%
        mutate(cancer_type=lin,algorithm=alg)
    }
    
  }
  
  
  
}


alg_colours <- fread("data/algorithm_colours.csv") %>% deframe()

if(any(!unique(stats$algorithm) %in% names(alg_colours))){
  # get the colour used for the selected de novo network method
  de_novo_col <- alg_colours[names(alg_colours) %in% c("CSN_DFVS",     "CSN_MDS",     "CSN_MMS",      "CSN_NCUA",
                                                "LIONESS_DFVS", "LIONESS_MDS", "LIONESS_MMS",  "LIONESS_NCUA",
                                                "SPCC_DFVS",    "SPCC_MDS",    "SPCC_MMS",     "SPCC_NCUA",
                                                "SSN_DFVS",     "SSN_MDS",     "SSN_MMS",      "SSN_NCUA"
                                                )]
  names(de_novo_col) <- "combined_de_novo_methods"
  alg_colours <- append(alg_colours,de_novo_col)
}


stats <- stats %>%
  left_join(cell_counts, by = c("cancer_type"="lineage")) %>%
  mutate(peak_VM_MiB=peak_VM_KiB/1000) %>%
  dplyr::select(-peak_VM_KiB) %>%
  pivot_longer(cols = c(runtime_sec,peak_VM_MiB), names_to = "measure", values_to = "value") %>%
  mutate(algorithm=factor(algorithm, levels = names(alg_colours))) #%>%
  
  
  ## Excluding some odd outliers
  
  #filter(!(algorithm=="PhenoDriverR" & cancer_type %in% c("Bone","Bowel","Breast","Bladder_Urinary_Tract"))) %>%
  #filter(!(algorithm=="PersonaDrive" & cancer_type %in% c("Bowel","Uterus")))
  


plot_breaks <- 10^(-10:10)
plot_minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))


ggplot(stats, aes(x=cell_count, y=value, colour=algorithm)) +
  geom_line(size=1.5) +
  scale_colour_manual(values = alg_colours) +
  #scale_y_continuous(trans = "log10") +
  scale_y_log10(breaks=plot_breaks,minor_breaks=plot_minor_breaks) +
  facet_wrap(~measure, scales="free",
             strip.position = "left", 
             labeller = as_labeller(c(peak_VM_MiB = "Peak Virtual Memory Use (MiB)", runtime_sec = "Total Runtime (Seconds)") )) +
  ylab(NULL) +
  xlab("Cohort Size") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  guides(colour=guide_legend(title="Algorithm (Colour)")) +
  annotation_logticks(sides = "l")
  
ggsave("plots/resource_utilisation.png", width = 30, height = 20, units = "cm", dpi = 300)

  






