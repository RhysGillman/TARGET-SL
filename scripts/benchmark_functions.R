# benchmark_functions.R


summarise_stats <- function(stats_df,alg_of_interest,min_gs,min_n) {
  n_samples <- length(unique(stats_df$sample_ID))
  stats_df_summarised <- stats_df %>%
    # Combining all badDriver simulations to one mean
    mutate(algorithm=ifelse(str_detect(algorithm,"randomDriver"), "randomDriver", algorithm)) %>%
    mutate(algorithm=ifelse(str_detect(algorithm,"randomDrug"), "randomDrug", algorithm)) %>%
    mutate(algorithm = ifelse(toupper(algorithm)=="CONSENSUS", "Consensus", algorithm)) %>%
    filter(toupper(algorithm) %in% toupper(alg_of_interest)) %>%
    # Only keeping results for samples with > min_gs gold standards
    filter(n_gs >= min_gs) %>%
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
           sample_size = ifelse(algorithm%in%c("randomDriver","randomDrug"),n_samples,sample_size))
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
  alg_linetype[alg_linetype=="randomDrug"] <- "dashed"
  alg_linetype[!alg_linetype%in%c("dotted","dashed")] <- "solid"
  
  p1 <- ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
    geom_line(data=plot_data %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "Consensus|random")), mapping=aes(alpha = sample_size), linewidth = 0.8) +
    geom_line(data=plot_data %>% filter(n<=N_max) %>% filter(algorithm=="Consensus"), linetype = "dotted", linewidth = 0.8) +
    geom_line(data=plot_data %>% filter(n<=N_max) %>% filter(algorithm%in%c("randomDriver","randomDrug")), linetype = "dashed", linewidth = 0.8) +
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
  
  
  p2 <- ggplot(plot_data %>% filter(!str_detect(algorithm, "Consensus|random")), aes(colour=sample_size, x=n,y=value)) +
    geom_line(alpha=0) +
    scale_color_gradient(high = "black",low = "white",
                         breaks = c(
                           min(plot_data %>% filter(!str_detect(algorithm, "Consensus|random")) %>% pull(sample_size)) %>% round(0),
                           max(plot_data %>% filter(!str_detect(algorithm, "Consensus|random")) %>% pull(sample_size)) %>% round(0)
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
  alg_linetype[alg_linetype=="randomDrug"] <- "dashed"
  alg_linetype[!alg_linetype%in%c("dotted","dashed")] <- "solid"
  
  ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
    geom_line(data=plot_data %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "Consensus|random")), linewidth = 0.8) +
    geom_line(data=plot_data %>% filter(n<=N_max) %>% filter(algorithm=="Consensus"), linetype = "dotted", linewidth = 0.8) +
    geom_line(data=plot_data %>% filter(n<=N_max) %>% filter(algorithm%in%c("randomDriver","randomDrug")), linetype = "dashed", linewidth = 0.8) +
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


top_n_plot_mean_CI <- function(title,stats_df,algs,N_max=10, measures="precision", comparisons=list(), y_pos=c(), y_pos_randomDriver){
  
  
  prediction_stats_trim <- stats_df %>%
    mutate(algorithm=ifelse(str_detect(algorithm,"randomDriver"), "randomDriver", algorithm)) %>%
    mutate(algorithm=ifelse(str_detect(algorithm,"randomDrug"), "randomDrug", algorithm)) %>%
    mutate(algorithm = ifelse(toupper(algorithm)=="CONSENSUS", "Consensus", algorithm)) %>%
    filter(toupper(algorithm) %in% toupper(algs)) %>%
    filter(n==N_max) %>%
    dplyr::select(algorithm,sample_ID,precision,recall,F1)
  
  
  pairwise_stats <- compare_means(precision~algorithm,data=prediction_stats_trim,p.adjust.method = "BH",method = "wilcox.test") %>%
    mutate(custom=ifelse(p.signif=="****", "<0.0001", "")) %>%
    mutate(custom=ifelse(p.signif=="***", "<0.001", custom)) %>%
    mutate(custom=ifelse(p.signif=="**", "<0.01", custom)) %>%
    mutate(custom=ifelse(p.signif=="*", "<0.05", custom)) %>%
    mutate(custom=ifelse(p.signif=="ns", "= ns",custom))
  
  vs_randomDriver <- pairwise_stats %>% filter(group1=="randomDriver"|group2=="randomDriver") %>% mutate(xpos=ifelse(group1=="randomDriver",group2,group1))

  specific_comparisons <- foreach(comp=comparisons, .combine = "rbind") %do% {
      pairwise_stats %>% filter(group1 %in% comp & group2 %in% comp)
    }
    
  
  
  plot_data <- prediction_stats_trim %>%
    group_by(algorithm) %>%
    summarise(
      n=n(),
      mean_precision=mean(precision, na.rm = T), 
      s_precision=sd(precision, na.rm=T),
      mean_recall=mean(recall, na.rm = T), 
      s_recall=sd(recall, na.rm=T),
      mean_F1=mean(F1, na.rm = T), 
      s_F1=sd(F1, na.rm=T)
    ) %>%
    ungroup() %>%
    mutate(
      ci95_precision=qt(0.975,df=n-2)*s_precision/sqrt(n),
      ci95_recall=qt(0.975,df=n-2)*s_recall/sqrt(n),
      ci95_F1=qt(0.975,df=n-2)*s_F1/sqrt(n)
    ) %>%
    mutate(algorithm = factor(algorithm, levels = names(alg_colours)))


  
  
  ggplot(plot_data, aes(x = algorithm, y = mean_precision, colour = algorithm)) +
    geom_point(position = position_dodge(width=0.75)) +
    geom_errorbar(aes(ymin=mean_precision-ci95_precision,ymax=mean_precision+ci95_precision),
                  position = position_dodge(width=0.75)) +
    scale_colour_manual(breaks = names(alg_colours),values = alg_colours) +
    {if(length(comparisons)>0)stat_pvalue_manual(specific_comparisons,label="p.signif", 
                       y.position = y_pos, color = "black", bracket.size = 0.5, tip.length = 0.01
    )} +
    stat_pvalue_manual(vs_randomDriver, x="xpos", y.position = y_pos_randomDriver, label = "p.signif") +
    ylab("Mean Precision +/- 95%CI") +
    xlab("") +
    theme_bw() +
    guides(colour="none") +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) +
    ggtitle(title) +
    coord_cartesian(ylim = c(0,NA))

}

traditional_benchmark <- function(predictions,gold_standard,n_max, type){
  
  
  
  if(toupper(type)=="DRUG"){
    predictions <- predictions %>% dplyr::rename(final_rank=drug_rank)
  }
  
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
        
        if(toupper(type)=="DRUG"){
          predicted <- tmp2 %>% filter(final_rank <= n) %>% pull(drug_ID)
        }else if(toupper(type)=="GENE"){
          predicted <- tmp2 %>% filter(final_rank <= n) %>% pull(target)
        }else{
          warning("Incorrect option for 'type' variable, use 'gene' or 'drug'")
        }
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
                   n_gs = length(gs))
        
        
      }
      
      
      
    }
    
  }
  
  return(result)
  
}
