read_driver_results <- function(mode,algorithms,cancertype){
  
  de_novo_collection <- c("CSN_DFVS",     "CSN_MDS",     "CSN_MMS",      "CSN_NCUA",
                          "LIONESS_DFVS", "LIONESS_MDS", "LIONESS_MMS",  "LIONESS_NCUA",
                          "SPCC_DFVS",    "SPCC_MDS",    "SPCC_MMS",     "SPCC_NCUA",
                          "SSN_DFVS",     "SSN_MDS",     "SSN_MMS",      "SSN_NCUA")
  
  suppressWarnings(rm(aggregated_results,specific_de_novo))
  for(alg in list.dirs(paste0("results/",mode,"/network_",network_choice), recursive = F)){
    
    
    # First, get names of algorithms that have been run
    alg = str_extract(alg,"(?<=/)[^/]*$")
    
    if(alg=="randomDrug"){next}
    
    # Check which algorithms are to be included
    
    # First check whether all are to be included, if so, then skip ahead to reading in files
    if(toupper(algorithms[1])!="ALL"){
      
      
      # Check if all of the names begin with "-" to indicate they are to be excluded
      if(all(substr(algorithms,1,1)=="-")){
        
        # If so, then skip any algorithms listed
        if(alg %in% gsub("-","",algorithms)){next}
        
      }else if(any(substr(algorithms,1,1)=="-")){
        warning("Invalid input for --algorithms. Cannot mix inclusions and exclusions in one statement.")
        stop() 
        # If the requested algorithm isn't present then skip it
      }else if(!alg %in% algorithms){
        
        # Unless it is a specific de novo method
        if(any(algorithms %in% de_novo_collection) & alg %in% de_novo_collection){
          specific_de_novo <- algorithms[algorithms %in% de_novo_collection]
          alg <- "combined_de_novo_methods"
        }else{
          message(paste0("Skipping ",alg,"..."))
          next
        }
      }
    }
    
    
    # read in consensus drivers
    
    if(alg=="consensus"){
      message(paste0("Reading result file for ", alg))
      indiv_result <- fread(paste0("results/",mode,"/network_",network_choice,"/consensus/consensus_drivers.csv")) %>%
        mutate(algorithm=alg)
      if("cell_ID" %in% colnames(indiv_result)){
        indiv_result <- indiv_result %>% dplyr::rename(sample_ID=cell_ID)
      }
      if("lineage" %in% colnames(indiv_result)){
        indiv_result <- indiv_result %>% dplyr::rename(cancer_type=lineage)
      }
      
      if(!exists("aggregated_results")){
        aggregated_results <- indiv_result
      }else{
        aggregated_results <- rbind(aggregated_results,indiv_result)
      }
    }else if(alg=="randomDriver"){
      # read in randomDriver results
      message(paste0("Reading result file for ", alg))
      for(result_dir in list.dirs(paste0("results/",mode,"/network_",network_choice,"/randomDriver"), recursive = F)){
        indiv_result <- foreach(result_file=list.files(result_dir), 
                                .combine = "rbind") %do% {
                                  fread(paste0(result_dir,"/",result_file))
                                }
      if(!exists("aggregated_results")){
        aggregated_results <- indiv_result
      }else{
        aggregated_results <- rbind(aggregated_results,indiv_result)
      }
      
      }
    }else{
      
    for(result_file in list.files(paste0("results/",mode,"/network_", network_choice,"/",alg), pattern = "*.csv", recursive = F)){
      
      
      if(!gsub(".csv","",result_file)%in%cancertype & toupper(cancertype)!= "ALL"){next}
      
      message(paste0("Reading result file for ", alg, "-", result_file))
      
      indiv_result <- data.table::fread(paste0("results/",mode,"/network_", network_choice,"/",alg,"/",result_file))
      
      if(!"algorithm" %in% colnames(indiv_result)){
        indiv_result <- indiv_result %>% mutate(algorithm=alg)
      }
      
      if("cell_ID" %in% colnames(indiv_result)){
        indiv_result <- indiv_result %>% dplyr::rename(sample_ID=cell_ID)
      }
      
      if("lineage" %in% colnames(indiv_result)){
        indiv_result <- indiv_result %>% dplyr::rename(cancer_type=lineage)
      }
      
      if(alg == "combined_de_novo_methods" & exists("specific_de_novo")){
        indiv_result <- indiv_result %>% filter(algorithm==specific_de_novo)
      }
      
      if(!exists("aggregated_results")){
        aggregated_results <- indiv_result
      }else{
        aggregated_results <- rbind(aggregated_results,indiv_result)
      }
    }
    }
  }
  
  return(aggregated_results)
  
}