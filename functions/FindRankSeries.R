##FUNCTION FOR PATHWAY SCORE
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/FindTruePathwayRank.R")

FindRanksSeries <- function(series_id,meta_frame, expression_frame, true_gene_set, randoms, random_df){
  #initialise all pos vector
  all_pos = NA
  #get samples for case and control 
  case_sub <- subset(meta_frame, meta_frame$condition == "case" & meta_frame$series == series_id)
  case_s <- as.vector(case_sub$sample)
  control_sub <- subset(meta_frame, meta_frame$condition == "control"  & meta_frame$series == series_id)
  control_s <- as.vector(control_sub$sample)
  #getting all samples for that series for expression  
  all_samples <- union(case_s, control_s)
  expression_series <- expression_frame[, colnames(expression_frame) %in% all_samples]
  
  #first find scores for true gene set using function FindTruePathwayRank
  res_true = FindTruePathway(case_s = case_s, control_s = control_s, true_gene_set = true_gene_set,expression_series = expression_series)
  
  #re-setting random_res
  random_res = random_df
  
  #looping over n numbers of random gene sets
  for (i in 1:length(randoms)){
    #getting one random gene set 
    name <- names(randoms)[i]
    random_geneSet <- GeneSet(randoms[[i]], setName = name)
    random_gs <- list(randoms[[i]])
    names(random_gs) <- name
    #name in df
    random_res[[1]][i] <- name
    #method = gsva
    gsva_series <- gsva(expression_series, random_gs , mx.diff = TRUE, abs.ranking = TRUE) 
    #getting medians of case and control group 
    median_gsva_case <- median(gsva_series[, colnames(gsva_series) %in% case_s])
    median_gsva_control <- median(gsva_series[, colnames(gsva_series) %in% control_s])
    #case median minus control median -> 
    delta_gsva <- abs(median_gsva_case - median_gsva_control)
    #in df
    random_res[[2]][i] = delta_gsva
    
    #method = ssGSEA
    ssgsea_series <-gsva(expression_series, random_gs , mx.diff = TRUE, abs.ranking = TRUE, method = "ssgsea") 
    #getting medians of case and control group 
    median_ssgsea_case <- median(ssgsea_series[, colnames(ssgsea_series) %in% case_s])
    median_ssgsea_control <- median(ssgsea_series[, colnames(ssgsea_series) %in% control_s])
    #case median minus control median -> 
    delta_ssgsea <- abs(median_ssgsea_case - median_ssgsea_control)
    #in df
    random_res[[3]][i] = delta_ssgsea
    
    #method = zscore
    zscore_series <- gsva(expression_series, random_gs , mx.diff = TRUE, abs.ranking = TRUE, method = "zscore") 
    #getting medians of case and control group 
    median_zscore_case <- median(zscore_series[, colnames(zscore_series) %in% case_s])
    median_zscore_control <- median(zscore_series[, colnames(zscore_series) %in% control_s])
    #case median minus control median -> 
    delta_zscore <- abs(median_zscore_case - median_zscore_control)
    #in df
    random_res[[4]][i] = delta_zscore
    
    #method = PLAGE
    plage_series <- gsva(expression_series, random_gs , mx.diff = TRUE, abs.ranking = TRUE, method = "plage" ) 
    #getting medians of case and control group 
    median_plage_case <- median(plage_series[, colnames(plage_series) %in% case_s])
    median_plage_control <- median(plage_series[, colnames(plage_series) %in% control_s])
    #case median minus control median -> 
    delta_plage <- abs(median_plage_case - median_plage_control)
    #in df
    random_res[[5]][i] = delta_plage
    
    #method = singscore
    rank_expression_series <- rankGenes(expression_series)
    singscore_series <- simpleScore(rank_expression_series, random_geneSet)
    #getting medians of case and control group 
    median_singscore_case <- median(singscore_series[1][rownames(singscore_series[1]) %in% case_s,])
    median_singscore_control <- median(singscore_series[1][rownames(singscore_series[1]) %in% control_s,])
    #case median minus control median -> 
    delta_singscore <- abs(median_singscore_case - median_singscore_control)
    #in df
    random_res[[6]][i] = delta_singscore
    
  }
  res_all = rbind(res_true, random_res)
  res_all <- na.omit(res_all)
  
  for (i in 2:6){
    ordered <- res_all[order(res_all[[i]]),]
    pos <- grep(names(true_gene_set), ordered$pathway)
    all_pos = c(all_pos, pos)
  }
  
  all_pos <- na.omit(all_pos)
  all_pos <- as.vector(all_pos)
  return(all_pos)
}


