##FUNCTION GET SCORE FOR TRUE PATHWAY

FindTruePathway <- function(case_s, control_s, true_gene_set, expression_series){
  #initiate data frame for all ranks
  curr_res = data.frame(pathway = names(true_gene_set), delta_gsva = integer(1),delta_ssgsea = integer(1), delta_zscore = integer(1), delta_plage = integer(1),   delta_singscore = integer(1))
  #true pathway GeneSet for singscore
  unlisted <- unlist(true_gene_set)
  true_geneSet <- GeneSet(unlisted, setName = names(true_gene_set))
  #go over all methods
  #method = gsva
  gsva_series <- gsva(expression_series, true_gene_set , mx.diff = TRUE, abs.ranking = TRUE) 
  #getting medians of case and control group 
  median_gsva_case <- median(gsva_series[, colnames(gsva_series) %in% case_s])
  median_gsva_control <- median(gsva_series[, colnames(gsva_series) %in% control_s])
  #case median minus control median -> 
  delta_gsva <- abs(median_gsva_case - median_gsva_control)
  #in df
  curr_res[[2]][1] = delta_gsva
  
  #method = ssGSEA
  ssgsea_series <-gsva(expression_series, true_gene_set , mx.diff = TRUE, abs.ranking = TRUE, method = "ssgsea") 
  #getting medians of case and control group 
  median_ssgsea_case <- median(ssgsea_series[, colnames(ssgsea_series) %in% case_s])
  median_ssgsea_control <- median(ssgsea_series[, colnames(ssgsea_series) %in% control_s])
  #case median minus control median -> 
  delta_ssgsea <- abs(median_ssgsea_case - median_ssgsea_control)
  #in df
  curr_res[[3]][1] = delta_ssgsea
  
  #method = zscore
  zscore_series <- gsva(expression_series, true_gene_set , mx.diff = TRUE, abs.ranking = TRUE, method = "zscore") 
  #getting medians of case and control group 
  median_zscore_case <- median(zscore_series[, colnames(zscore_series) %in% case_s])
  median_zscore_control <- median(zscore_series[, colnames(zscore_series) %in% control_s])
  #case median minus control median -> 
  delta_zscore <- abs(median_zscore_case - median_zscore_control)
  #in df
  curr_res[[4]][1] = delta_zscore
  
  #method = PLAGE
  plage_series <- gsva(expression_series, true_gene_set , mx.diff = TRUE, abs.ranking = TRUE, method = "plage" ) 
  #getting medians of case and control group 
  median_plage_case <- median(plage_series[, colnames(plage_series) %in% case_s])
  median_plage_control <- median(plage_series[, colnames(plage_series) %in% control_s])
  #case median minus control median -> 
  delta_plage <- abs(median_plage_case - median_plage_control)
  #in df
  curr_res[[5]][1] = delta_plage
  
  #method = singscore
  rank_expression_series <- rankGenes(expression_series)
  singscore_series <- simpleScore(rank_expression_series, true_geneSet)
  #getting medians of case and control group 
  median_singscore_case <- median(singscore_series[1][rownames(singscore_series[1]) %in% case_s,])
  median_singscore_control <- median(singscore_series[1][rownames(singscore_series[1]) %in% control_s,])
  #case median minus control median -> 
  delta_singscore <- abs(median_singscore_case - median_singscore_control)
  #in df
  curr_res[[6]][1] = delta_singscore
  
  return(curr_res)
}
