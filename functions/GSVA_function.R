##FUNCTION FOR GSVA
findRankPathway <- function(meta_frame, expression_frame, geneset, pathway_wanted, method_wanted, res_df = result_df){
  #get all series IDs
  all_pos <- NA
  unique_series <- unique(meta_frame$series)
  #loop for all series
  for (series_id in unique_series){
    #cases  
    #get samples for case and control 
    case_sub <- subset(meta_frame, meta_frame$condition == "case" & meta_frame$series == series_id)
    case_s <- as.vector(case_sub$sample)
    control_sub <- subset(meta_frame, meta_frame$condition == "control" & meta_frame$series == series_id)
    control_s <- as.vector(control_sub$sample)
    #getting all samples for that series for expression  
    all_samples <- union(case_s, control_s)
    expression_series <- expression_frame[, colnames(expression_frame) %in% all_samples]
    #gsva for both case and control
    gsva_series <- gsva(expression_series, geneset , mx.diff = TRUE, abs.ranking = TRUE, method = method_wanted) 
    #make data frame for deltas values between case and control 
    l_sets <- length(geneset)
    names <- names(geneset)
    for (i in 1:l_sets){
      res_df[[1]][i] <- names[i]
      #getting medians of case and control group 
      median_case <- median(gsva_series[i,colnames(gsva_series) %in% case_s])
      median_control <- median(gsva_series[i,colnames(gsva_series) %in% control_s])
      #case median minus control median -> 
      delta <- abs(median_case - median_control)
      res_df[[2]][i] <- delta 

    }
    #sort list ascending to p_value
    new <- res_df[order(res_df$delta),]
    res <- grep(pathway_wanted, new$pathway)
    a <- min(res)
    #print(a)
    all_pos <- c(all_pos, a)
  }
  all_pos <- na.omit(all_pos)
  return(all_pos)
}


##test with hypoxia 
res_hypoxia_gsva <- findRankPathway(meta_frame = hypoxia_res_5, expression_frame = expression_hypoxia_cutoff, geneset = gene_sets, pathway_wanted = "HALLMARK_HYPOXIA", method_wanted = "gsva")



