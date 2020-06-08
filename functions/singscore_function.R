##FUNCTION FOR SINGSCORE

#defining res_df
#length gene set (= number of pathways)
geneSet <- getGmt("/Users/michellemeier/Semesterproject/ARCHS4/raw_data/gene_sets/hallmark_gene_sets.gmt")
l_gs <- length(geneSet) 
result_df <- data.frame(pathway = integer(l_gs),delta = integer(l_gs))



findRankPathwaySingscore <- function(meta_frame, expression_frame, geneset, pathway_wanted, res_df = result_df){
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
    #rank gene expression
    ranked_expression <- rankGenes(expression_series)
    #getting multiscore for case and control
    score_singscore <- multiScore(ranked_expression, geneset)
    #make data frame for p deltas between case and control 
    l_sets <- length(geneset)
    names <- names(geneset)
    for (i in 1:l_sets){
      res_df[[1]][i] <- names[i]
      #getting medians of case and control group 
      u <- score_singscore[[1]]
      median_singscore_case <- median(u[i,colnames(u) %in% case_s])
      median_singscore_control <- median(u[i,colnames(u) %in% control_s])
      #case median minus control median -> 
      delta_singscore <- abs(median_singscore_case - median_singscore_control)
      res_df[[2]][i] <- delta_singscore
      

    }
    #sort list ascending to p_value
    new <- res_df[order(res_df$delta),]
    res <- grep(pathway_wanted, new$pathway)
    a <- min(res)
    all_pos <- c(all_pos, a)
  }
  all_pos <- na.omit(all_pos)

  return(all_pos)
}


##testing with hypoxia
#making sure only those with case and control

hypoxia_res_5 = hypoxia_res_4[hypoxia_res_4$case == TRUE,]

res_hypoxia <- findRankPathwaySingscore(meta_frame = hypoxia_res_5, expression_frame = expression_hypoxia_cutoff, geneset = geneSet, pathway_wanted = "HALLMARK_HYPOXIA")







