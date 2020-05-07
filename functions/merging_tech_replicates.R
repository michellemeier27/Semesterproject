##FUNCTION MERGING TECHNICAL REPLICATES

#defining function to merge expression values
merging_two_cols <- function(series_ids, expression_frame){
  exp = expression_frame
  for (i in series_ids){
    list = unlist(i)
    l = length(list)
    exp[,list[1]] <- rowMeans(exp[,list])
    exp = exp[,!colnames(exp) %in% list[2:l]]
  }
  return(exp)
}

#defining function to update meta data frame 
update <- function(series_ids, meta_df){
  ms = c()
  for (set in series_ids){
    list = unlist(set)
    l = length(list)
    ms = c(ms, list[2:l])
  }
  new_df <- meta_df[!meta_df$sample %in% ms,]
  return(new_df)
}





