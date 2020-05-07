##PCA FUNCTIONS
#loading libraries
library("factoextra")
library(DESeq2)

#defining PCA function for individual series
pca <- function(series_wanted, data_frame){
  subset <-subset(data_frame, data_frame$series == series_wanted)
  expression_subset = as.data.frame(expression(subset))
  transpose = t(expression_subset)
  res.pca <- prcomp(transpose, scale. = TRUE)
  res.var <- get_pca_var(res.pca)
  plot <-fviz_pca_ind(res.pca,repel = TRUE,habillage = subset$condition, invisible = "quali", ) + labs(title =series_wanted)
  show(plot)
  return(res.pca)
}

#defining PCA function for all individual series
pca_series <- function(data_frame,expression_frame, exp_string, path_wanted){
  unique_series = as.character(unique(data_frame$series))
  for (series_u in unique_series){
    df <- subset(data_frame, data_frame$series == series_u)
    samples_for_series_v = as.vector(df$sample)
    title = paste("PCA:",exp_string,series_u)
    expression_subset = expression_frame[,colnames(expression_frame) %in% samples_for_series_v]
    transpose = t(expression_subset)
    res.pca <- prcomp(transpose, scale. = FALSE)
    fviz_pca_ind(res.pca,repel = TRUE, habillage = df$condition, invisible = "quali") + labs(title =title)
    ggsave(filename = paste0(title,".png"), path = paste0("/Users/michellemeier/Semesterproject/ARCHS4/results/", path_wanted))
  }
}








