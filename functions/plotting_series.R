##FUNCTION PLOTTING TSNE 
#libraries
library(Rtsne)
library(ggplot2)

#defining function for each series
plotting_series <- function(data_frame,expression_frame, exp_string, path_wanted){
  unique_series = as.character(unique(data_frame$series))
  for (series_u in unique_series){
    samples_for_series <- subset(data_frame, data_frame$series == series_u)
    samples_for_series_v = as.vector(samples_for_series$sample)
    dim = length(samples_for_series_v)
    if (dim < 4){
      next #sample sizes < 3 will not yield any clustering in tSNE
    }
    expression_subset = expression_frame[,colnames(expression_frame) %in% samples_for_series_v]
    series_tsne <- tSNE(expression_subset,1)
    new = as.data.frame(series_tsne$Y)
    title = paste(exp_string,series_u)
    plot <- ggplot(new, aes(x= V1, y= V2, colour= samples_for_series$condition, label= samples_for_series$sample)) + ggtitle(title) + geom_point(size =0.5) + geom_text(size = 2, vjust =0, nudge_y = 5) + scale_color_discrete("condition")
    ggsave(filename = paste0(title,".png"), path = paste0("/Users/michellemeier/Semesterproject/ARCHS4/results/", path_wanted))
    
  }
}


