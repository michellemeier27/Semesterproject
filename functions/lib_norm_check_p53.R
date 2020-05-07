##FUNCTION TO CHECK LIBRARY SIZE NORMALISATION
#loading all libraries
library(ggplot2)

#defining function for all series (boxplot)
boxplot_for_series <- function(series_id, meta_df, expression_frame, condition, path_wanted){
  subset <- subset(meta_df, meta_df$series == series_id)
  subset_samples = subset$sample
  expression_subset = expression_frame[,colnames(expression_frame) %in% subset_samples]
  expression_subset = as.data.frame(expression_subset)
  g <- ggplot(stack(expression_subset), aes(x = ind, y= values)) + geom_boxplot(outlier.color = "red", outlier.shape = 10)
  g <- g +ggtitle(paste("Boxplot series", series_id, condition))
  show(g)
  title = paste0("boxplot_",series_id, "_", condition,".png")
  ggsave(filename = title, path =paste0("/Users/michellemeier/Semesterproject/ARCHS4/results/",path_wanted)) 
}

##test with single series GSE74493
#1: data log transformed and cut off (no clr)
expression_log = normalise_scale(expression_p53_raw, CLR = FALSE, scale = FALSE)
expression_log_coff = cutoff(expression_log, p53_res_4, cutoff_min = 8)
boxplot_for_series("GSE94980", p53_res_4, expression_log_coff, "not_norm")

#2: data fully normalised (with CLR)
boxplot_for_series("GSE74493", p53_res_4, expression_p53,"norm_clr")

#3:data fully normalised (with UQ)
expression_uq = UQ_FN(expression_log)
boxplot_for_series("GSE74493", p53_res_4, expression_uq, "norm_uq")

#for all series
unique = unique(p53_res_4$series)
for (series in unique){
  boxplot_for_series(series, p53_res_4, expression_p53_cutoff,"norm_clr")
}

for (series in unique){
  boxplot_for_series(series, p53_res_4, expression_uq, "norm_uq")
}






