##FUNCTION SSGSEA (QC)
#load libraries
library(GSVA)
library(ggbiplot)

#define function (biplot and boxplot) for each individual series
ssGSEA <- function(gene_set_list, data_frame, expression_frame, exp_str, pathway){
  unique_series = as.character(unique(data_frame$series))
  for (series_u in unique_series){
    df <- subset(data_frame, data_frame$series == series_u)
    #only if there are both cases and controls
    if (df$case == "FALSE"){
      next
    }
    #for boxplot
    controls <- subset(df$sample, df$condition == "control")
    case <- subset(df$sample, df$condition == "case")
    #for gsva
    samples_for_series_v = as.vector(df$sample)
    if (length(samples_for_series_v) <= 4){
      ell = FALSE
    }
    else {
      ell = TRUE 
    }
    cond = as.vector(df$condition)
    title = paste("ssGSEA:",exp_str,series_u)
    expression_subset = expression_frame[,colnames(expression_frame) %in% samples_for_series_v]
    df_es <- gsva(expression_subset,gene_set_list, method = "ssgsea")
    transpose = t(df_es)
    res.pca <- prcomp(transpose, scale. = TRUE)
    g <- ggbiplot(res.pca,obs.scale = 1, var.scale = 1 , groups = cond, ellipse = ell,
                  circle = TRUE)
    g <- g + scale_color_discrete(name = 'conditions')
    g <- g + theme(legend.direction = 'vertical', 
                   legend.position = 'right')
    g <-  g + ggtitle(title)
    print(g)
    ggsave(filename = paste0(title, ".png"), path = paste0("/Users/michellemeier/Semesterproject/ARCHS4/results/", exp_str,"/"))
    #boxplot
    con = df_es[pathway, colnames(df_es) %in% controls]
    ca = df_es[pathway, colnames(df_es) %in% case]
    a <- data.frame(group = "control", value = con)
    b <- data.frame(group = "case", value = ca)
    plot.data <- rbind(a,b)
    f <- ggplot(plot.data, aes(x = group, y = value, fill = group))+geom_boxplot()
    title = paste(exp_str,"_ES_for_", series_u)
    f <-  f +ggtitle(title)
    ggsave(plot = f, filename =paste0(title,".png") , path =paste0("/Users/michellemeier/Semesterproject/ARCHS4/results/",exp_str,"/"))
  }
}






