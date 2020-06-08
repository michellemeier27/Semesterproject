##FUNCTION FOR CELL CYCLE INHIBITOR PATHWAY SCORING
#loading all libraries
library(singscore)
library(ggplot2)
library(ggpubr)
library(GSVA)
#loading other source files
source("/Users/michellemeier/Semesterproject/ARCHS4/cell_cycle_inhibitors/runscript_drugs.R")
#defining genesets for singscore 
GeneSet1 <- GeneSet(unlist(true_ccp_gs), setName = "cell_cycle")
GeneSet2 <- GeneSet(unlist(true_hypoxia_gs), setName = "hypoxia")
GeneSet3 <- GeneSet(unlist(true_notch_gs), setName = "notch")
GeneSet4 <- GeneSet(unlist(true_met_gs), setName = "metastasis")
GeneSet5 <- GeneSet(unlist(true_p53_gs), setName = "p53")
GeneSetsCellCycleInhibitors <- GeneSetCollection(GeneSet1,GeneSet2,GeneSet3,GeneSet4,GeneSet5)
my_directory = c("/Users/michellemeier/Semesterproject/ARCHS4/results/cell_cycle_inhibitors")

list_avoid <- c("GSE99116", "GSE116187","GSE120298")
plus2 <- c("GSE116187")
plus4 <- c("GSE99116")
plabVSpalb <- summary_cell_cycle_inhibitors[!summary_cell_cycle_inhibitors$original_series %in% list_avoid,]
rest <- summary_cell_cycle_inhibitors[summary_cell_cycle_inhibitors$original_series %in% plus2,]
rest2 <- summary_cell_cycle_inhibitors[summary_cell_cycle_inhibitors$original_series %in% plus4,]

s <-summary_cell_cycle_inhibitors[summary_cell_cycle_inhibitors$original_series == "GSE130903",]



#defining function 
PathwayScoreAnalysis <- function(meta_frame=summary_cell_cycle_inhibitors, expression_frame=expression_drugs, pathway_sets=GeneSetsCellCycleInhibitors){
  #get all series
  unique_series <- unique(meta_frame$series)
  result <- list()
  k =0
  for (entry in unique_series){
    k = k+1
    sample_df <- meta_frame[meta_frame$series == entry,]
    sample_list <- sample_df$sample
    expression_series <- expression_frame[,colnames(expression_frame) %in% sample_list]
    ranked_series <- rankGenes(expression_series)
    res <- multiScore(ranked_series, pathway_sets)
    scores <- as.data.frame(t(res[[1]]))
    scores$samples <- rownames(scores)
    a <- merge(scores, sample_df, by.x = "samples", by.y = "sample")
    res <- a[,-c(8)]
    result[[k]] = res
  }
  return(result)
}


plottingPathwayScore <- function(pathway_scores =cell_cycle_inhibitors, directory= my_directory){
  for (i in 1:length(pathway_scores)){
    b <- pathway_scores[[i]]
    df <- data.frame(cell_cycle = b$cell_cycle,
                     hypoxia = b$hypoxia,
                     notch =b$notch,
                     metastasis = b$metastasis,
                     p53 = b$p53)
    res <- stack(df)
    colnames(res) <- c("scores","pathway")
    res$condition <- b$condition
    series <- unique(b$series)
    lb <- min(df) -  0.2
    ub <- max(df) + 0.3 * max(df)
    ly <- ub * 0.9
    #plotting
    p <- ggboxplot(res, x = "condition", y = "scores", color = "condition", palette = "jco", add = "jitter", facet.by = "pathway", short.panel.labs = T, ylim = c(lb, ub), xlab = "") +theme(axis.text=element_text(size=7))
    p <- p + rotate_x_text(angle = 45)+ ggtitle(series)+theme(legend.position= "bottom") +theme(axis.title.y= element_text(size=9))
    p <- p +theme(legend.text = element_text(size = 8)) +theme(legend.title = element_text(size = 8))+theme(plot.title=element_text(size=11))
    p <- p+ stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = ly, size = 3)
    ggsave(filename = paste0("pathway_score_", series,".png"), path = directory)
  }
}
  
p_s <- PathwayScoreAnalysis(meta_frame = s)
plottingPathwayScore(p_s)


#pathway scores 1
path_score <- PathwayScoreAnalysis(meta_frame = plabVSpalb)
plottingPathwayScore(path_score)
#pathway scores for the rest 
path_score_res <- PathwayScoreAnalysis(meta_frame = rest)
u_series <- unique(rest$series)


for (i in 1:length(path_score_res)){
  b <- path_score_res[[i]]
  df <- data.frame(cell_cycle = b$cell_cycle,
                   hypoxia = b$hypoxia,
                   notch =b$notch,
                   metastasis = b$metastasis,
                   p53 = b$p53)
  res <- stack(df)
  colnames(res) <- c("scores","pathway")
  res$condition <- b$condition
  series <- unique(b$series)
  lb <- min(df) -  0.2
  ub <- max(df) + 0.2
  ly <- ub * 0.9
  
  p <- ggboxplot(res, x = "condition", y = "scores", color = "condition", palette = "jco", add = "jitter", facet.by = "pathway", short.panel.labs = T, ylim = c(lb, ub), xlab = "") +theme(axis.text=element_text(size=7))
  p <- p + rotate_x_text(angle = 45)+ ggtitle(series)+theme(legend.position= "bottom") +theme(axis.title.y= element_text(size=9))
  p <- p +theme(legend.text = element_text(size = 8)) +theme(legend.title = element_text(size = 8))+theme(plot.title=element_text(size=11))
  ggsave(filename = paste0("pathway_score_", series,".png"), path = directory)
}

#zscore for pi3k 
#defining function 
PathwayZ_ScoreAnalysis <- function(meta_frame=plabVSpalb, expression_frame=expression_drugs, pathway_sets=true_pi3k_gs){
  #get all series
  unique_series <- unique(meta_frame$series)
  result <- list()
  k =0
  for (entry in unique_series){
    k = k+1
    sample_df <- meta_frame[meta_frame$series == entry,]
    sample_list <- sample_df$sample
    expression_series <- expression_frame[,colnames(expression_frame) %in% sample_list]
    res <- gsva(expression_series, pathway_sets, method = "zscore")
    scores <- as.data.frame(t(res))
    scores$samples <- rownames(scores)
    a <- merge(scores, sample_df, by.x = "samples", by.y = "sample")
    res <- a[,-c(8)]
    result[[k]] = res
  }
  return(result)
}
plottingPathwayZ_Score <- function(pathway_scores =result, directory= my_directory){
  for (i in 1:length(pathway_scores)){
    b <- pathway_scores[[i]]
    df <- data.frame(pi3k = b$pi3k)
    res <- stack(df)
    colnames(res) <- c("scores","pathway")
    res$condition <- b$condition
    series <- unique(b$series)
    lb <- min(df) -  0.2
    ub <- max(df) + 0.3 * max(df)
    ly <- ub * 0.9
    #plotting
    p <- ggboxplot(res, x = "condition", y = "scores", color = "condition", palette = "jco", add = "jitter", facet.by = "pathway", short.panel.labs = T, ylim = c(lb, ub), xlab = "") +theme(axis.text=element_text(size=7))
    p <- p + rotate_x_text(angle = 45)+ ggtitle(series)+theme(legend.position= "bottom") +theme(axis.title.y= element_text(size=9))
    p <- p +theme(legend.text = element_text(size = 8)) +theme(legend.title = element_text(size = 8))+theme(plot.title=element_text(size=11))
    p <- p+ stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = ly, size = 3)
    ggsave(filename = paste0("pathway_zscore_", series,".png"), path = directory)
  }
}
#for two group series 
result <- PathwayZ_ScoreAnalysis()
plottingPathwayZ_Score()
res_s <- PathwayZ_ScoreAnalysis(s)
plottingPathwayZ_Score(res_s)

#for 4 group series
result_4 <- PathwayZ_ScoreAnalysis(meta_frame = rest)
for (i in 1:length(result_4)){
  b <- result_4[[i]]
  df <- data.frame(pi3k = b$pi3k)
  res <- stack(df)
  colnames(res) <- c("scores","pathway")
  res$condition <- b$condition
  series <- unique(b$series)
  lb <- min(df) -  0.2
  ub <- max(df) + 0.5
  ly <- ub * 0.9
  
  p <- ggboxplot(res, x = "condition", y = "scores", color = "condition", palette = "jco", add = "jitter", facet.by = "pathway", short.panel.labs = T, ylim = c(lb, ub), xlab = "") +theme(axis.text=element_text(size=7))
  p <- p + rotate_x_text(angle = 45)+ ggtitle(series)+theme(legend.position= "bottom") +theme(axis.title.y= element_text(size=9))
  p <- p +theme(legend.text = element_text(size = 8)) +theme(legend.title = element_text(size = 8))+theme(plot.title=element_text(size=11))
  p <- p+ stat_compare_means(method = "t.test", label = "p.signif", label.y = ly, ref.group = "placebo", size = 3)
  ggsave(filename = paste0("pathway_zscore_", series,".png"), path = directory)
}

#for the ones with not enough samples 
result_not_enough <- PathwayZ_ScoreAnalysis(meta_frame = rest2)
for (i in 1:length(result_not_enough)){
  b <- result_not_enough[[i]]
  df <- data.frame(pi3k = b$pi3k)
  res <- stack(df)
  colnames(res) <- c("scores","pathway")
  res$condition <- b$condition
  series <- unique(b$series)
  lb <- min(df) -  0.2
  ub <- max(df) + 0.2
  ly <- ub * 0.9
  
  p <- ggboxplot(res, x = "condition", y = "scores", color = "condition", palette = "jco", add = "jitter", facet.by = "pathway", short.panel.labs = T, ylim = c(lb, ub), xlab = "") +theme(axis.text=element_text(size=7))
  p <- p + rotate_x_text(angle = 45)+ ggtitle(series)+theme(legend.position= "bottom") +theme(axis.title.y= element_text(size=9))
  p <- p +theme(legend.text = element_text(size = 8)) +theme(legend.title = element_text(size = 8))+theme(plot.title=element_text(size=11))
  p
  ggsave(filename = paste0("pathway_zscore_", series,".png"), path = directory)
}


