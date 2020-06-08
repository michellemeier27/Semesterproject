##clustering

subset_193 <- subset(hypoxia_res_4, hypoxia_res_4$series == "GSE104193")
dim(subset_193)
samples = as.vector(subset_193$sample)
cond = as.vector(subset_193$condition)
expression_193_draft = expression_hypoxia_cutoff[,colnames(expression_hypoxia_cutoff) %in% samples]



transpose = t(expression_193)
#hierarchical clustering
clusters <- hclust(dist(transpose))
plot(clusters, main = "hclust_GSE104193", labels = cond)
clusterCut <- cutree(clusters, 2)
table(clusterCut, cond)
png(filename = "test.png")
dev.off()
#ggsave(plot = clust,filename = "test.png", path= "/Users/michellemeier/Semesterproject/ARCHS4/results/hypoxia/")

#kmeans
set.seed(2)
kmeans <- kmeans(transpose,2, nstart = 2)
kmeans
table(kmeans$cluster, cond)
plot(kmeans$centers)


hier_clustering <- function(expression_frame, data_frame, exp_string){
  unique_series = as.character(unique(data_frame$series))
  for (series_u in unique_series){
    df <- subset(data_frame, data_frame$series == series_u)
    dim <- dim(df)
    if (dim[1] <=3){
      next
    }
    samples_for_series_v = as.vector(df$sample)
    cond = as.vector(df$condition)
    title = paste(exp_string, "hclust",series_u)
    expression_subset = expression_frame[,colnames(expression_frame) %in% samples_for_series_v]
    transpose = t(expression_subset)
    clusters <- hclust(dist(transpose))
    clust <- plot(clusters, main = title, labels = cond)
    clusterCut <- cutree(clusters, 2)
    table(clusterCut, cond)
    #ggsave(plot = clust, filename = paste0(title, ".png"), path = paste0("/Users/michellemeier/Semesterproject/ARCHS4/results/", path_wanted))
  }
}


hier_clustering(expression_hypoxia, hypoxia_res_4, "Hypoxia")









