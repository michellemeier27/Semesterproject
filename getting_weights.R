##weights for all pathways
source("/Users/michellemeier/Semesterproject/ARCHS4/ensembl2symbols.R")

#changing colnames for weights
colnames(weights_symb) <- c("gene_ids","pathway","gene_weight", "gene_id")

pathways <- c("hypoxia", "p53","pi3k","notch","cell_cycle","metastasis")

weights_list <- list()
for (entry in pathways){
  name <- paste0(entry,"_weights")
  weights <- subset(weights_symb[,-(1:2)], weights_symb$pathway == entry)
  weights_list[[name]] = weights
}
