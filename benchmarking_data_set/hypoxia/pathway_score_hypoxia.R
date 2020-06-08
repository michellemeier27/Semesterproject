##Pathway score p53

#adjust input data
#remove series of size 3
#find all series with only 3 samples
unique_series = unique(hypoxia_res_4$series)
v <- c("")
for (entry in unique_series){
  subset <- subset(hypoxia_res_4, hypoxia_res_4$series == entry)
  d <- dim(subset)
  if (d[1] < 4){
    v = c(v, entry)
  } 
}
#remove them
hypoxia_res_5 <- hypoxia_res_4[!hypoxia_res_4$series %in% v,]
#only those with both case and control
hypoxia_res_5 <- hypoxia_res_5[hypoxia_res_5$case == "TRUE",]
#only those in direction case
hypoxia_res_5 <- hypoxia_res_5[hypoxia_res_5$direction == "case",]

#generate random pathways
#generate true pathway
true_hypoxia_gs <- list(gene_set_list_hypoxia[[1]])
names(true_hypoxia_gs) <- c("hypoxia")
#length of wanted pathway 
l <- length(true_hypoxia_gs[[1]])

#get all gene sets and unlist
gene_sets <- read.gmt("/Users/michellemeier/Semesterproject/ARCHS4/raw_data/gene_sets/hallmark_gene_sets.gmt")
all_genes <- unique(unlist(gene_sets))
gene_set_random <- list()
for (i in 1:100){
  name <- paste0("random",i)
  random <- sample(all_genes, l)
  gene_set_random[[name]] <- random 
}

#generate random df
L <- length(gene_set_random)
random_res <- data.frame(pathway = integer(L), delta_gsva = integer(L),delta_ssgsea = integer(L), delta_zscore = integer(L), delta_plage = integer(L),   delta_singscore = integer(L))


#final data frame
unique_series <- unique(hypoxia_res_5$series)
N <- length(unique_series) 
pathway_pos = data.frame(gsva = integer(N), ssGSEA = integer(N), zscore = integer(N), PLAGE = integer(N), singscore = integer(N))


##starting pathway function
start.time <- Sys.time()
##test
data_res_hypoxia <- AllSeriesRank(hypoxia_res_5, expression_hypoxia_cutoff, true_hypoxia_gs, gene_set_random)
end.time <- Sys.time()
time.taken_hypoxia <- end.time - start.time


# #plot in the morning, see which one looks best
# ##plot boxplot 
# bp_hypoxia <- ggplot(data = stack(data_res_hypoxia), aes(x = ind, y= values)) + geom_boxplot()
# bp_hypoxia <- bp_hypoxia + labs(title ="hypoxia pathway score",
#                 subtitle = "100 random gene sets")
# show(bp_hypoxia)
# 
# ggsave(plot = bp_hypoxia, filename = "hypoxia_with_stat.png", path = "/Users/michellemeier/Semesterproject/ARCHS4/results/pathway_scoring/case_direction")
# 

