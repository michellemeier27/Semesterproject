##PATHWAY SCORING PI3K
#adjust input data
#remove series of size 3
#find all series with only 3 samples
unique_series <- unique(pi3k_res_4$series)
v <- c("")
for (entry in unique_series){
  subset <- subset(pi3k_res_4, pi3k_res_4$series == entry)
  d <- dim(subset)
  if (d[1] < 4){
    v = c(v, entry)
  } 
}
#remove them
pi3k_res_5 <- pi3k_res_4[!pi3k_res_4$series %in% v,]

#only those with both case and control
pi3k_res_5 <- pi3k_res_5[pi3k_res_5$case == "TRUE",]
#only those with both case and control
pi3k_res_5 <- pi3k_res_5[pi3k_res_5$direction== "case",]

#generate random pathways
#generate true pathway
true_pi3k_gs <- list(gene_set_list_pi3k[[1]])
names(true_pi3k_gs) <- c("pi3k")
#length of wanted pathway (in this case pi3k)
l <- length(true_pi3k_gs[[1]])

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
unique_series <- unique(pi3k_res_5$series)
N <- length(unique_series) 
pathway_pos = data.frame(gsva = integer(N), ssGSEA = integer(N), zscore = integer(N), PLAGE = integer(N), singscore = integer(N))


##starting pathway function
start.time <- Sys.time()
##test
data_res_pi3k <- AllSeriesRank(pi3k_res_5, expression_pi3k_cutoff, true_pi3k_gs, gene_set_random)
end.time <- Sys.time()
time.taken_pi3k <- end.time - start.time
time.taken

# 
# ##plot boxplot 
# bp <- ggplot(data = stack(res), aes(x = ind, y= values)) + geom_boxplot()
# #title <- c("Hypoxia pathway score")
# bp <- bp + labs(title ="pi3k pathway score",
#                 subtitle = "100 random gene sets")
# show(bp)
# 
# ggsave(plot = bp, filename = "pi3k_all.png", path = "/Users/michellemeier/Semesterproject/ARCHS4/results/pathway_scoring/all_directions")
# 




