##Pathway score notch

#adjust input data
#remove series of size 3
#find all series with only 3 samples
unique_series = unique(notch_res_4$series)
v <- c("")
for (entry in unique_series){
  subset <- subset(notch_res_4, notch_res_4$series == entry)
  d <- dim(subset)
  if (d[1] < 4){
    v = c(v, entry)
  } 
}
#remove them
notch_res_5 <- notch_res_4[!notch_res_4$series %in% v,]

#only those with both case and control
notch_res_5 <- notch_res_5[notch_res_5$case == "TRUE",]
notch_res_5 <- notch_res_5[notch_res_5$direction == "case",]

#generate random pathways
#generate true pathway
true_notch_gs <- list(gene_set_list_notch[[1]])
names(true_notch_gs) <- c("notch")
#length of wanted pathway (in this case pi3k)
l <- length(true_notch_gs[[1]])

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
unique_series <- unique(notch_res_5$series)
N <- length(unique_series) 
pathway_pos = data.frame(gsva = integer(N), ssGSEA = integer(N), zscore = integer(N), PLAGE = integer(N), singscore = integer(N))


##starting pathway function
start.time <- Sys.time()
##test
data_res_notch <- AllSeriesRank(notch_res_5, expression_notch_cutoff, true_notch_gs, gene_set_random)
end.time <- Sys.time()
time.taken_notch <- end.time - start.time
time.taken


# ##plot boxplot 
# bp <- ggplot(data = stack(res), aes(x = ind, y= values)) + geom_boxplot()
# bp <- bp + labs(title ="notch pathway score",
#                 subtitle = "100 random gene sets")
# show(bp)

ggsave(plot = bp, filename = "notch_all.png", path = "/Users/michellemeier/Semesterproject/ARCHS4/results/pathway_scoring/all_directions")

