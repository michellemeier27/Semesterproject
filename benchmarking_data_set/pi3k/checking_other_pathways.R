##look at mtor signaling for one particular series GSE94980

#reducing list to look at
list_to_check = c("GSE94980")
small_df <- pi3k_res_4[pi3k_res_4$series %in% list_to_check,]

#ssgsea
#generate new gene set mtor
gs_mtor <- read.table("raw_data/gene_sets/gene_set_mtor.txt",sep = ",", skip = 3)
colnames(gs_mtor)= c("genes")
#get random sample set
genes_used <- as.vector(rownames(expression_pi3k_cutoff))
random <- sample(genes_used, 150)
gene_set_list_pi3k_mtor <- list(as.vector(gs_mtor$genes), random)
names(gene_set_list_pi3k_mtor) = c("mtor" , "random")

#plot
ssGSEA(gene_set_list_pi3k_mtor, small_df, expression_pi3k_scale, "pi3k",1)


