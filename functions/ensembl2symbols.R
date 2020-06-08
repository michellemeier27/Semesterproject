##CONVERTING ENSEMBL 2 SYMBOLS
#loading libraries
library(biomaRt)
library(org.Hs.eg.db)

##translate ensembl ids to gene names (because ARCHS4 sucks)
#read in file from Natalie with weights for true pathway
weights <- read.delim(file = "/Users/michellemeier/Semesterproject/ARCHS4/generated_data/natalie/sspaths_michelle_sspaths_weights.tsv")
gene_ids <- as.character(weights$gene_ids)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=gene_ids,mart= mart)
weights_symb <-merge(weights,G_list,by.x="gene_ids",by.y="ensembl_gene_id")
weights_symb <- weights_symb[c("gene_ids","pathway", "gene_weight","hgnc_symbol")]


#read in file from Natalie with weights for random pathway
random_weights <- read.csv("/Users/michellemeier/Semesterproject/ARCHS4/generated_data/random_gene_sets/random_gene_sets_weights.csv")
gene_ids_random <- random_weights$gene_ids
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
RG_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=gene_ids_random,mart= mart)
weights_symb_random <-merge(random_weights,RG_list,by.x="gene_ids",by.y="ensembl_gene_id")
weights_symb_random <- weights_symb[c("gene_ids","pathway", "gene_weight","hgnc_symbol")]

#read in file with genes to normalise
norm_genes <- read.csv(file = "/Users/michellemeier/Semesterproject/ARCHS4/generated_data/natalie/genes_to_normalise.csv")
genes <- as.character(norm_genes$genes)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
A_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
genes_to_normalise <- A_list$hgnc_symbol



# ##converting all symbols to ensembl -> doesn't work because weird ARCHS4 annotation... 
# gene_symbols_hypoxia <- rownames(expression_hypoxia_cutoff)
# ensembl_hypoxia <- mapIds(org.Hs.eg.db, gene_symbols_hypoxia, "ENSEMBL", "SYMBOL",multiVals = "list")


