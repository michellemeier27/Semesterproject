##HALLMARKS PATHWAY SCORING 
#loading libraries
library(singscore)
library(qusage)
library(ggplot2)
library(GSVA)

#loading all source files
source("/Users/michellemeier/Semesterproject/ARCHS4/singscore_function.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/GSVA_function.R")
#not sure if i can just run this but that would be so lit
source("/Users/michellemeier/Semesterproject/ARCHS4/hypoxia/runfile_notch.R")


#get gene sets for both method types
gene_sets <- read.gmt("/Users/michellemeier/Semesterproject/ARCHS4/raw_data/gene_sets/hallmark_gene_sets.gmt")
geneSet <- getGmt("/Users/michellemeier/Semesterproject/ARCHS4/raw_data/gene_sets/hallmark_gene_sets.gmt")

#only keep the series with case and control 
notch_res_5 = notch_res_4[notch_res_4$case == TRUE,]

#only keep the series with direction case 
notch_res_5 = notch_res_5[notch_res_5$direction == "case",]

#make result df beforehand
#length gene set (= number of pathways)
geneSet <- getGmt("/Users/michellemeier/Semesterproject/ARCHS4/raw_data/gene_sets/hallmark_gene_sets.gmt")
l_gs <- length(geneSet) 
result_df <- data.frame(pathway = integer(l_gs),delta = integer(l_gs))


#all methods
#gsva
GSVA_rank <- findRankPathway(meta_frame = notch_res_5, expression_frame = expression_notch_cutoff, geneset = gene_sets, pathway_wanted = "HALLMARK_NOTCH_SIGNALING",
                             method_wanted = "gsva")
#ssGSEA
ssGSEA_rank <- findRankPathway(meta_frame = notch_res_5, expression_frame = expression_notch_cutoff, geneset = gene_sets, pathway_wanted = "HALLMARK_NOTCH_SIGNALING",
                               method_wanted = "ssgsea")
#zscore
zscore_rank <- findRankPathway(meta_frame = notch_res_5, expression_frame = expression_notch_cutoff, geneset = gene_sets, pathway_wanted = "HALLMARK_NOTCH_SIGNALING",
                               method_wanted = "zscore")
#PLAGE
plage_rank <- findRankPathway(meta_frame = notch_res_5, expression_frame = expression_notch_cutoff, geneset = gene_sets, pathway_wanted = "HALLMARK_NOTCH_SIGNALING",
                              method_wanted = "plage")
#singscore
singscore_rank <- findRankPathwaySingscore(meta_frame = notch_res_5, expression_frame = expression_notch_cutoff, geneset = geneSet, pathway_wanted = "HALLMARK_NOTCH_SIGNALING")

#make result overview data frame
all_res <- data.frame(GSVA = as.vector(GSVA_rank),
                      ssGSEA = as.vector(ssGSEA_rank),
                      zscore = as.vector(zscore_rank),
                      plage = as.vector(plage_rank),
                      singscore = as.vector(singscore_rank),
                      stringsAsFactors=FALSE)
##plot boxplot 
bp_1 <- ggplot(data = stack(all_res), aes(x = ind, y= values)) + geom_boxplot()
bp_1 <- bp_1 + labs(title ="notch pathway score (case)",
                subtitle = "hallmark gene sets")
show(bp_1)

ggsave(plot = bp_1, filename = "case_hallmarks_notch.png", path = "/Users/michellemeier/Semesterproject/ARCHS4/results/pathway_scoring/hallmarks/case")

