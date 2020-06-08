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
source("/Users/michellemeier/Semesterproject/ARCHS4/hypoxia/runfile_cell_cycle.R")


#get gene sets for both method types
gene_sets <- read.gmt("/Users/michellemeier/Semesterproject/ARCHS4/raw_data/gene_sets/hallmark_gene_sets.gmt")
geneSet <- getGmt("/Users/michellemeier/Semesterproject/ARCHS4/raw_data/gene_sets/hallmark_gene_sets.gmt")

#only keep the series with case and control 
ccp_res_5 = ccp_res_4[ccp_res_4$case == TRUE,]

#only keep the series with direction case  
ccp_res_5 = ccp_res_5[ccp_res_5$direction == "control",]

#make result df beforehand
#length gene set (= number of pathways)
geneSet <- getGmt("/Users/michellemeier/Semesterproject/ARCHS4/raw_data/gene_sets/hallmark_gene_sets.gmt")
l_gs <- length(geneSet) 
result_df <- data.frame(pathway = integer(l_gs),delta = integer(l_gs))


#all methods
#gsva
GSVA_rank <- findRankPathway(meta_frame = ccp_res_5, expression_frame = expression_ccp_cutoff, geneset = gene_sets, pathway_wanted = "HALLMARK_G2M_CHECKPOINT",
                             method_wanted = "gsva")
#ssGSEA
ssGSEA_rank <- findRankPathway(meta_frame = ccp_res_5, expression_frame = expression_ccp_cutoff, geneset = gene_sets, pathway_wanted = "HALLMARK_G2M_CHECKPOINT",
                               method_wanted = "ssgsea")
#zscore
zscore_rank <- findRankPathway(meta_frame = ccp_res_5, expression_frame = expression_ccp_cutoff, geneset = gene_sets, pathway_wanted = "HALLMARK_G2M_CHECKPOINT",
                               method_wanted = "zscore")
#PLAGE
plage_rank <- findRankPathway(meta_frame = ccp_res_5, expression_frame = expression_ccp_cutoff, geneset = gene_sets, pathway_wanted = "HALLMARK_G2M_CHECKPOINT",
                              method_wanted = "plage")
#singscore
singscore_rank <- findRankPathwaySingscore(meta_frame = ccp_res_5, expression_frame = expression_ccp_cutoff, geneset = geneSet, pathway_wanted = "HALLMARK_G2M_CHECKPOINT")

#make result overview data frame
all_res <- data.frame(GSVA = as.vector(GSVA_rank),
                      ssGSEA = as.vector(ssGSEA_rank),
                      zscore = as.vector(zscore_rank),
                      plage = as.vector(plage_rank),
                      singscore = as.vector(singscore_rank),
                      stringsAsFactors=FALSE)
##plot boxplot 
bp <- ggplot(data = stack(all_res), aes(x = ind, y= values)) + geom_boxplot()
bp_1 <- bp + labs(title ="cell cycle pathway score (control)",
                subtitle = "hallmark gene sets")
show(bp_1)

ggsave(plot = bp_1, filename = "control_hallmarks_cell_cycle.png", path = "/Users/michellemeier/Semesterproject/ARCHS4/results/pathway_scoring/hallmarks/case")



t_test <- t.test(all_res$singscore, all_res$GSVA)
show(t_test$p.value)

# Compute the analysis of variance
my_data <- PlantGrowth
res.aov <- aov(weight ~ group, data = my_data)
summary(res.aov)
