##RUNSCRIPT
#load all libraries
library(rhdf5)
library(tidyverse)
library(Rtsne)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)
library("factoextra")
library(scone)
library(dplyr)
library(scater)
library(edgeR)
library(GSVA)
library(ggbiplot)


#load all source files
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/search_function.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/condition.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/expression_data.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/tSNE.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/hypoxia/hypoxia_manual_stuff.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/ignoring_finalising.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/plotting_series.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/PCA.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/normalise_scale.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/ssGSEA_function.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/merging_tech_replicates.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/lib_norm_check_p53.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/plotting_tsne_series_overview.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/case_control.R")



# query h5 file 
blacklist_hyp = c("GSE79029", "GSE86095", "GSE95640","GSE112866", "GSE129307", "GSE135487","GSE99987", "GSE125511","GSE10146", "GSE70767","GSE118046", "GSE113353","GSE98060",
                  "GSE69599", "GSE89838","GSE100146","GSE118362","GSE123571","GSE59987","GSE95280","GSE112137","GSE129344","GSE136908")
hyp_search_str =  "hypoxia|hif1|hypoxic|vhl"
hypoxia_res_1 = get_res_table( hyp_search_str, blacklist_hyp, samp_desc, extract_protocol,
                             sample_title, sample_char, sample_source, experiment_str="hypoxia")

#add the condition vector based on condition string
condition_str_contr = "wt|control|contr|normoxia|0h|normal|scrambled|normoxy|t0"
hypoxia_res_2 = condition(hypoxia_res_1, condition_str_contr)

#manually annotate data 
#correcting control/case annotation (control and case lists in manual script)
for (entry in control_list){
  print(entry)
  hypoxia_res_2[hypoxia_res_2$sample_accession== entry, 4] = "control"
}

for (entry in case_list){
  print(entry)
  hypoxia_res_2[hypoxia_res_2$sample_accession== entry, 4] = "case"
}

#seperate series with multiple cell lines into subseries 

for (entry in samples_1){
  hypoxia_res_2[hypoxia_res_2$sample_accession == entry , 8] = 2
  series = hypoxia_res_2[hypoxia_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_2")
  hypoxia_res_2$series_id <- addLevel(hypoxia_res_2$series_id, series_new)
  hypoxia_res_2[hypoxia_res_2$sample_accession == entry , 1] = series_new
}

for (entry in samples_2){
  hypoxia_res_2[hypoxia_res_2$sample_accession == entry , 8] = 3
  series = hypoxia_res_2[hypoxia_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_3")
  hypoxia_res_2$series_id <- addLevel(hypoxia_res_2$series_id, series_new)
  hypoxia_res_2[hypoxia_res_2$sample_accession == entry , 1] = series_new
}

for (entry in samples_3){
  hypoxia_res_2[hypoxia_res_2$sample_accession == entry , 8] = 4
  series = hypoxia_res_2[hypoxia_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_4")
  hypoxia_res_2$series_id <- addLevel(hypoxia_res_2$series_id, series_new)
  hypoxia_res_2[hypoxia_res_2$sample_accession == entry , 1] = series_new
}

#ignoring wrong samples
#read samples from generated csv file 
wrong_samples = read.csv("/Users/michellemeier/Semesterproject/ARCHS4/hypoxia/wrong_samples_hypoxia.csv")
wrong_samples = as.vector(wrong_samples[,1])
wrong_samples = c(wrong_samples, add_to_df)
hypoxia_res_3 <- hypoxia_res_2[!hypoxia_res_2$sample_accession %in% wrong_samples,]


#noting if there is a control for cases in each series
for (entry in list_of_series){
  hypoxia_res_3$one_control <- addLevel(hypoxia_res_3$one_control, "FALSE")
  hypoxia_res_3[hypoxia_res_3$series_id == entry , 7] = "FALSE"
}

#finalising 
hypoxia_res_4 = finalising_df(hypoxia_res_3)


#getting expression data
expression_hypoxia_raw = expression(hypoxia_res_4)

#loging hypoxia expression results for cutoff
expression_hypoxia_prep = log(expression_hypoxia_raw +1)

#cutoff for genes to only have significantly expressed genes
expression_hypoxia_cutoff = cutoff(expression_hypoxia_prep, hypoxia_res_4, cutoff_min = 8)

#lib size normalisation 
expression_hypoxia_norm = UQ_FN(expression_hypoxia_cutoff)

#scaling data before tSNE overview:
expression_hypoxia = normalise_scale(expression_hypoxia_norm, log = FALSE, CLR = FALSE)

# write tsv file. 
write_file(hypoxia_res_4, "results/hypoxia/FINAL/hypoxia_results.tsv")

# plots
tsne_hypoxia = tSNE(expression_hypoxia)


#plots
#plotting series (original) and series (subseries) overview
plotting_overview(tsne_hypoxia, hypoxia_res_4, "hypoxia")

#plotting all series tSNE
plotting_series(hypoxia_res_4, expression_hypoxia, "Hypoxia", "hypoxia/")

#pca for all series 
pca_series(hypoxia_res_4,expression_hypoxia, "Hypoxia", "hypoxia/")

##ssGSEA
#for all series individually: same data as for PCA and tSNE 
ssGSEA(gene_set_list_hypoxia, hypoxia_res_4, expression_hypoxia, "hypoxia",1)





