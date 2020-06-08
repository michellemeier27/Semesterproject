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
source("/Users/michellemeier/Semesterproject/ARCHS4/cell_cycle/manual_cell_cycle.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/ignoring_finalising.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/plotting_series.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/PCA.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/normalise_scale.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/ssGSEA_function.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/merging_tech_replicates.R")

##query
search_str_ccp = c("cell cycle|cdkn|ccnd|ccne1|cdk2|cdk4|cdk6|rb1|e2f|p16|wee1|cdc20|apc/c|cdc25|ccne|ccna|ccnm")
ccp_res_1 = get_res_table( search_str_ccp, blacklist_ccp, samp_desc, extract_protocol,
                           sample_title, sample_char, sample_source, experiment_str="cell_cylce")

#add condition 
condition_str_contr_ccp = "control|contr|normal|scrambled|0h|none|non-treated|ctrl|ctl|ctr|wild-type|untreated"
ccp_res_2 = condition(ccp_res_1, condition_str_contr_ccp)

#manually annotate data 
#correcting control/case annotation (control and case lists in manual script)
for (entry in control_list_ccp){
  ccp_res_2[ccp_res_2$sample_accession== entry, 4] = "control"
}

for (entry in case_list_ccp){
  ccp_res_2[ccp_res_2$sample_accession== entry, 4] = "case"
}

#seperate series with multiple cell lines into subseries 
for (entry in subseries_1_ccp){
  ccp_res_2[ccp_res_2$sample_accession == entry , 8] = 1
  series = ccp_res_2[ccp_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_1")
  ccp_res_2$series_id <- addLevel(ccp_res_2$series_id, series_new)
  ccp_res_2[ccp_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_2_ccp){
  ccp_res_2[ccp_res_2$sample_accession == entry , 8] = 2
  series = ccp_res_2[ccp_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_2")
  ccp_res_2$series_id <- addLevel(ccp_res_2$series_id, series_new)
  ccp_res_2[ccp_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_3_ccp){
  ccp_res_2[ccp_res_2$sample_accession == entry , 8] = 3
  series = ccp_res_2[ccp_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_3")
  ccp_res_2$series_id <- addLevel(ccp_res_2$series_id, series_new)
  ccp_res_2[ccp_res_2$sample_accession == entry , 1] = series_new
}

#ignoring wrong samples
ccp_res_3 <- ccp_res_2[!ccp_res_2$sample_accession %in% wrong_samples_ccp,]

#noting if there is a control for cases in each series
for (entry in no_control_ccp){
  ccp_res_3$one_control <- addLevel(ccp_res_3$one_control, "FALSE")
  ccp_res_3[ccp_res_3$series_id == entry , 7] = "FALSE"
}

#noting directionality of enrichment score in ssGSEA
for (entry in directionality_ccp){
  ccp_res_3$direction_ES <- addLevel(ccp_res_3$direction_ES, "control")
  ccp_res_3[ccp_res_3$series_id == entry, 9] = "control"
}


#finalising 
ccp_res_4 = finalising_df(ccp_res_3)

# write tsv file. 
write_file(ccp_res_4, "results/cell_cycle/FINAL/ccp_results.tsv")

##getting expression data
expression_ccp_raw_1 = expression(ccp_res_4)
#removing technical duplicates from data
expression_ccp_raw = merging_two_cols(merge_list_ccp, expression_ccp_raw_1)

#updating meta dataframe
ccp_res_4 = update(merge_list_ccp, ccp_res_4)

#log transform data for cutoff
expression_ccp_prep = log(expression_ccp_raw +1)

#filter
expression_ccp_cutoff = cutoff(expression_ccp_prep, ccp_res_4, cutoff_min = 7)

#lib size normalisation 
expression_ccp = UQ_FN(expression_ccp_cutoff)

#check normalisation 
#for all series
#unique = unique(ccp_res_4$series)
#for (series in unique){
#boxplot_for_series(series, ccp_res_4, expression_ccp_cutoff,"not_norm","cell_cycle/")
#}

#for (series in unique){
#boxplot_for_series(series, ccp_res_4, expression_ccp, "norm_uq","cell_cycle/")
#}

#scale for tSNE
expression_ccp_scale = normalise_scale(expression_ccp, CLR = FALSE, log = FALSE)

##tSNE
#overview: data log, clr, scale and cutoff
#calculations:
tsne_ccp = tSNE(expression_ccp_scale)
#plotting series (original) and series (subseries) overview
plotting_overview(tsne_ccp, ccp_res_4, "cell_cycle")

#tSNE for each series: data log, clr, scale and cutoff
plotting_series(ccp_res_4, expression_ccp_scale, "cell_cycle", "cell_cycle/")

##PCA
#for all series individually: data log, clr, scale and cutoff
pca_series(ccp_res_4, expression_ccp_scale, "cell_cycle", "cell_cycle/")

##ssGSEA
#for all series individually: same data as for PCA and tSNE 
ssGSEA(gene_set_list_ccp, ccp_res_4, expression_ccp_scale, "cell_cycle",1)




