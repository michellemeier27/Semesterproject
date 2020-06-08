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
source("/Users/michellemeier/Semesterproject/ARCHS4/notch/manual_notch.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/ignoring_finalising.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/plotting_series.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/PCA.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/normalise_scale.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/ssGSEA_function.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/merging_tech_replicates.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/lib_norm_check_p53.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/plotting_tsne_series_overview.R")


##query
search_str_notch = c("notch|creb|ep300|jag|cntn6|kdm5a|fbxw7|t-all|dll|hes1|nicd|lunar|mk0752|mk0752|ro4929097|mrk-003|mk-0752|pf03084014")
notch_res_1 = get_res_table( search_str_notch, blacklist_notch, samp_desc, extract_protocol,
                           sample_title, sample_char, sample_source, experiment_str="notch")

#add condition 
condition_str_contr_notch = "control|contr|normal|scrambled|0h|none|non-treated|ctrl|ctl|ctr|wild-type|untreated|wt"
notch_res_2 = condition(notch_res_1, condition_str_contr_notch)

#manually annotate data 
#correcting control/case annotation (control and case lists in manual script)
for (entry in control_list_notch){
  notch_res_2[notch_res_2$sample_accession== entry, 4] = "control"
}

for (entry in case_list_notch){
  notch_res_2[notch_res_2$sample_accession== entry, 4] = "case"
}

#seperate series with multiple cell lines into subseries 
for (entry in subseries_1_notch){
  notch_res_2[notch_res_2$sample_accession == entry , 8] = 1
  series = notch_res_2[notch_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_1")
  notch_res_2$series_id <- addLevel(notch_res_2$series_id, series_new)
  notch_res_2[notch_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_2_notch){
  notch_res_2[notch_res_2$sample_accession == entry , 8] = 2
  series = notch_res_2[notch_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_2")
  notch_res_2$series_id <- addLevel(notch_res_2$series_id, series_new)
  notch_res_2[notch_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_3_notch){
  notch_res_2[notch_res_2$sample_accession == entry , 8] = 3
  series = notch_res_2[notch_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_3")
  notch_res_2$series_id <- addLevel(notch_res_2$series_id, series_new)
  notch_res_2[notch_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_4_notch){
  notch_res_2[notch_res_2$sample_accession == entry , 8] = 4
  series = notch_res_2[notch_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_4")
  notch_res_2$series_id <- addLevel(notch_res_2$series_id, series_new)
  notch_res_2[notch_res_2$sample_accession == entry , 1] = series_new
}

#ignoring wrong samples
notch_res_3 <- notch_res_2[!notch_res_2$sample_accession %in% wrong_samples_notch,]

#noting if there is a control for cases in each series
for (entry in no_control_notch){
  notch_res_3$one_control <- addLevel(notch_res_3$one_control, "FALSE")
  notch_res_3[notch_res_3$series_id == entry , 7] = "FALSE"
}

#noting directionality of enrichment score in ssGSEA
for (entry in control_direction_notch){
  notch_res_3$direction_ES <- addLevel(notch_res_3$direction_ES, "control")
  notch_res_3[notch_res_3$series_id == entry, 9] = "control"
}

#finalising 
notch_res_4 = finalising_df(notch_res_3)
# write tsv file. 
write_file(notch_res_4, "results/notch/FINAL/notch_results.tsv")


##getting expression data
expression_notch_raw_1 = expression(notch_res_4)
#removing technical duplicates from data
expression_notch_raw = merging_two_cols(merge_list_notch, expression_notch_raw_1)

#updating meta dataframe
notch_res_4 = update(merge_list_notch, notch_res_4)

#log transform data for cutoff
expression_notch_prep = log(expression_notch_raw +1)

#filter
expression_notch_cutoff = cutoff(expression_notch_prep, notch_res_4, cutoff_min = 7)

#lib size normalisation 
expression_notch = UQ_FN(expression_notch_cutoff)

#check normalisation 
#for all series
# unique = unique(notch_res_4$series)
# for (series in unique){
# boxplot_for_series(series, notch_res_4, expression_notch_cutoff,"not_norm","notch/")
# }
# 
# for (series in unique){
# boxplot_for_series(series, notch_res_4, expression_notch, "norm_uq","notch/")
# }

#scale for tSNE
expression_notch_scale = normalise_scale(expression_notch, CLR = FALSE, log = FALSE)

##tSNE
#overview: data log, clr, scale and cutoff
#calculations:
tsne_notch = tSNE(expression_notch_scale)
#plotting series (original) and series (subseries) overview
plotting_overview(tsne_notch, notch_res_4, "notch")

#tSNE for each series: data log, clr, scale and cutoff
plotting_series(notch_res_4, expression_notch_scale, "notch", "notch/")

##PCA
#for all series individually: data log, clr, scale and cutoff
pca_series(notch_res_4, expression_notch_scale, "notch", "notch/")

##ssGSEA
#for all series individually: same data as for PCA and tSNE 
ssGSEA(gene_set_list_notch, notch_res_4, expression_notch_scale, "notch",1)


list_to_check = c("GSE59810","GSE77308","GSE97541")
small_df <- notch_res_4[notch_res_4$original_series %in% list_to_check,]
#tsne
plotting_series(small_df,expression_notch_scale, "notch", "notch/")
#pca
pca_series(small_df,expression_notch_scale, "notch", "notch/")
#ssgsea
ssGSEA(gene_set_list_notch, small_df, expression_notch_scale, "notch",1)





