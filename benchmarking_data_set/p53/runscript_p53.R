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
source("/Users/michellemeier/Semesterproject/ARCHS4/p53/manual_p53.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/ignoring_finalising.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/plotting_series.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/PCA.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/normalise_scale.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/ssGSEA_function.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/merging_tech_replicates.R")

# query h5 file 
blacklist_p53 = c("GSE90013","GSE100292","GSE76656","GSE63234","GSE67109","GSE103934","GSE80609","GSE104729","GSE114502",
                  "GSE109326","GSE102796","GSE114850","GSE120741","GSE95057","GSE134441","GSE120396","GSE80467","GSE71876",
                  "GSE97356", "GSE110677","GSE101577","GSE129851","GSE118654","GSE99133","GSE81626","GSE108084","GSE118538",
                  "GSE132331","GSE132333","GSE109373","GSE79249","GSE94980","GSE99909","GSE120534","GSE68248","GSE81593",
                  "GSE81601","GSE95169","GSE110387","GSE123029","GSE119654","GSE86190","GSE58507")
#maybe GSE132331, GSE132333, GSE109373
p53_search_str =  c("p53|tp53|mdm2")
p53_search_str_drugs = c("cp-31398|prima-1|mira-1|rita|stictic acid|pk7088|chetomin|sch529074|wr-1065|sv40|mdm4") 
#-> gives one series hit that has nothing to do with p53...
p53_res_1 = get_res_table( p53_search_str, blacklist_p53, samp_desc, extract_protocol,
                               sample_title, sample_char, sample_source, experiment_str="p53")

#add condition 
condition_str_contr_p53 = "control|contr|normal|scrambled|0h|none|non-treated|ctrl|ctl|p53+/+|ctr|wild-type"
p53_res_2 = condition(p53_res_1, condition_str_contr_p53)

#manually annotate data 
#correcting control/case annotation (control and case lists in manual script)
for (entry in control_list_p53){
  print(entry)
  p53_res_2[p53_res_2$sample_accession== entry, 4] = "control"
}

for (entry in case_list_p53){
  print(entry)
  p53_res_2[p53_res_2$sample_accession== entry, 4] = "case"
}

#seperate series with multiple cell lines into subseries 
for (entry in subseries_1){
  p53_res_2[p53_res_2$sample_accession == entry , 8] = 2
  series = p53_res_2[p53_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_1")
  p53_res_2$series_id <- addLevel(p53_res_2$series_id, series_new)
  p53_res_2[p53_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_2){
  p53_res_2[p53_res_2$sample_accession == entry , 8] = 3
  series = p53_res_2[p53_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_2")
  p53_res_2$series_id <- addLevel(p53_res_2$series_id, series_new)
  p53_res_2[p53_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_3){
  p53_res_2[p53_res_2$sample_accession == entry , 8] = 4
  series = p53_res_2[p53_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_3")
  p53_res_2$series_id <- addLevel(p53_res_2$series_id, series_new)
  p53_res_2[p53_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_4){
  p53_res_2[p53_res_2$sample_accession == entry , 8] = 4
  series = p53_res_2[p53_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_4")
  p53_res_2$series_id <- addLevel(p53_res_2$series_id, series_new)
  p53_res_2[p53_res_2$sample_accession == entry , 1] = series_new
}

#ignoring wrong samples
p53_res_3 <- p53_res_2[!p53_res_2$sample_accession %in% wrong_samples_p53,]

#noting if there is a control for cases in each series
for (entry in no_control_p53){
  p53_res_3$one_control <- addLevel(p53_res_3$one_control, "FALSE")
  p53_res_3[p53_res_3$series_id == entry , 7] = "FALSE"
}

#noting directionality of enrichment score in ssGSEA
for (entry in case_direction){
  p53_res_3$direction_ES <- addLevel(p53_res_3$direction_ES, "case")
  p53_res_3[p53_res_3$series_id == entry, 9] = "case"
}
for (entry in control_direction_p53){
  p53_res_3$direction_ES <- addLevel(p53_res_3$direction_ES, "control")
  p53_res_3[p53_res_3$series_id == entry, 9] = "control"
}

#finalising 
p53_res_4 = finalising_df(p53_res_3)

# write tsv file. 
write_file(p53_res_4, "results/p53/FINAL/p53_results.tsv")


#getting expression data
expression_p53_raw_1 = expression(p53_res_4)

#removing technical duplicates from data + get missing samples
expression_p53_raw = merging_two_cols(merge_list, expression_p53_raw_1)

#updating meta dataframe
p53_res_4 = update(merge_list, p53_res_4)

# write tsv file. 
write_file(p53_res_4, "results/p53/FINAL/p53_results.tsv")

#log: not necessarily good in the beginning, but because i am only doing downstream analysis that needs
#log transforms is okay
expression_p53_prep = log(expression_p53_raw +1)

#cutoff for genes to only have significantly expressed genes
expression_p53_cutoff = cutoff(expression_p53_prep, p53_res_4, cutoff_min = 6)

#lib size normalisation 
expression_p53 = UQ_FN(expression_p53_cutoff)

#scaling data before tSNE overview:
expression_p53_scale = normalise_scale(expression_p53_cutoff, CLR = FALSE, log = FALSE)

##tSNE
#overview: data log, clr, scale and cutoff
#calculations:
tsne_p53 = tSNE(expression_p53_scale)
#plotting series (original) and series (subseries) overview
plotting_overview(tsne_p53, p53_res_4, "p53")
#tSNE for each series: data log, clr, scale and cutoff
plotting_series(p53_res_4, expression_p53_scale, "p53", "p53/")

##PCA
#for all series individually: data log, clr, scale and cutoff
pca_series(p53_res_4, expression_p53_scale, "p53", "p53/")

##ssGSEA
#for all series individually: same data as for PCA and tSNE 
ssGSEA(gene_set_list_p53, p53_res_4, expression_p53_scale, "p53",1)

#for all series individually: data lonly logged and cutoff (as recommended) ->
#not doing that because i want my data to have the same processing steps for all analysis 
ssGSEA(gene_set_list_p53, p53_res_4, expression_p53_cutoff, "p53")



