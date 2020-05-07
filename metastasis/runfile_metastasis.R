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
source("/Users/michellemeier/Semesterproject/ARCHS4/metastasis/manual_metastasis.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/ignoring_finalising.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/plotting_series.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/PCA.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/normalise_scale.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/ssGSEA_function.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/merging_tech_replicates.R")

##query
search_str_met = c("metastasis|catenin|nm23|brms1|kai1|klf17|gas1|sdpr|kiss1|neoangiogenesis|emt|e-cadherin|metastases")
met_res_1 = get_res_table( search_str_met, blacklist_met, samp_desc, extract_protocol,
                           sample_title, sample_char, sample_source, experiment_str="metastasis")


#add condition 
condition_str_contr_met = "control|contr|normal|scrambled|0h|none|non-treated|ctrl|ctl|ctr|wild-type|primary|parent"
met_res_2 = condition(met_res_1, condition_str_contr_met)

#manually annotate data 
#correcting control/case annotation (control and case lists in manual script)
for (entry in control_list_met){
  met_res_2[met_res_2$sample_accession== entry, 4] = "control"
}

for (entry in case_list_met){
  met_res_2[met_res_2$sample_accession== entry, 4] = "case"
}

#seperate series with multiple cell lines into subseries 
for (entry in subseries_1_met){
  met_res_2[met_res_2$sample_accession == entry , 8] = 1
  series = met_res_2[met_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_1")
  met_res_2$series_id <- addLevel(met_res_2$series_id, series_new)
  met_res_2[met_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_2_met){
  met_res_2[met_res_2$sample_accession == entry , 8] = 2
  series = met_res_2[met_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_2")
  met_res_2$series_id <- addLevel(met_res_2$series_id, series_new)
  met_res_2[met_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_3_met){
  met_res_2[met_res_2$sample_accession == entry , 8] = 3
  series = met_res_2[met_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_3")
  met_res_2$series_id <- addLevel(met_res_2$series_id, series_new)
  met_res_2[met_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_4_met){
  met_res_2[met_res_2$sample_accession == entry , 8] = 4
  series = met_res_2[met_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_4")
  met_res_2$series_id <- addLevel(met_res_2$series_id, series_new)
  met_res_2[met_res_2$sample_accession == entry , 1] = series_new
}


#ignoring wrong samples
met_res_3 <- met_res_2[!met_res_2$sample_accession %in% wrong_samples_met,]

#noting if there is a control for cases in each series
for (entry in no_controls){
  met_res_3$one_control <- addLevel(met_res_3$one_control, "FALSE")
  met_res_3[met_res_3$series_id == entry , 7] = "FALSE"
}

#noting directionality of enrichment score in ssGSEA
for (entry in control_direction_met){
  met_res_3$direction_ES <- addLevel(met_res_3$direction_ES, "control")
  met_res_3[met_res_3$series_id == entry, 9] = "control"
}

#finalising 
met_res_4 = finalising_df(met_res_3)
# write tsv file. 
write_file(met_res_4, "results/metastasis/FINAL/met_results.tsv")

##getting expression data
expression_met_raw_1 = expression(met_res_4)
#removing technical duplicates from data + get missing samples
expression_met_raw = merging_two_cols(merge_list_met, expression_met_raw_1)

#updating meta dataframe
met_res_4 = update(merge_list_met, met_res_4)

#log transform data for cutoff
expression_met_prep = log(expression_met_raw +1)

#filter
expression_met_cutoff = cutoff(expression_met_prep, met_res_4)

#lib size normalisation 
expression_met = UQ_FN(expression_met_cutoff)

#check normalisation 
#for all series
#unique = unique(met_res_4$series)
#for (series in unique){
  #boxplot_for_series(series, met_res_4, expression_met_cutoff,"not_norm","metastasis/")
#}

#for (series in unique){
  #boxplot_for_series(series, met_res_4, expression_met, "norm_uq","metastasis/")
#}

#scale for tSNE
expression_met_scale = normalise_scale(expression_met, CLR = FALSE, log = FALSE)

##tSNE
#overview: data log, clr, scale and cutoff
#calculations:
tsne_met = tSNE(expression_met_scale)
#plotting series (original) and series (subseries) overview
plotting_overview(tsne_met, met_res_4, "metastasis")
#tSNE for each series: data log, clr, scale and cutoff
plotting_series(met_res_4, expression_met_scale, "metastasis", "metastasis/")

##PCA
#for all series individually: data log, clr, scale and cutoff
pca_series(met_res_4, expression_met_scale, "metastasis", "metastasis/")


##ssGSEA
#for all series individually: same data as for PCA and tSNE 
ssGSEA(gene_set_list_met, met_res_4, expression_met_scale, "metastasis",1)



