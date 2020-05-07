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
source("/Users/michellemeier/Semesterproject/ARCHS4/pi3k/manual_pi3k.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/ignoring_finalising.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/plotting_series.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/PCA.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/normalise_scale.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/ssGSEA_function.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/merging_tech_replicates.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/lib_norm_check_p53.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/plotting_tsne_series_overview.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/case_control.R")


##query
search_str_pi3k = c("pi3k|pten|pik3|inpp4b|tsc|stk11|rheb|ppp2r1a|rictor|mtor|wortmannin|ly294002|akt")
pi3k_res_1 = get_res_table( search_str_pi3k, blacklist_pi3k, samp_desc, extract_protocol,
                             sample_title, sample_char, sample_source, experiment_str="pi3k")

#add condition 
condition_str_contr_pi3k = "control|contr|normal|scrambled|0h|none|non-treated|ctrl|ctl|ctr|wild-type|untreated|wt"
pi3k_res_2 = condition(pi3k_res_1, condition_str_contr_pi3k)

#manually annotate data 
#correcting control/case annotation (control and case lists in manual script)
pi3k_res_2 = case_control(case_list_pi3k, control_list_pi3k, pi3k_res_2)

#seperate series with multiple cell lines into subseries 
#run source file because there are 6 subseries and that would be the most terrible mess 
source("/Users/michellemeier/Semesterproject/ARCHS4/pi3k/subseries_pi3k_annotation.R")

#ignoring wrong samples
pi3k_res_3 <- pi3k_res_2[!pi3k_res_2$sample_accession %in% wrong_samples_pi3k,]

#noting if there is a control for cases in each series
for (entry in no_control_pi3k){
  pi3k_res_3$one_control <- addLevel(pi3k_res_3$one_control, "FALSE")
  pi3k_res_3[pi3k_res_3$series_id == entry , 7] = "FALSE"
}

#noting directionality of enrichment score in ssGSEA
for (entry in control_direction){
  pi3k_res_3$direction_ES <- addLevel(pi3k_res_3$direction_ES, "control")
  pi3k_res_3[pi3k_res_3$series_id == entry, 9] = "control"
}

#finalising 
pi3k_res_4 = finalising_df(pi3k_res_3)
# write tsv file. 
write_file(pi3k_res_4, "results/pi3k/FINAL/pi3k_results.tsv")


##getting expression data
expression_pi3k_raw_1 = expression(pi3k_res_4)
#removing technical duplicates from data
expression_pi3k_raw = merging_two_cols(merge_list_pi3k, expression_pi3k_raw_1)

#updating meta dataframe
pi3k_res_4 = update(merge_list_pi3k, pi3k_res_4)

#log transform data for cutoff
expression_pi3k_prep = log(expression_pi3k_raw +1)

#filter
expression_pi3k_cutoff = cutoff(expression_pi3k_prep, pi3k_res_4, cutoff_min = 7)

#lib size normalisation
expression_pi3k = UQ_FN(expression_pi3k_cutoff)

#check normalisation
#for all series
# unique = unique(pi3k_res_4$series)
# for (series in unique){
# boxplot_for_series(series, pi3k_res_4, expression_pi3k_cutoff,"not_norm","pi3k/")
# }
# 
# for (series in unique){
# boxplot_for_series(series, pi3k_res_4, expression_pi3k, "norm_uq","pi3k/")
# }

#scale for tSNE
expression_pi3k_scale = normalise_scale(expression_pi3k, CLR = FALSE, log = FALSE)

##tSNE
#overview: data log, clr, scale and cutoff
#calculations:
tsne_pi3k = tSNE(expression_pi3k_scale)
#plotting series (original) and series (subseries) overview
plotting_overview(tsne_pi3k, pi3k_res_4, "pi3k")

#tSNE for each series: data log, clr, scale and cutoff
plotting_series(pi3k_res_4, expression_pi3k_scale, "pi3k", "pi3k/")

##PCA
#for all series individually: data log, clr, scale and cutoff
pca_series(pi3k_res_4, expression_pi3k_scale, "pi3k", "pi3k/")

##ssGSEA
#for all series individually: same data as for PCA and tSNE 
ssGSEA(gene_set_list_pi3k, pi3k_res_4, expression_pi3k_scale, "pi3k",1)

#reducing list to look at
list_to_check = c("GSE110610")
small_df <- pi3k_res_4[pi3k_res_4$series %in% list_to_check,]
#tsne
plotting_series(small_df,expression_pi3k_scale, "pi3k", "pi3k/")
#pca
pca_series(small_df,expression_pi3k_scale, "pi3k", "pi3k/")
#ssgsea
ssGSEA(gene_set_list_pi3k, small_df, expression_pi3k_scale, "pi3k",1)



