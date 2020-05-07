##loading all manual data into R for metastasis pathway
# -> download irrelevant samples from file
wrong_samples_notch_df <- read.csv(file = "/Users/michellemeier/Semesterproject/ARCHS4/notch/removed_series_sample_notch.csv", skip = 1)
wrong_samples_notch = as.vector(wrong_samples_notch_df$Samples)
blacklist_notch= as.vector(wrong_samples_notch_df$Series)

#downloading all manual case control annotation corrections
case_control_df <- read.csv(file = "/Users/michellemeier/Semesterproject/ARCHS4/notch/case_control_notch.csv", skip = 1)
case_list_notch = as.vector(case_control_df$case)
control_list_notch = as.vector(case_control_df$control)

#downloading all subseries
subseries_df <- read.csv(file= "/Users/michellemeier/Semesterproject/ARCHS4/notch/subseries_notch.csv", skip = 1)
subseries_1_notch = as.vector(subseries_df$Series1)
subseries_2_notch = as.vector(subseries_df$Series2)
subseries_3_notch = as.vector(subseries_df$Series3)
subseries_4_notch = as.vector(subseries_df$Series4)

#merging tech replicates 
tech_replicates_df <- read.csv(file= "/Users/michellemeier/Semesterproject/ARCHS4/notch/merging_tech_reps_notch.csv",na.strings = "", stringsAsFactors = FALSE)
merge_list_notch = as.list(tech_replicates_df)
merge_list_notch = lapply(merge_list_notch, function(x) x[!is.na(x)])
#merge_list_met = Filter(length, merge_list_met)

#series that don't have a control 
no_control_notch = as.vector(wrong_samples_notch_df$no_case_control)

#directionality 
control_direction_notch = c("GSE57982_2","GSE77108_1","GSE77108_2","GSE77108_3","GSE77108_4","GSE77308_1","GSE110142","GSE113753_1")

##gene set list for ssGSEA
#read gene sets from txt files
gs_notch <- read.table("raw_data/gene_sets/gene_set_notch.txt", sep = ",", skip = 3)
colnames(gs_notch) = c("genes")
gs_apoptosis <- read.table("raw_data/gene_sets/gene_set_apoptosis.txt", sep = ",", skip = 3)
colnames(gs_apoptosis)= c("genes")
gs_growth <- read.table("raw_data/gene_sets/gene_set_cell_growth.txt",sep = ",", skip = 3)
colnames(gs_growth)= c("genes")
gene_set_list_notch <- list(as.vector(gs_notch$genes), as.vector(gs_apoptosis$genes), as.vector(gs_growth$genes))
names(gene_set_list_notch) = c("notch" ,"apoptosis","cell_growth")







