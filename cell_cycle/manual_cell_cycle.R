##manual file for cell cycle pathway
# -> download irrelevant samples from file
wrong_samples_ccp_df <- read.csv(file = "/Users/michellemeier/Semesterproject/ARCHS4/cell_cycle/removed_series_samples_ccp.csv", header = TRUE, skip = 1)
wrong_samples_ccp = as.vector(wrong_samples_ccp_df$Sample)
blacklist_ccp= as.vector(wrong_samples_ccp_df$Series)

#downloading all manual case control annotation corrections
case_control_df <- read.csv(file = "/Users/michellemeier/Semesterproject/ARCHS4/cell_cycle/case_control_ccp.csv", skip = 1)
case_list_ccp = as.vector(case_control_df$case)
control_list_ccp = as.vector(case_control_df$control)

#downloading all subseries
subseries_df <- read.csv(file= "/Users/michellemeier/Semesterproject/ARCHS4/cell_cycle/subseries_ccp.csv", skip = 1)
subseries_1_ccp = as.vector(subseries_df$subseries1)
subseries_2_ccp = as.vector(subseries_df$subseries2)
subseries_3_ccp = as.vector(subseries_df$subseries3)

#merging tech replicates 
tech_replicates_df <- read.csv(file= "/Users/michellemeier/Semesterproject/ARCHS4/cell_cycle/merging_tech_reps_ccp.csv",na.strings = "", stringsAsFactors = FALSE, header = TRUE)
merge_list_ccp= as.list(tech_replicates_df)
#merge_list_ccp = lapply(merge_list_met, function(x) x[!is.na(x)])
#merge_list_ccp = Filter(length, merge_list_met)

#series that don't have a control 
no_control_ccp = as.vector(subseries_df$no_control)

##gene set list for ssGSEA
#read gene sets from txt files 
gs_cell_cycle <- read.table("raw_data/gene_sets/gene_set_cell_cycle.txt", sep = ",", skip = 3)
colnames(gs_cell_cycle) = c("genes")
gs_g2m <- read.table("raw_data/gene_sets/gene_set_g2m_checkpoint.txt", sep = ",", skip = 3)
colnames(gs_g2m)= c("genes")
gs_mit_spindle <- read.table("raw_data/gene_sets/gene_set_mitotic_spindle.txt",sep = ",", skip = 3)
colnames(gs_mit_spindle)= c("genes")

gene_set_list_ccp <- list(as.vector(gs_cell_cycle$genes), as.vector(gs_g2m$genes), as.vector(gs_mit_spindle$genes))
names(gene_set_list_ccp) = c("cell_cycle" ,"g2m","mitotic_spindle") 

#directionality of enrichment score 
directionality_ccp = c("GSE74620","GSE77260","GSE84706","GSE94570_1","GSE94570_2","GSE94570_3","GSE101577","GSE109724_1",
                       "GSE109724_2","GSE121258_1","GSE122927_1","GSE122927_2","GSE133660")

