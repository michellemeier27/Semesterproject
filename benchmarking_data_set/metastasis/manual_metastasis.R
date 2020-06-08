##loading all manual data into R for metastasis pathway
# -> download irrelevant samples from file
wrong_samples_met_df <- read.csv(file = "/Users/michellemeier/Semesterproject/ARCHS4/metastasis/removed_samples_series_metastasis.csv", skip = 2)
wrong_samples_met = as.vector(wrong_samples_met_df$Samples)
blacklist_met= as.vector(wrong_samples_met_df$Series)


#downloading all manual case control annotation corrections
case_control_df <- read.csv(file = "/Users/michellemeier/Semesterproject/ARCHS4/metastasis/case_control_metastasis.csv", skip = 2)
case_list_met = as.vector(case_control_df$cases)
control_list_met = as.vector(case_control_df$controls)

#downloading all subseries
subseries_df <- read.csv(file= "/Users/michellemeier/Semesterproject/ARCHS4/metastasis/subseries_metastasis.csv", skip = 2)
subseries_1_met = as.vector(subseries_df$subseries.1)
subseries_2_met = as.vector(subseries_df$subseries.2)
subseries_3_met = as.vector(subseries_df$subseries.3)
subseries_4_met = as.vector(subseries_df$subseries.4)

#merging tech replicates 
tech_replicates_df <- read.csv(file= "/Users/michellemeier/Semesterproject/ARCHS4/metastasis/merging_tech_reps_metastasis.csv",na.strings = "", stringsAsFactors = FALSE)
merge_list_met = as.list(tech_replicates_df)
merge_list_met = lapply(merge_list_met, function(x) x[!is.na(x)])
merge_list_met = Filter(length, merge_list_met)

##gene set list for ssGSEA
#read gene sets from txt files 
gs_emt <- read.table("raw_data/gene_sets/gene_set_emt.txt", sep = ",", skip = 3)
colnames(gs_emt) = c("genes")
gs_focal <- read.table("raw_data/gene_sets/gene_set_focal_adhesion.txt", sep = ",", skip = 3)
colnames(gs_focal)= c("genes")
gs_a_junction <- read.table("raw_data/gene_sets/gene_set_adherens_junction.txt",sep = ",", skip = 3)
colnames(gs_a_junction)= c("genes")

gene_set_list_met <- list(as.vector(gs_emt$genes), as.vector(gs_focal$genes), as.vector(gs_a_junction$genes))
names(gene_set_list_met) = c("emt" ,"focal_adhesion","adherens_junction") 

no_controls <- c("GSE100066")
control_direction_met <- c("GSE59020","GSE63764","GSE71854_2","GSE74369","GSE80963_1","GSE80963_2","GSE80963_3","GSE80963_4",
                           "GSE110626_1","GSE112855_2")


