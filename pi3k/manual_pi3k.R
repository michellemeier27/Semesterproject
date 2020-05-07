##loading all manual data into R for metastasis pathway
# -> download irrelevant samples from file
wrong_samples_pi3k_df <- read.csv(file = "/Users/michellemeier/Semesterproject/ARCHS4/pi3k/removed_series_samples_pi3k.csv", skip = 1)
wrong_samples_pi3k = as.vector(wrong_samples_pi3k_df$Sample)
blacklist_pi3k= as.vector(wrong_samples_pi3k_df$Series)

#downloading all manual case control annotation corrections
case_control_df <- read.csv(file = "/Users/michellemeier/Semesterproject/ARCHS4/pi3k/case_control_pi3k.csv", skip = 1)
case_list_pi3k = as.vector(case_control_df$Case)
control_list_pi3k = as.vector(case_control_df$Control)

#downloading all subseries
subseries_df <- read.csv(file= "/Users/michellemeier/Semesterproject/ARCHS4/pi3k/subseries_pi3k.csv", skip = 1)
subseries_1_pi3k = as.vector(subseries_df$Series1)
subseries_2_pi3k = as.vector(subseries_df$Series2)
subseries_3_pi3k = as.vector(subseries_df$Series3)
subseries_4_pi3k = as.vector(subseries_df$Series4)
subseries_5_pi3k = as.vector(subseries_df$Series5)
subseries_6_pi3k = as.vector(subseries_df$Series6)

#merging tech replicates 
tech_replicates_df <- read.csv(file= "/Users/michellemeier/Semesterproject/ARCHS4/pi3k/merging_tech_reps_pi3k.csv",na.strings = "", stringsAsFactors = FALSE)
merge_list_pi3k = as.list(tech_replicates_df)


#series that don't have a control 
no_control_pi3k = as.vector(case_control_df$no_case_control)

#gene set list for ssGSEA
#read gene sets from txt files
gs_pi3k <- read.table("raw_data/gene_sets/gene_set_pi3k.txt", sep = ",", skip = 3)
colnames(gs_pi3k) = c("genes")
gs_growth <- read.table("raw_data/gene_sets/gene_set_cell_growth.txt",sep = ",", skip = 3)
colnames(gs_growth)= c("genes")
gs_mtor <-read.table("raw_data/gene_sets/gene_set_mtor.txt",sep = ",", skip = 3)
colnames(gs_mtor) = c("genes")


gene_set_list_pi3k <- list(as.vector(gs_pi3k$genes), as.vector(gs_growth$genes), gs_mtor)
names(gene_set_list_pi3k) = c("pi3k" ,"cell_growth", "mtor")

#control direction
control_direction = c("GSE107707_3")









