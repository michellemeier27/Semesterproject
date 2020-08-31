##Correlation study
#loading library
library(readr)
#read in files 
singscore_all_samples <- read_tsv(file = "/Users/michellemeier/Semesterproject/ARCHS4/generated_data/natalie/singscore_michelle_scores.tsv")
raw_meta <- read.csv(file = "/Users/michellemeier/Semesterproject/ARCHS4/raw_data/clinical.cases_selection/meta_data_tcga.csv", na.strings =c("'--","not reported"))
#only getting relevant columns
columns_wanted <- c("case_id","case_submitter_id",	"project_id",	"age_at_index","days_to_birth","ethnicity",	"gender","race","year_of_birth","age_at_diagnosis",
                    "ajcc_pathologic_m",	"ajcc_pathologic_n",	"ajcc_pathologic_stage",	"ajcc_pathologic_t","primary_diagnosis","prior_malignancy","prior_treatment",
                    "site_of_resection_or_biopsy","tissue_or_organ_of_origin","treatment_or_therapy",	"treatment_type")
meta_data_tcga <- raw_meta[,colnames(raw_meta) %in% columns_wanted]
#delete every second row because they come in double
toDelete <- seq(1, nrow(meta_data_tcga), 2)
meta_data_tcga_final <- meta_data_tcga[toDelete,]

#getting same tcga sample id for both data frames (adjusting singscore_tcga)
d <- dim(singscore_all_samples)
for (index in 1281:d[1]){
  new_id = character()
  a <- as.character(singscore_all_samples[index,1])
  b <- strsplit(a, "")[[1]]
  for (i in 1:12){
    new_id <- paste0(new_id, b[i])
  }
  singscore_all_samples[index,1] = new_id
}

#only tcga samples, not gtex
singscore_tcga <- singscore_all_samples[1281:d[1],]

#meta data for only the samples i have scores for 
temp_meta_score <- merge(singscore_tcga, meta_data_tcga_final, by.x = "sample_id", by.y = "case_submitter_id")




