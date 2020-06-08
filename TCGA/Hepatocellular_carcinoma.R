##Hepatocellular carcinoma, NOS
#libraries
library(ggplot2)
library(jtools)
#early stages
#liver (early stages)
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "treatment_or_therapy", "ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stages","prior_malignancy")

liver_early <- temp_meta_score[temp_meta_score$primary_diagnosis == "Hepatocellular carcinoma, NOS" & temp_meta_score$ajcc_pathologic_stage %in% early_stages ,]
summary(liver_early)


liver_early <- temp_meta_score[temp_meta_score$primary_diagnosis == "Hepatocellular carcinoma, NOS" & temp_meta_score$ajcc_pathologic_stage %in% early_stages , colnames(temp_meta_score) %in% colnames_wanted]
#pi3k
liver_early_fit1 <- lm(pi3k~ . , data =liver_early)
summary(liver_early_fit)
name <- "Hepatocellular carcinoma (early stage)\npi3k (N=91)"
t <- plot_summs(liver_early_fit1, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=67)" = "gender","caucasian (N=48)" = "racewhite","black or african american(N=3)"= "raceblack or african american", "not hispanic or latino (N=83)" = "ethnicity",
                                                                                                    "prior malignancy (N=5)" = "prior_malignancy"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(liver_early_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Hepatocellular_carcinoma/checked_model")
#p53
liver_early_fit2 <- lm(p53~ . , data =liver_early)
summary(liver_early_fit)

name <- "Hepatocellular carcinoma (early stage)\np53 (N=91)"
t <- plot_summs(liver_early_fit2, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=67)" = "gender","caucasian (N=48)" = "racewhite","black or african american(N=3)"= "raceblack or african american", "not hispanic or latino (N=83)" = "ethnicity",
                                                                                                    "prior malignancy (N=5)" = "prior_malignancy"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(liver_early_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Hepatocellular_carcinoma/checked_model")
#metastasis
liver_early_fit3 <- lm(metastasis~ . , data =liver_early)
summary(liver_early_fit)

name <- "Hepatocellular carcinoma (early stage)\nmetastasis (N=91)"
t <- plot_summs(liver_early_fit3, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs =  c("male (N=67)" = "gender","caucasian (N=48)" = "racewhite","black or african american(N=3)"= "raceblack or african american", "not hispanic or latino (N=83)" = "ethnicity",
                                                                                                     "prior malignancy (N=5)" = "prior_malignancy"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(liver_early_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Hepatocellular_carcinoma/checked_model")
#notch
liver_early_fit4 <- lm(notch~ . , data =liver_early)
summary(liver_early_fit)

name <- "Hepatocellular carcinoma (early stage)\nnotch (N=91)"
t <- plot_summs(liver_early_fit4, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=67)" = "gender","caucasian (N=48)" = "racewhite","black or african american(N=3)"= "raceblack or african american", "not hispanic or latino (N=83)" = "ethnicity",
                                                                                                    "prior malignancy (N=5)" = "prior_malignancy"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(liver_early_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Hepatocellular_carcinoma/checked_model")

#hypoxia
liver_early_fit5 <- lm(hypoxia~ . , data =liver_early)
summary(liver_early_fit)

name <- "Hepatocellular carcinoma (early stage)\nhypoxia (N=91)"
t <- plot_summs(liver_early_fit5, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=67)" = "gender","caucasian (N=48)" = "racewhite","black or african american(N=3)"= "raceblack or african american", "not hispanic or latino (N=83)" = "ethnicity",
                                                                                                    "prior malignancy (N=5)" = "prior_malignancy"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(liver_early_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Hepatocellular_carcinoma/checked_model")

#cell_cycle
liver_early_fit6 <- lm(cell_cycle~ . , data =liver_early)
summary(liver_early_fit)

name <- "Hepatocellular carcinoma (early stage)\ncell_cycle (N=91)"
t <- plot_summs(liver_early_fit6, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=67)" = "gender","caucasian (N=48)" = "racewhite","black or african american(N=3)"= "raceblack or african american", "not hispanic or latino (N=83)" = "ethnicity",
                                                                                                    "prior malignancy (N=5)" = "prior_malignancy"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(liver_early_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Hepatocellular_carcinoma/checked_model")



#liver late_stages
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "ajcc_pathologic_n","prior_malignancy")

liver_all <- temp_meta_score[temp_meta_score$primary_diagnosis == "Hepatocellular carcinoma, NOS" & temp_meta_score$ajcc_pathologic_stage %in% late_stages ,colnames(temp_meta_score) %in% colnames_wanted]

#pi3k
liver_all_fit1 <- lm(pi3k~ . , data =liver_all)
summary(liver_all_fit1)
name <- "Hepatocellular carcinoma (late_stage)\npi3k (N=36)"
t <- plot_summs(liver_all_fit1, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=24)" = "gender","caucasian (N=21)" = "racewhite","black or african american(N=1)"= "raceblack or african american", "not hispanic or latino (N=31)" = "ethnicity",
                                                                                                    "prior malignancy (N=2)" = "prior_malignancy"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Hepatocellular_carcinoma/checked_model")

#p53
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                      "ajcc_pathologic_m","ajcc_pathologic_n")
liver_all <- temp_meta_score[temp_meta_score$primary_diagnosis == "Hepatocellular carcinoma, NOS" & temp_meta_score$ajcc_pathologic_stage %in% late_stages ,colnames(temp_meta_score) %in% colnames_wanted]

liver_all_fit2 <- lm(p53~ . , data =liver_all)
summary(liver_all_fit2)

name <- "Hepatocellular carcinoma (late_stage)\np53 (N=36)"
t <- plot_summs(liver_all_fit2, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=24)" = "gender","caucasian (N=21)" = "racewhite","black or african american(N=1)"= "raceblack or african american", "not hispanic or latino (N=31)" = "ethnicity",
                                                                                                   "prior malignancy (N=2)" = "prior_malignancy"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Hepatocellular_carcinoma/checked_model")

#metastasis
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "treatment_or_therapy", "ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stages","prior_malignancy")
liver_all <- temp_meta_score[temp_meta_score$primary_diagnosis == "Hepatocellular carcinoma, NOS" & temp_meta_score$ajcc_pathologic_stage %in% late_stages ,colnames(temp_meta_score) %in% colnames_wanted]

liver_all_fit3 <- lm(metastasis~ . , data =liver_all)
summary(liver_all_fit3)

name <- "Hepatocellular carcinoma (late_stage)\nmetastasis (N=36)"
t <- plot_summs(liver_all_fit3, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=24)" = "gender","caucasian (N=21)" = "racewhite","black or african american(N=1)"= "raceblack or african american", "not hispanic or latino (N=31)" = "ethnicity",
                                                                                                   "prior malignancy (N=2)" = "prior_malignancy"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Hepatocellular_carcinoma/checked_model")

#notch
liver_all_fit4 <- lm(notch~ . , data =liver_all)
summary(liver_all_fit4)

name <- "Hepatocellular carcinoma (late_stage)\nnotch (N=36)"
t <- plot_summs(liver_all_fit4, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=24)" = "gender","caucasian (N=21)" = "racewhite","black or african american(N=1)"= "raceblack or african american", "not hispanic or latino (N=31)" = "ethnicity",
                                                                                                   "prior malignancy (N=2)" = "prior_malignancy"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Hepatocellular_carcinoma/checked_model")

#hypoxia
liver_all_fit5 <- lm(hypoxia~ . , data =liver_all)
summary(liver_all_fit5)

name <- "Hepatocellular carcinoma (late_stage)\nhypoxia (N=36)"
t <- plot_summs(liver_all_fit5, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=24)" = "gender","caucasian (N=21)" = "racewhite","black or african american(N=1)"= "raceblack or african american", "not hispanic or latino (N=31)" = "ethnicity",
                                                                                                   "prior malignancy (N=2)" = "prior_malignancy"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Hepatocellular_carcinoma/checked_model")

#cell_cycle
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "ajcc_pathologic_n","prior_malignancy")

liver_all <- temp_meta_score[temp_meta_score$primary_diagnosis == "Hepatocellular carcinoma, NOS" & temp_meta_score$ajcc_pathologic_stage %in% late_stages ,colnames(temp_meta_score) %in% colnames_wanted]

liver_all_fit6 <- lm(cell_cycle~ . , data =liver_all)
summary(liver_all_fit6)

name <- "Hepatocellular carcinoma (late_stage)\ncell_cycle (N=36)"
t <- plot_summs(liver_all_fit6, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=24)" = "gender","caucasian (N=21)" = "racewhite","black or african american(N=1)"= "raceblack or african american", "not hispanic or latino (N=31)" = "ethnicity",
                                                                                                   "prior malignancy (N=2)" = "prior_malignancy"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Hepatocellular_carcinoma/checked_model")


export_summs(liver_early_fit1, liver_early_fit2, liver_early_fit3,liver_early_fit4,liver_early_fit5,liver_early_fit6, scale = T, to.file = "docx", file.name = "LIHC_early_coefficients.docx",
             statistics = summ.lm, model.names = pathways)

export_summs(liver_all_fit1, liver_all_fit2, liver_all_fit3,liver_all_fit4,liver_all_fit5,liver_all_fit6, scale = T, to.file = "docx", file.name = "LIHC_late_coefficients.docx",
             statistics = summ.lm, model.names = pathways)


