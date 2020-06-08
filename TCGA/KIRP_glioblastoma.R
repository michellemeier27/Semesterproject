##kidney cancers
library(ggplot2)
library(jtools)
library(tidyverse)

#kidney (early stages)
colnames_wanted <- c("age_at_index","days_to_birth","gender","race","treatment_type","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis")
kidney_early <- temp_meta_score[temp_meta_score$project_id == "TCGA-KIRP" & temp_meta_score$ajcc_pathologic_stage %in% early_stages, colnames(temp_meta_score) %in% colnames_wanted]
#pi3k
kidney_early_fit_1 <- lm(pi3k~. , data =kidney_early)
name <- "Renal papillary cell carcinoma (early stage)\npi3k (N=77)"
t <- plot_summs(kidney_early_fit_1, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" ="age_at_index","male (N=57)" = "gender","caucasian (N=56)" = "racewhite","black or african american(N=17)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Pappillary_adenocarcinoma(kidney)/checked_model")

#p53  
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "treatment_or_therapy")
kidney_early <- temp_meta_score[temp_meta_score$project_id == "TCGA-KIRP" & temp_meta_score$ajcc_pathologic_stage %in% early_stages, colnames(temp_meta_score) %in% colnames_wanted]
kidney_early_fit2 <- lm(p53~. , data =kidney_early)

name <- "Renal papillary cell carcinoma (early stage)\np53 (N=77)"
t <- plot_summs(kidney_early_fit2, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" ="age_at_index","male (N=57)" = "gender","caucasian (N=56)" = "racewhite","black or african american(N=17)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Pappillary_adenocarcinoma(kidney)/checked_model")

#metastasis
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis",
                     "treatment_or_therapy", "ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stages","prior_malignancy")
kidney_early <- temp_meta_score[temp_meta_score$project_id == "TCGA-KIRP" & temp_meta_score$ajcc_pathologic_stage %in% early_stages, colnames(temp_meta_score) %in% colnames_wanted]
kidney_early_fit3 <- lm(metastasis~ . , data =kidney_early)
name <- "Renal papillary cell carcinoma (early stage)\nmetastasis (N=77)"
t <- plot_summs(kidney_early_fit3, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" ="age_at_index","male (N=57)" = "gender","caucasian (N=56)" = "racewhite","black or african american(N=17)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Pappillary_adenocarcinoma(kidney)/checked_model")

#notch
colnames_wanted <- c("age_at_index","days_to_birth","gender","race","treatment_type","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis")
kidney_early <- temp_meta_score[temp_meta_score$project_id == "TCGA-KIRP" & temp_meta_score$ajcc_pathologic_stage %in% early_stages, colnames(temp_meta_score) %in% colnames_wanted]
kidney_early_fit4 <- lm(notch~ . , data =kidney_early)

name <- "Renal papillary cell carcinoma (early stage)\nnotch (N=77)"
t <- plot_summs(kidney_early_fit4, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" ="age_at_index","male (N=57)" = "gender","caucasian (N=56)" = "racewhite","black or african american(N=17)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Pappillary_adenocarcinoma(kidney)/checked_model")

#hypoxia
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis",
                     "treatment_or_therapy", "ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stages","prior_malignancy")
kidney_early <- temp_meta_score[temp_meta_score$project_id == "TCGA-KIRP" & temp_meta_score$ajcc_pathologic_stage %in% early_stages, colnames(temp_meta_score) %in% colnames_wanted]
kidney_early_fit5 <- lm(hypoxia~ . , data =kidney_early)

name <- "Renal papillary cell carcinoma (early stage)\nhypoxia (N=77)"
t <- plot_summs(kidney_early_fit5, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" = "age_at_index","male (N=57)" = "gender","caucasian (N=56)" = "racewhite","black or african american(N=17)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Pappillary_adenocarcinoma(kidney)/checked_model")
#cell_cycle
colnames_wanted <- c("age_at_index","days_to_birth","gender","race","treatment_type","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis")
kidney_early <- temp_meta_score[temp_meta_score$project_id == "TCGA-KIRP" & temp_meta_score$ajcc_pathologic_stage %in% early_stages, colnames(temp_meta_score) %in% colnames_wanted]
kidney_early_fit6 <- lm(cell_cycle~ . , data =kidney_early)

name <- "Renal papillary cell carcinoma (early stage)\ncell_cycle (N=77)"
t <- plot_summs(kidney_early_fit6, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" ="age_at_index","male (N=57)" = "gender","caucasian (N=56)" = "racewhite","black or african american(N=17)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Pappillary_adenocarcinoma(kidney)/checked_model")
pathways <- c("pi3k","p53","metastasis","notch","hypoxia","cell_cycle")
summ.lm = c(N_obs = "nobs",  adj.R2 = "adj.r.squared")
export_summs(kidney_early_fit_1, kidney_early_fit2, kidney_early_fit3,kidney_early_fit4,kidney_early_fit5,kidney_early_fit6, scale = T, to.file = "docx", file.name = "KIRP_coefficients.docx",
            statistics = summ.lm, model.names = pathways)


#glioblastoma manual
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis",
                     "treatment_or_therapy")
#pi3k
glio_early <- temp_meta_score[temp_meta_score$project_id == "TCGA-GBM",colnames(temp_meta_score) %in% colnames_wanted ]
glio_test <- glio_early[-c(13,17),] #try removing outlier from leverage plot 
glio_early_fit1 <- lm(pi3k~. , data =glio_early)
glio_test_fit1 <- lm(pi3k ~ ., data = glio_test)
name <- "Glioblastoma mulitforme (stage unknown)\npi3k (N=65)"
t <- plot_summs(glio_early_fit1, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=45)" = "gender","caucasian (N=60)" = "racewhite","black or african american(N=2)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Glioblastoma/checked_model")

#p53
glio_early_fit2 <- lm(p53~. , data =glio_early)
glio_early_test_fit2 <- lm(p53~. , data =glio_test)
summary(glio_early_test_fit2)
summary(glio_early_fit2)
name <- "Glioblastoma mulitforme (stage unknown)\np53 (N=65)"
t <- plot_summs(glio_early_fit2, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=45)" = "gender","caucasian (N=60)" = "racewhite","black or african american(N=2)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Glioblastoma/checked_model")

#metastasis
glio_early_fit3 <- lm(metastasis ~. , data =glio_early)
summary(glio_early_fit)
name <- "Glioblastoma mulitforme (stage unknown)\nmetastasis (N=65)"
t <- plot_summs(glio_early_fit3, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=45)" = "gender","caucasian (N=60)" = "racewhite","black or african american(N=2)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Glioblastoma/checked_model")

#notch
glio_early_fit4 <- lm(notch ~. , data =glio_early)
glio_early_test_fit4 <- lm(notch ~. , data =glio_test)
summary(glio_early_fit4)
summary(glio_early_test_fit4)

name <- "Glioblastoma mulitforme (stage unknown)\nnotch (N=65)"
t <- plot_summs(glio_early_fit4, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=45)" = "gender","caucasian (N=60)" = "racewhite","black or african american(N=2)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Glioblastoma/checked_model")

#hypoxia
glio_early_fit5 <- lm(hypoxia ~. , data =glio_early)
summary(glio_early_fit)
name <- "Glioblastoma mulitforme (stage unknown)\nhypoxia (N=65)"
t <- plot_summs(glio_early_fit5, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=45)" = "gender","caucasian (N=60)" = "racewhite","black or african american(N=2)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Glioblastoma/checked_model")

#cell cycle
glio_early_fit6 <- lm(cell_cycle ~. , data =glio_early)
summary(glio_early_fit)
name <- "Glioblastoma mulitforme (stage unknown)\ncell_cycle (N=65)"
t <- plot_summs(glio_early_fit6, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("male (N=45)" = "gender","caucasian (N=60)" = "racewhite","black or african american(N=2)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Glioblastoma/checked_model")

export_summs(glio_early_fit1, glio_early_fit2, glio_early_fit3,glio_early_fit4,glio_early_fit5,glio_early_fit6, scale = T, to.file = "docx", file.name = "GBM_coefficients.docx",
             statistics = summ.lm, model.names = pathways)
