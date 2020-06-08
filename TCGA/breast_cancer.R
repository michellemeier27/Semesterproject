#breast cancer 
library(ggplot2)
library(jtools)
#breast cancer (early stages)
colnames_wanted <- c("age_at_index","days_to_birth","gender","race","treatment_type","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis")
breast <- temp_meta_score[temp_meta_score$tissue_or_organ_of_origin == "Breast, NOS" & temp_meta_score$ajcc_pathologic_stage %in% early_stages,colnames(temp_meta_score) %in% colnames_wanted]
#pi3k
breast_fit1 <- lm(pi3k~ . , data =breast)
summary(breast_fit)
name <- "Breast cancer (early stage)\npi3k (N=337)"
t <- plot_summs(breast_fit1, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" = "age_at_index","male (N=3)" = "gender","caucasian (N=251)" = "racewhite","black or african american(N=52)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(breast_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/breast_cancer/checked_models")

#p53
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "treatment_or_therapy", "ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stages","prior_malignancy")
breast <- temp_meta_score[temp_meta_score$tissue_or_organ_of_origin == "Breast, NOS" & temp_meta_score$ajcc_pathologic_stage %in% early_stages,colnames(temp_meta_score) %in% colnames_wanted]
breast_fit2 <- lm(p53~ . , data =breast)
summary(breast_fit)
name <- "Breast cancer  (early stage)\np53 (N=337)"
t <- plot_summs(breast_fit2, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" = "age_at_index","male (N=3)" = "gender","caucasian (N=251)" = "racewhite","black or african american(N=52)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/breast_cancer/checked_models")

#metastasis
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "treatment_or_therapy", "ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stages","prior_malignancy")

breast <- temp_meta_score[temp_meta_score$tissue_or_organ_of_origin == "Breast, NOS" & temp_meta_score$ajcc_pathologic_stage %in% early_stages,colnames(temp_meta_score) %in% colnames_wanted]
breast_fit3 <- lm(metastasis~ . , data =breast)
summary(breast_fit)
name <- "Breast cancer  (early stage)\nmetastasis (N=337)"
t <- plot_summs(breast_fit3, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" = "age_at_index","male (N=3)" = "gender","caucasian (N=251)" = "racewhite","black or african american(N=52)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/breast_cancer/checked_models")

#notch
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "treatment_or_therapy", "ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stages","prior_malignancy")

breast <- temp_meta_score[temp_meta_score$tissue_or_organ_of_origin == "Breast, NOS" & temp_meta_score$ajcc_pathologic_stage %in% early_stages,colnames(temp_meta_score) %in% colnames_wanted]
breast_fit4 <- lm(notch~ . , data =breast)
summary(breast_fit)
name <- "Breast cancer  (early stage)\nnotch (N=337)"
t <- plot_summs(breast_fit4, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" = "age_at_index","male (N=3)" = "gender","caucasian (N=251)" = "racewhite","black or african american(N=52)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(breast_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/breast_cancer/checked_models")
#hypoxia
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "treatment_or_therapy", "ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stages","prior_malignancy")

breast <- temp_meta_score[temp_meta_score$tissue_or_organ_of_origin == "Breast, NOS" & temp_meta_score$ajcc_pathologic_stage %in% early_stages,colnames(temp_meta_score) %in% colnames_wanted]
breast_fit5 <- lm(hypoxia~ . , data =breast)
summary(breast_fit5)
name <- "Breast cancer  (early stage)\nhypoxia (N=337)"
t <- plot_summs(breast_fit5, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" = "age_at_index","male (N=3)" = "gender","caucasian (N=251)" = "racewhite","black or african american(N=52)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(breast_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/breast_cancer/checked_models")
#cell_cycle
colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "treatment_or_therapy", "ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stages","prior_malignancy")

breast <- temp_meta_score[temp_meta_score$tissue_or_organ_of_origin == "Breast, NOS" & temp_meta_score$ajcc_pathologic_stage %in% early_stages,colnames(temp_meta_score) %in% colnames_wanted]
breast_fit6 <- lm(cell_cycle~ . , data =breast)
summary(breast_fit)
name <- "Breast cancer  (early stage)\ncell_cycle (N=337)"
t <- plot_summs(breast_fit6, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" = "age_at_index","male (N=3)" = "gender","caucasian (N=251)" = "racewhite","black or african american(N=52)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(breast_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/breast_cancer/checked_models")


export_summs(breast_fit1, breast_fit2, breast_fit3,breast_fit4,breast_fit5,breast_fit6, scale = T, to.file = "docx", file.name = "BRCA_early_coefficients.docx",
             statistics = summ.lm, model.names = pathways)

#breast cancer (late stage)

colnames_wanted <- c("age_at_index","days_to_birth","race","treatment_type","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis")

breast <- temp_meta_score[temp_meta_score$tissue_or_organ_of_origin == "Breast, NOS" & temp_meta_score$ajcc_pathologic_stage %in% late_stages, colnames(temp_meta_score) %in% colnames_wanted]
#pi3k
breast_fit1 <- lm(pi3k~ . , data =breast)
summary(breast_fit)
name <- "Breast cancer (late stage)\npi3k (N=115)"
t <- plot_summs(breast_fit1, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" = "age_at_index","caucasian (N=90)" = "racewhite","black or african american(N=11)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(breast_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/breast_cancer/checked_models")
#p53
colnames_wanted <- c("age_at_index","days_to_birth","race","treatment_type","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis")

breast <- temp_meta_score[temp_meta_score$tissue_or_organ_of_origin == "Breast, NOS" & temp_meta_score$ajcc_pathologic_stage %in% late_stages, colnames(temp_meta_score) %in% colnames_wanted]
breast_fit2 <- lm(p53~ . , data =breast)
summary(breast_fit)
name <- "Breast cancer  (late stage)\np53 (N=115)"
t <- plot_summs(breast_fit2, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" = "age_at_index","caucasian (N=90)" = "racewhite","black or african american(N=11)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(breast_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/breast_cancer/checked_models")
#metastasis
colnames_wanted <- c("age_at_index","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "treatment_or_therapy", "ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stages","prior_malignancy")

breast <- temp_meta_score[temp_meta_score$tissue_or_organ_of_origin == "Breast, NOS" & temp_meta_score$ajcc_pathologic_stage %in% late_stages, colnames(temp_meta_score) %in% colnames_wanted]
breast_fit3 <- lm(metastasis~ . , data =breast)
summary(breast_fit)
name <- "Breast cancer  (late stage)\nmetastasis (N=115)"
t <- plot_summs(breast_fit3, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" = "age_at_index","caucasian (N=90)" = "racewhite","black or african american(N=11)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(breast_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/breast_cancer/checked_models")
#notch
colnames_wanted <- c("age_at_index","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "treatment_or_therapy", "ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stages","prior_malignancy")

breast <- temp_meta_score[temp_meta_score$tissue_or_organ_of_origin == "Breast, NOS" & temp_meta_score$ajcc_pathologic_stage %in% late_stages, colnames(temp_meta_score) %in% colnames_wanted]
breast_fit4 <- lm(notch~ . , data =breast)
summary(breast_fit)
name <- "Breast cancer  (late stage)\nnotch (N=115)"
t <- plot_summs(breast_fit4, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" = "age_at_index","caucasian (N=90)" = "racewhite","black or african american(N=11)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(breast_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/breast_cancer/checked_models")
#hypoxia
colnames_wanted <- c("age_at_index","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "treatment_or_therapy", "ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stages","prior_malignancy")

breast <- temp_meta_score[temp_meta_score$tissue_or_organ_of_origin == "Breast, NOS" & temp_meta_score$ajcc_pathologic_stage %in% late_stages, colnames(temp_meta_score) %in% colnames_wanted]
breast_fit5 <- lm(hypoxia~ . , data =breast)
summary(breast_fit5)
name <- "Breast cancer  (late stage)\nhypoxia (N=115)"
t <- plot_summs(breast_fit5, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" = "age_at_index","caucasian (N=90)" = "racewhite","black or african american(N=11)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
summary(breast_fit)
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/breast_cancer/checked_models")
#cell_cycle
colnames_wanted <- c("age_at_index","days_to_birth","race","treatment_type","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis")
breast <- temp_meta_score[temp_meta_score$tissue_or_organ_of_origin == "Breast, NOS" & temp_meta_score$ajcc_pathologic_stage %in% late_stages, colnames(temp_meta_score) %in% colnames_wanted]
breast_fit6 <- lm(cell_cycle~ . , data =breast)
summary(breast_fit)
name <- "Breast cancer  (late stage)\ncell_cycle (N=115)"
t <- plot_summs(breast_fit6, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("age" = "age_at_index","caucasian (N=90)" = "racewhite","black or african american(N=11)"= "raceblack or african american"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
summary(breast_fit)
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/breast_cancer/checked_models")



export_summs(breast_fit1, breast_fit2, breast_fit3,breast_fit4,breast_fit5,breast_fit6, scale = T, to.file = "docx", file.name = "BRCA_late_coefficients.docx",
             statistics = summ.lm, model.names = pathways)

