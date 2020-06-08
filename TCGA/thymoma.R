#project_id thym
library(ggplot2)
library(jtools)

colnames_wanted <- c("age_at_index","gender","race","hypoxia","notch", "p53","pi3k","cell_cycle","metastasis","ethnicity",
                     "treatment_or_therapy")
#pi3k
thym <- temp_meta_score[temp_meta_score$project_id == "TCGA-THYM", colnames(temp_meta_score) %in% colnames_wanted]
thym_fit1 <- lm(pi3k~ . , data =thym)
summary(thym_fit)
#n asian = 3
name <- paste0("Thymoma (stage unknown)\npi3k (N=38)")
t <- plot_summs(thym_fit1, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("hypoxia","notch","p53","metastasis","cell_cycle") )
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Thymoma/checked_models")

#p53
thym <- temp_meta_score[temp_meta_score$project_id == "TCGA-THYM", colnames(temp_meta_score) %in% colnames_wanted]
thym_fit2 <- lm(p53~ . , data =thym)
summary(thym_fit)
#n asian = 3
name <- paste0("Thymoma (stage unknown)\np53 (N=38)")
t <- plot_summs(thym_fit2, plot.distributions = T, scale = T, inner_ci_level = 0.9,coefs =  c("hypoxia","notch","pi3k","metastasis","cell_cycle"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Thymoma/checked_models")

#metastasis
thym <- temp_meta_score[temp_meta_score$project_id == "TCGA-THYM", colnames(temp_meta_score) %in% colnames_wanted]
thym_fit3 <- lm(metastasis~ . , data =thym)
summary(thym_fit)
#n asian = 3
name <- paste0("Thymoma (stage unknown)\nmetastasis (N=38)")
t <- plot_summs(thym_fit3, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("hypoxia","notch","pi3k","p53","cell_cycle"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Thymoma/checked_models")

#notch
thym <- temp_meta_score[temp_meta_score$project_id == "TCGA-THYM", colnames(temp_meta_score) %in% colnames_wanted]
thym_fit4 <- lm(notch~ . , data =thym)
summary(thym_fit)
#n asian = 3
name <- paste0("Thymoma (stage unknown)\nnotch (N=38)")
t <- plot_summs(thym_fit4, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("hypoxia","metastasis","pi3k","p53","cell_cycle"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Thymoma/checked_models")

#hypoxia 
thym <- temp_meta_score[temp_meta_score$project_id == "TCGA-THYM", colnames(temp_meta_score) %in% colnames_wanted]
thym_fit5 <- lm(hypoxia~ . , data =thym)
summary(thym_fit)
#n asian = 3
name <- paste0("Thymoma (stage unknown)\nhypoxia (N=38)")
t <- plot_summs(thym_fit5, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("notch","metastasis","pi3k","p53","cell_cycle"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Thymoma/checked_models")

#cell cycle
thym <- temp_meta_score[temp_meta_score$project_id == "TCGA-THYM", colnames(temp_meta_score) %in% colnames_wanted]
thym_fit6 <- lm(cell_cycle~ . , data =thym)
summary(thym_fit)
#n asian = 3
name <- paste0("Thymoma (stage unknown)\ncell_cycle (N=38)")
t <- plot_summs(thym_fit6, plot.distributions = T, scale = T, inner_ci_level = 0.9, coefs = c("notch","metastasis","pi3k","p53","hypoxia"))
t <- t + ggtitle(name)+ theme(plot.title = element_text(size = 11))+ theme(plot.title.position = "plot")
t
ggsave(filename = paste0(name, ".png"), path = "/Users/michellemeier/Semesterproject/ARCHS4/results/TCGA/Thymoma/checked_models")

export_summs(thym_fit1, thym_fit2, thym_fit3,thym_fit4,thym_fit5,thym_fit6, scale = T, to.file = "docx", file.name = "THYM_coefficients.docx",
             statistics = summ.lm, model.names = pathways)










