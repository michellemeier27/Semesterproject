##visualising pathway results
##for palbiciclib
total_p = c("GSE116187","GSE133568","GSE110397_6","GSE110397_5","GSE110397_4","GSE110397_3","GSE110397_2",
            "GSE110397_1","GSE130903","GSE120920_2","GSE115976","GSE74620") 
cell_cycle_dys_p <- c("GSE116187","GSE133568","GSE110397_6","GSE110397_5","GSE110397_2","GSE120920_2","GSE115976",
                      "GSE74620","GSE130903")
notch_dys_p <- c("GSE110397_5")
hypoxia_dys_p <- c("GSE133568","GSE110397_6","GSE120920_2","GSE74620","GSE130903")
metastasis_dys_p <- c("GSE133568","GSE110397_6","GSE110397_2","GSE74620")#"GSE110397_4","GSE110397_1"
p53_dys_p <- c("GSE133568","GSE120920_2","GSE74620")#"GSE110397_4"
pi3k_dys_p <- c("GSE110397_5","GSE110397_2")

##for resistance
total_r = c("GSE130437","GSE128056") 
l_r <- length(total_r)
cell_cycle_dys_r <- c("GSE130437")
lc_r <- length(cell_cycle_dys_r)
notch_dys_r <- c("GSE130437","GSE128056")
hypoxia_dys_r <- c()
metastasis_dys_r <- c("GSE128056")
p53_dys_r <- c("GSE130437")

l_p <- length(total_p)
lc_p <- length(cell_cycle_dys_p)
l_r <- length(total_r)
lc_r <- length(cell_cycle_dys_r)

##overview plot for cell cycle dysregulation 
cell_cycle_comp_df <- data.frame(not_significant = c(l_p- lc_p,l_r- lc_r),
                                 significant = c(lc_p,lc_r),
                                 type = c("sensitive","resistant"))
t <-stack(cell_cycle_comp_df)
t$type = c("sensitive\n(N=12)","resistance\n(N=2)")
colnames(t) <- c("score","sig","type")

g <- ggbarplot(data = t, x = "type", y = "score", fill = "sig", position = position_dodge(0.8), xlab = "", ylab = "#series")
g <- g + ggtitle("differential cell cycle pathway regulation\n(N=14)") +theme(legend.position = "right") +theme(legend.title = element_blank())+rotate_x_text(angle =45)
g <- g+ theme(axis.text = element_text(size = 11)) + theme(plot.title = element_text(size = 13))
g <- g+ annotate("text", x = 0.8, y = 2.5, label ="25%")
g <- g+ annotate("text", x = 1.2, y = 8.5, label ="75%")
g <- g+ annotate("text", x = 1.8, y = 0.5, label ="50%")
g <- g+ annotate("text", x = 2.2, y = 0.5, label ="50%")
g
ggsave(filename = "cell_cycle_dysregulation_overview.png", path = my_directory)


##out of all sensitive strains pathway distributions
#calculating percentages
combined_pathway_sensitive <- data.frame(percentage = c(6,6,4,1,2),
                                         pathway = c("hypoxia","metastasis","p53","notch","pi3k"))

bp_1 <- ggbarplot(data = combined_pathway_sensitive, x = "pathway", y = "percentage", fill = "navyblue", xlab = "", ylab = "sample count")
bp_1 <- bp_1 + ggtitle("Pathway dysregulation in reported\ndrug sensitive samples (N = 12)")+theme(title = element_text(size =11))+theme(axis.text = element_text(size = 9))
bp_1
ggsave(filename = "pathway_dysregulation_overview_sensitive.png", path = my_directory)


##out of all resistant strains pathway distributions
l_r <- length(total_r)
hypoxia_np <- length(hypoxia_dys_r)
metastasis_np <- length(metastasis_dys_r)
p53_np <- length(p53_dys_r)
notch_np <- length(notch_dys_r)

combined_pathway_resistance <- data.frame(percentage = c(hypoxia_np,metastasis_np,p53_np,notch_np,0),
                                         pathway = c("hypoxia","metastasis","p53","notch","pi3k"))

bp_1 <- ggbarplot(data = combined_pathway_resistance, x = "pathway", y = "percentage", fill = "#5ab4ac", xlab = "", ylab = "sample count")
bp_1 <- bp_1 + ggtitle("Pathway dysregulation in reported\ndrug resistant samples (N = 2)")+theme(title = element_text(size =11))+theme(axis.text = element_text(size = 9))
bp_1
ggsave(filename = "pathway_dysregulation_overview_resistant.png", path = my_directory)


#new plot for better overview 
comparison_pathway <- data.frame(percentage = c(50,50, 30 ,8.3,16.7, 0,50,50,100,0),
                                 condition = c("sensitive","sensitive","sensitive","sensitive","sensitive","resistant","resistant","resistant","resistant","resistant"),
                                 pathway = c("hypoxia","metastasis","p53","notch","pi3k","hypoxia","metastasis","p53","notch","pi3k"))

bp <- ggbarplot(data =comparison_pathway, x= "pathway", y= "percentage", fill = "condition",position = position_dodge(0.8), xlab = "", ylim = c(1,115))
bp <- bp + labs(fill = "drug sensitivity") +theme(legend.position = "right") +theme(legend.title = element_text(size = 9))+ annotate ("text", x= 1, y = 105, label ="ns",size = 5)
bp <- bp + theme(axis.text = element_text(size = 10)) +ggtitle("Pathway dysregulation in all samples\n(N=14)") +theme(title = element_text(size =11))
bp <- bp + theme(axis.title = element_text(size = 10)) + rotate_x_text(angle = 45) + scale_fill_discrete(labels = c("resistant (N=2)", "sensitive (N=12)"))
bp <- bp + annotate ("text", x= 2, y = 105, label ="ns",size = 5)+ annotate ("text", x= 3, y = 105, label ="*",size = 5)+ annotate ("text", x= 4, y = 105, label ="ns",size = 5)
bp <- bp + annotate ("text", x= 1.2, y = 115, label ="Fisher's exact Test:", size = 5) +annotate ("text", x= 5, y = 105, label ="ns",size = 5)
bp <- bp +scale_y_continuous(breaks=c(20,40, 60, 80,100))
bp

ggsave(filename = "pathway_dysregulation_overview_all_samples.png", path = my_directory)






