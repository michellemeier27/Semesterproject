##BOXPLOTS WITH STATISTICS
#libraries
library(ggbiplot)
library(ggpubr)
#all with reference to singscore

#####
#hypoxia
#####
#transform data for plotting
stack.hypoxia <- stack(data_res_hypoxia)
colnames(stack.hypoxia) <- c("rank", "method")
#plotting
bp_hypoxia <- ggboxplot(stack.hypoxia, x = "method", y = "rank", fill = "method", ylim = c(0,105)) +theme(axis.text=element_text(size=9)) + geom_hline(yintercept= mean(stack.hypoxia$rank), linetype =2) + theme(axis.title.x=element_blank())+ rotate_x_text(angle = 45)+
  ggtitle("Pathway score hypoxia (cases)", subtitle = "100 random pathways")+ theme(legend.position= "right")+theme(legend.text = element_text(size = 8)) +theme(legend.title = element_text(size = 8))+theme(plot.title=element_text(size=13))+ theme(plot.subtitle=element_text(size=11))+
          stat_compare_means( method = "t.test", ref.group = "singscore", hide.ns = TRUE,label = "p.signif", size = 5, label.y = 103) 
ggsave(filename = "hypoxia_stats.png",path = "/Users/michellemeier/Semesterproject/ARCHS4/results/pathway_scoring/case_direction", plot = bp_hypoxia)

#####
#p53
#####
#transform data for plotting
stack.p53 <- stack(data_res_p53)
colnames(stack.p53) <- c("rank", "method")
#plotting
bp_p53 <- ggboxplot(stack.p53, x = "method", y = "rank", fill = "method", ylim = c(0,105)) +theme(axis.text=element_text(size=9)) + geom_hline(yintercept= mean(stack.p53$rank), linetype =2) + theme(axis.title.x=element_blank())+ rotate_x_text(angle = 45)+
  ggtitle("Pathway score p53 (controls)", subtitle = "100 random pathways")+ theme(legend.position= "right")+theme(legend.text = element_text(size = 8)) +theme(legend.title = element_text(size = 8))+theme(plot.title=element_text(size=13))+ theme(plot.subtitle=element_text(size=11))+
  stat_compare_means( method = "t.test", ref.group = "singscore", hide.ns = TRUE,label = "p.signif", size = 5, label.y = 103) 
show(bp_p53)
ggsave(filename = "p53_stats.png",path = "/Users/michellemeier/Semesterproject/ARCHS4/results/pathway_scoring/case_direction", plot = bp_p53)


#####
#pi3k
#####
#transform data for plotting
stack.pi3k <- stack(data_res_pi3k)
colnames(stack.pi3k) <- c("rank", "method")
#plotting
bp_pi3k <- ggboxplot(stack.pi3k, x = "method", y = "rank", fill = "method", ylim = c(0,105)) +theme(axis.text=element_text(size=9)) + geom_hline(yintercept= mean(stack.pi3k$rank), linetype =2) + theme(axis.title.x=element_blank())+ rotate_x_text(angle = 45)+
  ggtitle("Pathway score pi3k (cases)", subtitle = "100 random pathways")+ theme(legend.position= "right")+theme(legend.text = element_text(size = 8)) +theme(legend.title = element_text(size = 8))+theme(plot.title=element_text(size=13))+ theme(plot.subtitle=element_text(size=11))+
  stat_compare_means( method = "t.test", ref.group = "singscore", hide.ns = TRUE,label = "p.signif", size = 5, label.y = 103) 
show(bp_pi3k)
ggsave(filename = "pi3k_stats.png",path = "/Users/michellemeier/Semesterproject/ARCHS4/results/pathway_scoring/case_direction", plot = bp_pi3k)

#####
#notch
#####
#transform data for plotting
stack.notch <- stack(data_res_notch)
colnames(stack.notch) <- c("rank", "method")
#plotting
bp_notch <- ggboxplot(stack.notch, x = "method", y = "rank", fill = "method", ylim = c(0,105)) +theme(axis.text=element_text(size=9)) + geom_hline(yintercept= mean(stack.notch$rank), linetype =2) + theme(axis.title.x=element_blank())+ rotate_x_text(angle = 45)+
  ggtitle("Pathway score notch (cases)", subtitle = "100 random pathways")+ theme(legend.position= "right")+theme(legend.text = element_text(size = 8)) +theme(legend.title = element_text(size = 8))+theme(plot.title=element_text(size=13))+ theme(plot.subtitle=element_text(size=11))+
  stat_compare_means( method = "t.test", ref.group = "singscore", hide.ns = TRUE,label = "p.signif", size = 5, label.y = 103) 
show(bp_notch)
ggsave(filename = "notch_stats.png",path = "/Users/michellemeier/Semesterproject/ARCHS4/results/pathway_scoring/case_direction", plot = bp_notch)

#####
#cell cycle
#####
#transform data for plotting
stack.ccp <- stack(data_res_ccp)
colnames(stack.ccp) <- c("rank", "method")
#plotting
bp_ccp <- ggboxplot(stack.ccp, x = "method", y = "rank", fill = "method", ylim = c(0,105)) +theme(axis.text=element_text(size=9)) + geom_hline(yintercept= mean(stack.ccp$rank), linetype =2) + theme(axis.title.x=element_blank())+ rotate_x_text(angle = 45)+
  ggtitle("Pathway score cell cycle (cases)", subtitle = "100 random pathways")+ theme(legend.position= "right")+theme(legend.text = element_text(size = 8)) +theme(legend.title = element_text(size = 8))+theme(plot.title=element_text(size=13))+ theme(plot.subtitle=element_text(size=11))+
  stat_compare_means( method = "t.test", ref.group = "singscore", hide.ns = TRUE,label = "p.signif", size = 5, label.y = 103) 
show(bp_ccp)
ggsave(filename = "cell_cycle_stats.png",path = "/Users/michellemeier/Semesterproject/ARCHS4/results/pathway_scoring/case_direction", plot = bp_ccp)

#####
#metastasis
#####
#transform data for plotting
stack.met <- stack(data_res_met)
colnames(stack.met) <- c("rank", "method")
#plotting
bp_met <- ggboxplot(stack.met, x = "method", y = "rank", fill = "method", ylim = c(0,105)) +theme(axis.text=element_text(size=9)) + geom_hline(yintercept= mean(stack.met$rank), linetype =2) + theme(axis.title.x=element_blank())+ rotate_x_text(angle = 45)+
  ggtitle("Pathway score metastasis (cases)", subtitle = "100 random pathways")+ theme(legend.position= "right")+theme(legend.text = element_text(size = 8)) +theme(legend.title = element_text(size = 8))+theme(plot.title=element_text(size=13))+ theme(plot.subtitle=element_text(size=11))+
  stat_compare_means( method = "t.test", ref.group = "singscore", hide.ns = TRUE,label = "p.signif", size = 5, label.y = 103) 
show(bp_met)
ggsave(filename = "metastasis_stats.png",path = "/Users/michellemeier/Semesterproject/ARCHS4/results/pathway_scoring/case_direction", plot = bp_met)



