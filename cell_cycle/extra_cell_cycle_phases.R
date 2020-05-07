##extra cell cycle:
##analysis of G2 and M phase cells 
#is there differential expression between G2 and M?
#series: GSE64016, GSE73565
list_of_series <- c("GSE64016","GSE73565")

phases_df <- read.csv(file = "/Users/michellemeier/Semesterproject/ARCHS4/cell_cycle/phases_cell_cycle.csv", header = FALSE, stringsAsFactors = FALSE)
g1_list <- as.vector(phases_df$V1)
s_list <- as.vector(phases_df$V2)
g2_list <- as.vector(phases_df$V3)
m_list <- as.vector(phases_df$V4)

for (entry in g1_list){
  ccp_res_4$condition <- addLevel(ccp_res_4$condition, "G1")
  ccp_res_4[ccp_res_4$sample== entry, 4] = "G1"
}

for (entry in s_list){
  ccp_res_4$condition <- addLevel(ccp_res_4$condition, "S")
  ccp_res_4[ccp_res_4$sample== entry, 4] = "S"
}

for (entry in g2_list){
  ccp_res_4$condition <- addLevel(ccp_res_4$condition, "G2")
  ccp_res_4[ccp_res_4$sample== entry, 4] = "G2"
}

for (entry in m_list){
  ccp_res_4$condition <- addLevel(ccp_res_4$condition, "M")
  ccp_res_4[ccp_res_4$sample== entry, 4] = "M"
}

small_df <- ccp_res_4[ccp_res_4$series %in% list_of_series,]


df <- subset(ccp_res_4, series == "GSE64016")
dim = dim(df)
samples_for_series <- subset(ccp_res_4, ccp_res_4$series == "GSE64016")
samples_for_series_v = as.vector(samples_for_series$sample)
expression_ccp_s = expression_ccp_scale[,colnames(expression_ccp_scale) %in% samples_for_series_v]
#scale
series_tsne <- tSNE(expression_ccp_s,1)
new = as.data.frame(series_tsne$Y)
plot1 <- ggplot(new, aes(x= V1, y= V2,colour= df$condition, label= df$sample)) + ggtitle("GSE64016") + geom_point(size =0.5) + geom_text(size = 2, vjust =0, nudge_y = 5) + scale_color_discrete("condition")
show(plot1)
ggsave(plot = plot1,filename ="Cell_cycle_GSE64016.png", path = "/Users/michellemeier/Semesterproject/ARCHS4/results/cell_cycle/")

df <- subset(ccp_res_4, series == "GSE73565")
dim = dim(df)
samples_for_series <- subset(ccp_res_4, ccp_res_4$series == "GSE73565")
samples_for_series_v = as.vector(samples_for_series$sample)
expression_ccp_s = expression_ccp_scale[,colnames(expression_ccp_scale) %in% samples_for_series_v]
series_tsne <- tSNE(expression_ccp_s,1)
new = as.data.frame(series_tsne$Y)
plot1 <- ggplot(new, aes(x= V1, y= V2, colour= df$condition, label= df$sample)) + ggtitle("GSE73565") + geom_point(size =0.5) + geom_text(size = 2, vjust =0, nudge_y = 5) + scale_color_discrete("condition")
show(plot1)
ggsave(plot = plot1,filename ="Cell_cycle_GSE73565.png", path = "/Users/michellemeier/Semesterproject/ARCHS4/results/cell_cycle/")


#ssGSEA
ssGSEA(gene_set_list_ccp, small_df, expression_ccp_scale, "cell_cycle",1)
