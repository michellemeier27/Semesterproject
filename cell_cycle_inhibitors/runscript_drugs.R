## finding CDK4/6 inhibitor samples
search_str_inhibitors <- c("palbociclib|abemaciclib|ribociclib|cdk4|cdk6")
blacklist_drugs = c("")
drugs_res_1 = get_res_table( search_str_inhibitors, blacklist_drugs, samp_desc, extract_protocol,
                           sample_title, sample_char, sample_source, experiment_str="drugs")

#read in sample annotation
source("/Users/michellemeier/Semesterproject/ARCHS4/cell_cycle_inhibitors/manual_drugs.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/cell_cycle_inhibitors/seperate_series_drugs.R")

#manual samples 
manual <- read.csv(file = "/Users/michellemeier/Semesterproject/ARCHS4/cell_cycle_inhibitors/manual_samples_cell_cycle_inhibitors.csv", skip = 1, na.strings = "")
df_1 <- manual[,1:2]
df_2 <- data.frame(manual[,4])
df_2 <- na.omit(df_2)
df_2$condition = "placebo"
df_3 <- data.frame(manual[,5])
df_3 <- na.omit(df_2)
df_3$condition = "Palbociclib"
df_4 <- rbind(df_2,df_3)
colnames(df_4) <- c("samples", "condition")
manual_df <- merge(df_1,df_4, by.x = "samples", by.y = "samples")
done_manual_df <- data.frame(series = manual_df$series,
                             sample = manual_df$samples,
                             original_series = manual_df$series,
                             condition = manual_df$condition)
#make data frame
placebo <- sample_annotation$placebo
placebo <- na.omit(placebo)
BSJ.03.123 <- sample_annotation$BSJ.03.123
BSJ.03.123 <- na.omit(BSJ.03.123)
YKL.06.102 <- sample_annotation$YKL.06.102
YKL.06.102 <- na.omit(YKL.06.102)
palbo <- sample_annotation$Palbociclib
palbo <- na.omit(palbo)
ribo <- sample_annotation$ribociclib
ribo <- na.omit(ribo)
abema <- sample_annotation$abemaciclib
abema <- na.omit(abema)
resistant_p <- sample_annotation$Palbociclib.resistant
resistant_p <- na.omit(resistant_p)
all_samples_drugs <- c(placebo,BSJ.03.123,YKL.06.102,palbo,ribo,abema,resistant_p)
l <- length(all_samples_drugs)
df <- drugs_res_1[drugs_res_1$sample_accession %in% all_samples_drugs,]
summary_cell_cycle_inhibitors <- df[, -c(2,4,5,6,7,8)]
status <- integer(l)
summary_cell_cycle_inhibitors <- add_column(summary_cell_cycle_inhibitors, status)
summary_cell_cycle_inhibitors[summary_cell_cycle_inhibitors$sample_accession %in% placebo,4] = "placebo"
summary_cell_cycle_inhibitors[summary_cell_cycle_inhibitors$sample_accession %in% BSJ.03.123,4] = "BSJ.03.123"
summary_cell_cycle_inhibitors[summary_cell_cycle_inhibitors$sample_accession %in% YKL.06.102,4] = "YKL.06.102"
summary_cell_cycle_inhibitors[summary_cell_cycle_inhibitors$sample_accession %in% palbo,4] = "Palbociclib"
summary_cell_cycle_inhibitors[summary_cell_cycle_inhibitors$sample_accession %in% ribo,4] = "ribociclib"
summary_cell_cycle_inhibitors[summary_cell_cycle_inhibitors$sample_accession %in% abema,4] = "abemaciclib"
summary_cell_cycle_inhibitors[summary_cell_cycle_inhibitors$sample_accession %in% resistant_p,4] = "Palbociclib.resistant"
colnames(summary_cell_cycle_inhibitors) <- c("series", "sample", "original_series", "condition")
#update with manual data
summary_cell_cycle_inhibitors = rbind(summary_cell_cycle_inhibitors, done_manual_df)

#getting expression data
expression_drugs_raw <- expression(summary_cell_cycle_inhibitors)
#log transforms 
expression_drugs_prep = log(expression_drugs_raw +1)

#cutoff for genes to only have significantly expressed genes
expression_drugs_cutoff = cutoff(expression_drugs_raw, summary_cell_cycle_inhibitors)

#lib size normalisation 
expression_drugs = UQ_FN(expression_drugs_cutoff)

#scaling data before tSNE overview:
expression_drugs_scale = normalise_scale(expression_drugs, CLR = FALSE, log = FALSE)

##tSNE
#overview: data log, clr, scale and cutoff
#calculations:
tsne_drugs = tSNE(expression_drugs_scale)
#plotting series (original) and series (subseries) overview
plotting_overview(tsne_drugs, summary_cell_cycle_inhibitors, "cell_cycle_inhibitors")

#tSNE for each series: data log, clr, scale and cutoff
plotting_series(summary_cell_cycle_inhibitors, expression_drugs_scale, "cell_cycle_inhibitors", "cell_cycle_inhibitors/")

##PCA
#for all series individually: data log, clr, scale and cutoff
pca_series(summary_cell_cycle_inhibitors, expression_drugs_scale, "cell_cycle_inhibitors", "cell_cycle_inhibitors/")

# write tsv file. 
write_file(summary_cell_cycle_inhibitors, "results/cell_cycle_inhibitors/results.tsv")







