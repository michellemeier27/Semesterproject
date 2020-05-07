## PATHWAY SCORINGFUNCTION FOR ALL SERIES
#loading libraries
library(GSVA)
library(singscore)
library(ggplot2)

#loading other files 
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/FindRankSeries.R")
source("/Users/michellemeier/Semesterproject/ARCHS4/functions/FindTruePathwayRank.R")

#-> testing for hypoxia
#generate random pathways
#get all gene sets and unlist
gene_sets <- read.gmt("/Users/michellemeier/Semesterproject/ARCHS4/raw_data/gene_sets/hallmark_gene_sets.gmt")
all_genes <- unique(unlist(gene_sets))
#length of wanted pathway (in this case hypoxia)
l <- length(true_hypoxia_gs[[1]])

gene_set_random <- list()

for (i in 1:100){
  name <- paste0("random",i)
  random <- sample(all_genes, l)
  gene_set_random[[name]] <- random 
}

#generate true pathway
true_hypoxia_gs <- list(gene_set_list_hypoxia[[1]])
names(true_hypoxia_gs) <- c("hypoxia")

#generate random df
L <- length(gene_set_random)
random_res <- data.frame(pathway = integer(L), delta_gsva = integer(L),delta_ssgsea = integer(L), delta_zscore = integer(L), delta_plage = integer(L),   delta_singscore = integer(L))


#final data frame
unique_series <- unique(hypoxia_res_5$series)
N <- length(unique_series) 
pathway_pos = data.frame(gsva = integer(N), ssGSEA = integer(N), zscore = integer(N), PLAGE = integer(N), singscore = integer(N), ssPATHS = integer(N))


#add to input: randoms for sspaths, res_df for sspaths, pathway string, weights 
# -> maybe makes more sense to just add it after this df is complete
AllSeriesRank <- function(meta_frame,expression_frame, true_gene_set, random_list = gene_set_random, final_df = pathway_pos, random_df = random_res){
  #initialise df to store all positions 
  #get all series
  unique_series <- unique(meta_frame$series)
  #for all series
  i = 0
  for (series_id in unique_series){
    print(series_id)
    i = i+1
    pos_series <- FindRanksSeries(series_id = series_id, meta_frame = meta_frame, expression_frame = expression_frame, true_gene_set =true_gene_set, randoms = random_list, random_df = random_df)
    #pos_sspaths <- ssPATHS (...) -> adjust input of this function
    print(pos_series)
    final_df[[1]][i] = pos_series[1]
    final_df[[2]][i] = pos_series[2]
    final_df[[3]][i] = pos_series[3]
    final_df[[4]][i] = pos_series[4]
    final_df[[5]][i] = pos_series[5]
    #final_df[[6]][i] = pos_sspaths -> final_df -> bigger dimension by 1 col
  
  }
  final_df <- na.omit(final_df)
  return(final_df)
}


  