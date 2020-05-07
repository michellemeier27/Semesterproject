##FUNCTION TO GET EXPRESSION DATA FROM H5 FILE

#loading libraries
library(rhdf5)

#h5 file
destination_file = "/Users/michellemeier/Semesterproject/ARCHS4/raw_data/human_matrix.h5" # this is the location of your h5 file
accession_list = h5read(destination_file, "meta/Sample_geo_accession")
genes_list = h5read(destination_file, "meta/genes")

#define function to get expression values
expression <- function(data_frame, accession = accession_list, genes = genes_list, file = destination_file, cut_off = 10){
  list = c()
  all_sample_geo = data_frame[,2]
  sample_geo = which(accession %in% all_sample_geo)
  expression_results = h5read(file, "data/expression", index=list(1:length(genes), sample_geo))
  h5closeAll()
  rownames(expression_results) = genes
  colnames(expression_results) = accession[sample_geo]
  vectorised = as.vector(data_frame$sample)
  expression_results = expression_results[, vectorised]
  #expression = expression_results[!(apply(expression_results, 1, function(y) any(y < cut_off))),]
  return(expression_results)
}








