##FUNCTION TSNE
#loading library
library(Rtsne)
#tSNE function for individual series 
tSNE <- function(expression_results, perplexity_wanted=5){
  transpose = t(expression_results)
  tsne_results <- Rtsne(transpose, perplexity = perplexity_wanted , check_duplicates = FALSE)
}


