##FUNCTION TO NORMALISE AND SCALE
#loading libraries
library(scone)
library(edgeR)

#defining function to scale, centre and normalise
#later note: only using this to scale, TPM still not possible due to ambiguous annotation in ARCHS4
normalise_scale <- function(expression_matrix,CLR = FALSE, log = TRUE, scale = TRUE, UQ = TRUE){
  #CLR
  if (CLR == TRUE){
    expression_matrix = CLR_FN(expression_matrix)
  }
  #UQ
  if (UQ == TRUE){
    expression_matrix = UQ_FN(expression_matrix)
  }
  #TPM
  #expression_matrix = calculateTPM(expression_matrix, gene_lengths$length)
  #scale
  if (scale == TRUE){
    expression_matrix = scale(expression_matrix, scale = TRUE)
  }
  if (log == TRUE){
    expression_matrix = log(expression_matrix+1)
  }
  return(expression_matrix)
}

#defining function for filtering out lowly expressed genes (to about 12'000 genes)
cutoff <- function(expression_frame, meta_frame, cutoff_min=10, cutoff_min_total = 15){
  dgelist <- DGEList(counts = expression_frame, samples = meta_frame, group = meta_frame$series)
  keep <- filterByExpr(dgelist, min.count = cutoff_min, min.total.count = cutoff_min_total)
  new_dgelist <- dgelist[keep, , keep.lib.sizes=FALSE]
  return(new_dgelist$counts)
}








