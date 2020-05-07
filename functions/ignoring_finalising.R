#FUNCTIONS TO FINALISE META DATA FRAME

#defining function to correct wrong control case annotation
correcting <- function(data_frame, controls, cases){
  for (entry in controls){
    data_frame[data_frame$sample_accession== entry, 4] = "control"
  }
  for (entry in cases){
    data_frame[data_frame$sample_accession== entry, 4] = "case"
  }
}


#defining function to integrate everything neatly in a dataframe
finalising_df <- function(data_frame){
  result = data.frame(series = data_frame$series_id,
                      sample = data_frame$sample_accession,
                      pathway = data_frame$experiment,
                      condition = data_frame$condition_vector,
                      case = data_frame$one_control,
                      direction = data_frame$direction_ES,
                      original_series = data_frame$original_series_id)
  return(result)
}



#defining function to write file to export all samples with meta data 
write_file <- function(data_frame, file_name){
  write.table(data_frame, file = file_name, sep ="\t", quote = FALSE)
  print(paste0("Expression file was created at ", getwd(), "/", file_name))
}

#defining function to add factor level
addLevel <- function(x, newlevel=NULL) {
  if(is.factor(x)) {
    if (is.na(match(newlevel, levels(x))))
      return(factor(x, levels=c(levels(x), newlevel)))
  }
  return(x)
}


