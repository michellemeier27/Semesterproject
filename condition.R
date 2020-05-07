##FUNCTION TO ANNOTATE CASE AND CONTROLS

#loading libraries
library(rhdf5)
library(tidyverse)

#getting necessary stuff ready
destination_file = "/Users/michellemeier/Semesterproject/ARCHS4/raw_data/human_matrix.h5" # this is the location of your h5 file
characteristics = h5read(destination_file, "meta/Sample_characteristics_ch1")
characteristics<- enc2utf8(characteristics)
characteristics = tolower(characteristics)
title = h5read(destination_file, "meta/Sample_title")
title<- enc2utf8(title)
title = tolower(title)
accession = h5read(destination_file, "meta/Sample_geo_accession")

#defining function that adds case/control vector
condition <- function(data_frame, condition_string){
  sample_ids = data_frame[,3]
  titles = data_frame[,2]
  char = data_frame [,5]
  control_title = grep(condition_string, titles,useBytes = TRUE)
  control_char = grep(condition_string, char, useBytes = TRUE)
  all_controls = union(control_title, control_char)
  
  condition_vector = rep("case", length(sample_ids))
  for (control in all_controls){
    condition_vector[control] = 'control'
  }
  
  result = add_column(data_frame, condition_vector, .after = 3)
  
  return(result)
}





  
