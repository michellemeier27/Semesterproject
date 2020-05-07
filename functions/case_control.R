##FUNCTION TO CORRECT CASE CONTROL ANNOTATION

#defining function to correct for case control annotation 
case_control <- function(case_list, control_list, data_frame){
  for (entry in control_list){
    data_frame[data_frame$sample_accession== entry, 4] = "control"
  }
  
  for (entry in case_list){
    data_frame[data_frame$sample_accession== entry, 4] = "case"
  }
  return(data_frame)
}
