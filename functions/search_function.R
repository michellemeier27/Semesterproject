##QUERY FUNCTION
#load all libraries
library(rhdf5)

#adapted from Natalie's script
get_res_table <- function(search_str,blacklist,  samp_desc , extract_protocol, sample_title , sample_char, sample_source , experiment_str){
  # grep through all the meta data to find related stuff
  desc = grep(search_str, samp_desc, useBytes = TRUE)
  proto = grep(search_str, extract_protocol,useBytes = TRUE)
  title = grep(search_str, sample_title,useBytes = TRUE)
  samp_char = grep(search_str, sample_char,useBytes = TRUE)
  samp_source = grep(search_str, sample_source, useBytes = TRUE)
  all_ids = union(desc, proto) #why not all at once?
  all_ids = union(all_ids, title)
  all_ids = union(all_ids, samp_char)
  all_ids = union(all_ids, samp_source)
  curr_series_id = unique(series_id[all_ids]) 
  a = strsplit(curr_series_id, "Xx-xX") 
  a = lapply(a, function(x) x[1])
  curr_series_id = unlist(a)
  curr_series_id = curr_series_id[!curr_series_id %in% blacklist]  
  #now make a meta_data table that we will annotate by hand later
  a = NA
  for(curr_series in curr_series_id){
    sample_idx = which(series_id %in% curr_series)
    curr_res = data.frame(series_id = curr_series,
                          sample_title=sample_title[sample_idx],
                          sample_accession=sample_accession[sample_idx],
                          experiment=experiment_str,
                          charac = sample_char[sample_idx],
                          one_control = "TRUE",
                          part = 1,
                          direction_ES = "case",
                          original_series_id = curr_series)
    a = rbind(a, curr_res)
  }
  res = na.omit(a)
  return(res)
  
}
# get all the info we need from the h5 file
destination_file = "/Users/michellemeier/Semesterproject/ARCHS4/raw_data/human_matrix.h5" # this is the location of your h5 file
samp_desc = h5read(destination_file, "meta/Sample_data_processing")
samp_desc <- enc2utf8(samp_desc)
samp_desc = tolower(samp_desc)
extract_protocol = h5read(destination_file, "meta/Sample_extract_protocol_ch1")
extract_protocol <- enc2utf8(extract_protocol)
extract_protocol = tolower(extract_protocol)
series_id = h5read(destination_file, "meta/Sample_series_id") #some of the series are connected with Xx-xX
a = strsplit(series_id, "Xx-xX") #only seperates string, still in same number
a = lapply(a, function(x) x[1]) #only first element for all in a 
series_id = unlist(a) #converts to vector
sample_title = h5read(destination_file, "meta/Sample_title") 
sample_title <- enc2utf8(sample_title)
sample_title = tolower(sample_title)
sample_accession = h5read(destination_file, "meta/Sample_geo_accession")
sample_char = h5read(destination_file, "meta/Sample_characteristics_ch1") 
sample_char <- enc2utf8(sample_char)
sample_char = tolower(sample_char)
sample_source = h5read(destination_file, "meta/Sample_source_name_ch1")
sample_source <- enc2utf8(sample_source)
sample_source = tolower(sample_source)



