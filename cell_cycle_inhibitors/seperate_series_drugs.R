#seperate series with multiple cell lines into subseries 
for (entry in subseries_drugs_1){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 2
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_1")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

for (entry in subseries_drugs_2){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 3
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_2")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

for (entry in subseries_drugs_3){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 4
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_3")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

for (entry in subseries_drugs_4){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 5
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_4")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

for (entry in subseries_drugs_5){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 6
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_5")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

for (entry in subseries_drugs_6){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 7
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_6")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

for (entry in subseries_drugs_7){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 6
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_7")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

for (entry in subseries_drugs_8){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 7
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_8")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

for (entry in subseries_drugs_9){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 8
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_9")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

for (entry in subseries_drugs_10){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 9
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_10")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

for (entry in subseries_drugs_11){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 10
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_11")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

for (entry in subseries_drugs_12){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 11
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_12")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

for (entry in subseries_drugs_13){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 12
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_13")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

for (entry in subseries_drugs_14){
  drugs_res_1[drugs_res_1$sample_accession == entry , 7] = 13
  series = drugs_res_1[drugs_res_1$sample_accession == entry , 1]
  series_new = paste0(series, "_14")
  drugs_res_1$series_id <- addLevel(drugs_res_1$series_id, series_new)
  drugs_res_1[drugs_res_1$sample_accession == entry , 1] = series_new
}

