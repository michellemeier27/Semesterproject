#file to devide series into subseries for PI3K 
for (entry in subseries_1_pi3k){
  pi3k_res_2[pi3k_res_2$sample_accession == entry , 8] = 1
  series = pi3k_res_2[pi3k_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_1")
  pi3k_res_2$series_id <- addLevel(pi3k_res_2$series_id, series_new)
  pi3k_res_2[pi3k_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_2_pi3k){
  pi3k_res_2[pi3k_res_2$sample_accession == entry , 8] = 2
  series = pi3k_res_2[pi3k_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_2")
  pi3k_res_2$series_id <- addLevel(pi3k_res_2$series_id, series_new)
  pi3k_res_2[pi3k_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_3_pi3k){
  pi3k_res_2[pi3k_res_2$sample_accession == entry , 8] = 3
  series = pi3k_res_2[pi3k_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_3")
  pi3k_res_2$series_id <- addLevel(pi3k_res_2$series_id, series_new)
  pi3k_res_2[pi3k_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_4_pi3k){
  pi3k_res_2[pi3k_res_2$sample_accession == entry , 8] = 4
  series = pi3k_res_2[pi3k_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_4")
  pi3k_res_2$series_id <- addLevel(pi3k_res_2$series_id, series_new)
  pi3k_res_2[pi3k_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_5_pi3k){
  pi3k_res_2[pi3k_res_2$sample_accession == entry , 8] = 5
  series = pi3k_res_2[pi3k_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_5")
  pi3k_res_2$series_id <- addLevel(pi3k_res_2$series_id, series_new)
  pi3k_res_2[pi3k_res_2$sample_accession == entry , 1] = series_new
}

for (entry in subseries_6_pi3k){
  pi3k_res_2[pi3k_res_2$sample_accession == entry , 8] = 6
  series = pi3k_res_2[pi3k_res_2$sample_accession == entry , 1]
  series_new = paste0(series, "_6")
  pi3k_res_2$series_id <- addLevel(pi3k_res_2$series_id, series_new)
  pi3k_res_2[pi3k_res_2$sample_accession == entry , 1] = series_new
}