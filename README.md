This is the github repository for my semesterproject. It contains the R scripts used to generate all the data presented in the final report, intermediate results for the linear regression and the cell cycle inhibitor analysis and the final benchmarking data set for all pathways. Please note that in general, all scripts contain local paths that must be adjusted before using them. I also added my finished report, in case someone needs something to read during this pandemic.

# Benchmarking data set
The folder benchmarking_data_sets contains 6 subfolders, one for each pathway. In each folder, you can find the runfile which combines all functions, a list of samples removed from further analysis, a list of defined subseries, a list of case/control annotations, a R script combining technical replicates and a R script called manual which was used for miscellaneous manual curation of the data. Running the runfile should produce the data set for that pathway as a tsv file. It also produces all t-SNE, PCA and ssGSEA plots that were used to assess the quality of the series/sample. Please note that you will need to download the human expression matrix (gene level) from the ARCHS4 website and, if you want to plot ssGSEA biplots, you will also need to download the corresponding gene sets from MSigDB.
Additionally, each folder contains a hallmark script which will produce boxplots comparing the performance of different pathway scoring methods (GSVA, ssGSEA, zscore, PLAGE, singscore) in identifying the correct hallmark gene set out of all 50 gene sets. Please note that you will need to download all hallmark gene sets from MSigDB. 
Lastly, the R script called pathway score compares the performance of the pathway scoring methods mentioned above in identifying the correct gene set out of 100 randomly generated gene sets. These results were plotted using a new script called plotting_pathway_scores.R in the benchmarking_data_sets folder. 

# Cell cycle inhibitors
Similarly to how the benchmarking data sets were generated, the runscript in the cell_cycle_inhibitors folder produces the samples that were used for this study. You can also find an excel sheet with an overview of these samples in the folder. The two scripts manual_drugs.R and seperate_series_drugs.R were used to split series into subseries after manual curation. Also, there is a list of manually added samples/series in a csv file. Further the script called cell_cycle_inhibitor_pathway_score.R was used to calculate the pathway score for cell cycle in these samples and then visualised using the visualizing_plots script.

# TCGA Linear regression analysis
The folder called TCGA contains one R script for each cancer that was studied. In each script, linear regression analysis is performed separately for all pathways and you should be able to get the results by simply running this script. The prepping_data script merges the metadata (meta_data_tcga) and the calculated scores (singscore_TCGA) for further analysis. Thus, this script must be run first.

# Functions
In order to (at least attempt to) keep my project tidy, each function usually has its own designated script. The following table gives you a rough overview of what each script/function does:

| script | What it does | 
| --- | --- |
| AllSeriesRanks | combines the ranks of the true gene set out of 100 random gene sets for all series and scoring methods (using FindRankSeries) |
|condition | semi-automatic case/control annotation based on the occurence of specific words |
| expression_data | gets expression data for a list of sample accessions|
| FindRankSeries |  finds the ranks of the true gene set out of 100 random gene sets for one series for all scoring methods (using FindTruePathwayRank |
| FindTruePathwayRank | finds the rank of the true gene set for all scoring methods | 
| GSVA_function | finds the rank of the true hallmark gene set out of all gene sets using the scoring methods available in the GSVA package | 
| ignoring_finalising | misc. file with multiple functions: 1. case/control annotation in the benchmarking data set, needs list of samples to correct 2. makes a nice final data frame 3. writes tsv file with results to path wanted |
| lib_norm_check_p53 | checks if library normalisation worked by plotting expression data in boxplot for a series |
| merging_tech_replicates | merges technical replicates, needs list of replicates | normalise_scale | two functions: 1. library normalisation (upper quantile UQ or centred-log-ratio CLR), log transformation and gene length correction (did not work on this data set) 2. defines expression cutoff | 
| PCA | plots PCA for series and/or samples |
| plotting_series | produces t-SNE plots for each series (using tSNE) | 
| plotting_tsne_series_overview | produces t-SNE plot for all series | 
| search_function | searches the h5 file downloaded from ARCHS4, needs search terms | 
| singscore_function |  finds the rank of the true hallmark gene set out of all gene sets using singscore | 
| ssGSEA_function | produces ssGSEA biplots for each series, needs gene sets | 
| tSNE | transposes data and gives tSNE results |

# Contact
If you have any further questions, feel free to contact me via email: mimeier@biol.ethz.ch
