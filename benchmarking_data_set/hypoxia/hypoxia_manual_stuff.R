#updating final table: run whole script to update table and create tsv with results
library(tidyverse)
#indexes of samples to ignore -> careful, indices not a good idea only run this when you know they are true
wrong = c(4: 9,  12:17,  20: 25,  28,  29,  32,  33,  34,  35,  37, 
 42,  44,  47,  50,  51,  53,  54,  55,  57,  58,  59 , 63,  65 , 66,  68,  69,  71, 87, 90, 94, 95,  96, 97,105, 106,109:111, 113, 122, 125, 127, 129, 130, 131,145:153,240:245,
 264:272,275, 276,279, 280,283, 284, 287, 288, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344,
 345, 349, 350, 351, 352, 353 ,354, 355, 356, 364, 365, 366, 367, 368, 369, 371, 379:381,388:390,391, 392, 393, 394, 395, 396,
 400, 401, 402, 420, 421, 422, 426, 427, 428, 441, 442, 451, 452, 454, 455, 456, 458, 459, 460, 468, 469, 470, 471, 472,
 473, 474, 475, 476, 498, 499, 501, 502, 503, 504, 507, 508,511, 512, 521, 522, 523, 530, 531, 532, 535, 536, 537, 538,
 539, 540, 549, 550,553, 554, 557,558,561, 562, 565,566, 569,570, 573,574, 577,578, 581,583, 585,586, 589,590, 641,642, 644,636, 638, 639)
#correct the indices for my stupid mistake
#wrong_corrected = vapply(wrong, function(x) x-1, numeric(1))
#get geo accessions of indices
#wrong_samples = c()
#for (i in wrong_corrected){
        #sam = c(hypoxia_res_2[i,3])
        #wrong_samples = c(wrong_samples, sam)
#}

add_to_df = c("GSM1295106","GSM1295107", "GSM1295108","GSM2876104","GSM2876106","GSM2876108","GSM2876110","GSM2876112","GSM2876114","GSM2501257","GSM2501258","GSM2501259","GSM2501260", "GSM2501261",
              "GSM2501261","GSM2501263","GSM2501264","GSM2501265","GSM2501266","GSM2791588", "GSM2791589", "GSM3273231", "GSM1654546", "GSM1654544","GSM1654540","GSM1654542","GSM1862653","GSM1862663")


#write csv file to access later (do not run this again!!)
#write_csv(wrong_samples, "/Users/michellemeier/Semesterproject/ARCHS4/wrong_samples_hypoxia.csv")

#control and case list from manual annotation 
control_list = c("GSM3094351","GSM3094352","GSM3094353","GSM2501230","GSM2501231", "GSM2501232","GSM2048413", "GSM2048414", "GSM2048415", "GSM1704487","GSM2466352")
case_list = c("GSM3711368", "GSM3711367","GSM2048422", "GSM2048423", "GSM2048424", "GSM4145788", "GSM4145789","GSM4145790", "GSM1463076", "GSM2390684","GSM2390694","GSM2390698","GSM2390702","GSM2033147",
              "GSM2033148","GSM2033149","GSM3058375","GSM3058376","GSM3058387","GSM3058388","GSM3058399","GSM3058400","GSM3058411","GSM3058412")

##list for split series
#series = c("GSE72437")
subseries_1 = c("GSM1862651","GSM1862652", "GSM1862653", "GSM1862654","GSM1862655", "GSM1862656","GSM1862657","GSM1862658", "GSM1862659", "GSM1862660")
subseries_2a = c("GSM1862661","GSM1862662", "GSM1862663", "GSM1862664","GSM1862665", "GSM1862666","GSM1862667","GSM1862668", "GSM1862669", "GSM1862670")

#series = c("GSE77307")
subseries_1 = c("GSM2048413","GSM2048414","GSM2048415","GSM2048416", "GSM2048417", "GSM2048418")
subseries_2b = c("GSM2048419","GSM2048420", "GSM2048421", "GSM2048422","GSM2048423", "GSM2048424")

#series = c("GSE95280")
subseries_1 = c("GSM2501230", "GSM2501231", "GSM2501232", "GSM2501233", "GSM2501234","GSM2501235", "GSM2501236", "GSM2501237", "GSM2501238", "GSM2501239","GSM2501240","GSM2501241")
subseries_2c = c("GSM2501242","GSM2501243","GSM2501244","GSM2501245", "GSM2501246","GSM2501247")
subseries_3c = c( "GSM2501248", "GSM2501249","GSM2501250","GSM2501251","GSM2501252","GSM2501253","GSM2501254","GSM250155","GSM2501256")
                 

#series = c("GSE106305")
subseries_1 = c("GSM3145511", "GSM3145512","GSM3145513","GSM3145514")
subseries_2d = c("GSM3145517", "GSM3145518", "GSM3145521", "GSM3145522")

#series = ("GSE120886")
subseries_1 = c("GSM3417844","GSM3417845","GSM3417846","GSM3417847","GSM3417848","GSM3417849")
subseries_2e = c("GSM3417850","GSM3417851","GSM3417852","GSM3417853","GSM3417854","GSM3417855")
subseries_3 = c("GSM3417856","GSM3417857","GSM3417858", "GSM3417859","GSM3417860","GSM3417861")

#series = c("GSE132624")
subseries_1 = c("GSM3882337","GSM3882338","GSM3882339","GSM3882340","GSM3882341","GSM3882342","GSM3882343","GSM3882344","GSM3882345","GSM3882346","GSM3882347","GSM3882348")
subseries_2f = c("GSM3882349","GSM3882350","GSM3882351","GSM3882352","GSM3882353","GSM3882354","GSM3882355","GSM3882356","GSM3882357","GSM3882358","GSM3882359","GSM3882360")

subseries_3b = c("GSM3882361","GSM3882362","GSM3882363","GSM3882364","GSM3882365",
                "GSM3882366","GSM3882367","GSM3882368","GSM3882369","GSM3882370","GSM3882371","GSM3882372")
subseries_4 = c("GSM3417856", "GSM3417857","GSM3417858","GSM3417859","GSM3417860","GSM3417861")

samples_1 = c(subseries_2a, subseries_2b, subseries_2c, subseries_2d, subseries_2e, subseries_2f)
samples_2 = c(subseries_3, subseries_3b, subseries_3c)
samples_3 = subseries_4

#series with only case and no controls
list_of_series = c("GSE93988","GSE107692","GSE139963")

#gene set list for hypoxia
gs_hypoxia <- read.table("raw_data/gene_sets/gene_set_hypoxia.txt", sep = ",", skip = 3)
colnames(gs_hypoxia) = c("genes")
gs_tgfb <- read.table("raw_data/gene_sets/gene_set_tgfbeta.txt", sep= ",", skip = 3)
colnames(gs_tgfb) = c("genes")
gene_set_list_hypoxia <- list(as.vector(gs_hypoxia$genes), as.vector(gs_tgfb$genes))
names(gene_set_list_hypoxia) = c("hypoxia" ,"tgfbeta")



