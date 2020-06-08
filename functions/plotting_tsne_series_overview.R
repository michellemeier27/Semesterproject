##FUNCTION PLOT SERIES OVERVIEW WITH TSNE
#load libraries
library(cowplot)
library(ggplot2)
library(Rtsne)

#defining function 
plotting_overview <- function(tsne_res = tsne_pi3k, meta_frame = pi3k_res_4, exp_str = "pi3k"){
  #transforming tSNE results
  new = as.data.frame(tsne_res$Y)
  #avoiding paste0 in ggsave
  title = paste0( "Overview for adapted series: ",exp_str)
  title_o = paste0("Overview for original series: ",exp_str)
  path_wanted_pic_sub = paste0("/Users/michellemeier/Semesterproject/ARCHS4/results/",exp_str, "/FINAL/",exp_str,"_subseries_overview.png")
  path_wanted_pic_o = paste0("/Users/michellemeier/Semesterproject/ARCHS4/results/",exp_str, "/FINAL/",exp_str,"_orseries_overview.png")
  path_wanted_leg_sub = paste0("/Users/michellemeier/Semesterproject/ARCHS4/results/",exp_str, "/FINAL/",exp_str,"_subseries_overview_legend.png")
  path_wanted_leg_o= paste0("/Users/michellemeier/Semesterproject/ARCHS4/results/",exp_str, "/FINAL/",exp_str,"_orseries_overview_legend.png")
  #plot and legend for subseries
  plot_s <- ggplot(new, aes(x= V1, y= V2, colour= meta_frame$series)) + geom_point(size = 0.5) +ggtitle(title) + labs(color='adapted series')
  plot_1s <- ggplot(new, aes(x= V1, y= V2, colour= meta_frame$series)) + geom_point(size = 0.5) +ggtitle(title) + theme(legend.position = "none")
  ggsave(path_wanted_pic_sub, plot = plot_1s)
  legend <- cowplot::get_legend(plot_s)
  grid.newpage()
  grid.draw(legend)
  ggsave(path_wanted_leg_sub, plot = legend)
  #plot and legend for orginal series 
  plot_o <- ggplot(new, aes(x= V1, y= V2, colour= meta_frame$original_series)) + geom_point(size = 0.5) +ggtitle(title_o) + labs(color = "original series")
  plot_1o <- ggplot(new, aes(x= V1, y= V2, colour= meta_frame$original_series)) + geom_point(size = 0.5) +ggtitle(title_o) + theme(legend.position = "none")
  ggsave(path_wanted_pic_o, plot = plot_1o)
  legend <- cowplot::get_legend(plot_o)
  grid.newpage()
  grid.draw(legend)
  ggsave(path_wanted_leg_o, plot = legend)
}
