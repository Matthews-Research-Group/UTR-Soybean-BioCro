# Clear workspace
rm(list=ls())

load("SoyFACE_results_and_measurements.RData")
source('plot_partitioning.R')
co2_opt <- '_ambient_'
source('../Data/Soybean-BioCro_Parameters/soybean_parameters.R')
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')
source('../ParameterOptimization/soybean_parameter_expansion.R')
# Load shared data
load('../Data/Soybean-BioCro_Parameters/full_soybean_ld11.RData')

# Update UTR parameters
fitted.utr.params <- optim_params_conversion(optim_params_short)
names(fitted.utr.params) <- arg_names
parameters <-c(parameters, fitted.utr.params)[!duplicated(c(names(parameters), 
                                                            names(fitted.utr.params)), 
                                                          fromLast = TRUE)]


years <- c('2002', '2004', '2005', '2006')
figs <- list()
for (i in 1:length(results)){
  figs[[i]] <- plot_partitioning(results[[i]], years[i])
}

library(grid)
library(gridExtra)
g_legend <-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# common_legend <- g_legend(figs[[1]])
# combined_graph.allocation <- grid.arrange(arrangeGrob(textGrob('Allocation %', rot = 90, gp=gpar(fontsize=12))),
#                                           arrangeGrob(arrangeGrob(figs[[1]] + theme(legend.position="none"),
#                                                                   figs[[2]] + theme(legend.position="none"),
#                                                                   figs[[3]] + theme(legend.position="none"),
#                                                                   figs[[4]] + theme(legend.position="none"),
#                                                                   ncol = 4),
#                                                       ncol = 1),
#                                           common_legend, 
#                                           ncol=3, widths=c(0.3, 5, 1.1))

combined_graph.allocation <- grid.arrange(arrangeGrob(textGrob('Remoblized C %', rot = 90, gp=gpar(fontsize=12))),
                                          arrangeGrob(arrangeGrob(figs[[1]] + theme(legend.position="none"),
                                                                  figs[[2]] + theme(legend.position="none"),
                                                                  figs[[3]] + theme(legend.position="none"),
                                                                  figs[[4]] + theme(legend.position="none"),
                                                                  ncol = 4),
                                                      ncol = 1),
                                          ncol=2, widths=c(0.3, 5))
