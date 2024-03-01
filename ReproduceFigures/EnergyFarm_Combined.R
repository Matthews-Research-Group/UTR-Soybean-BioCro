library(BioCro)
library(ggplot2)
library(grid)
library(gridExtra)

# Clear workspace
rm(list=ls())
# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# # Source 2021, 2022 files
# source('2021_LD11_graphs.R')
# source('2022_LD11_graphs.R')
# source('2022_LD10_graphs.R')
# Load Data
load('organ_biomass_plot_2021_ld11.RData')
load('organ_biomass_plot_2022_ld11.RData')
load('organ_biomass_plot_2022_ld10.RData')

load('organ_biomass_plot_2021_ld11_partitioning.RData')
load('organ_biomass_plot_2022_ld11_partitioning.RData')
load('organ_biomass_plot_2022_ld10_partitioning.RData')

years <- c('2021', '2022')
figs.thornley <- list()
figs.partitioning <- list()
plot_biomass <- function(year, cultivar, biocro_organ_biomass_tall, field_organ_biomass_tall, model = "partitioning"){
  # Reorganize data
  biocro_organ_biomass_tall$Source = "Simulated"
  field_organ_biomass_tall$Source = "Observed"
  r.all <- rbind(biocro_organ_biomass_tall, field_organ_biomass_tall)
  
  size.title <- 10
  size.axislabel <- 10
  size.axis <- 10
  size.legend <- 12
  
  # Colorblind friendly color palette (https://personal.sron.nl/~pault/)
  col.palette.muted <- c( "#117733", "#999933",  "#882255", "#332288")
  
  f <- ggplot() + theme_classic()
  
  f <- f + geom_point(data=r.all, aes(x=time, y=biomass,
                                      color=Organ,
                                      size = Source, shape = Source),
                      show.legend = TRUE, stroke=0.5) +
    scale_shape_manual(values = c(15, 16)) +
    scale_size_manual(values = c(2, 0.25)) +
    scale_color_manual(values = col.palette.muted)
  
  # change the plot labels and theme
  f <- f + labs(title=cultivar, x=paste0('Day of Year (',year,')'),y=NULL)
  if (cultivar == 'LD11, Training'){
    f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
                   axis.text=element_text(size=size.axis),
                   axis.title=element_text(size=size.axislabel),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), 
                   panel.background = element_rect(fill = "transparent",colour = NA),
                   plot.background = element_rect(fill = "grey90", colour = NA))
  }else{
    f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
                 axis.text=element_text(size=size.axis),
                 axis.title=element_text(size=size.axislabel),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
                 plot.background = element_rect(fill = "transparent", colour = NA))
  }
  f <- f + scale_x_continuous(breaks = seq(150,280,30))
  

  return (f)
}

figs.thornley[[1]] <- plot_biomass(years[[1]], 'LD11, Testing',
                                   biocro_organ_biomass_tall_2021_ld11, 
                                   field_organ_biomass_tall_2021_ld11,
                                   model = 'thornley')
figs.thornley[[2]] <- plot_biomass(years[[2]], 'LD11, Training',
                                   biocro_organ_biomass_tall_2022_ld11, 
                                   field_organ_biomass_tall_2022_ld11,
                                   model = 'thornley')
figs.thornley[[3]] <- plot_biomass(years[[2]], 'LD10, Testing',
                                   biocro_organ_biomass_tall_2022_ld10, 
                                   field_organ_biomass_tall_2022_ld10,
                                   model = 'thornley')

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

common_legend <- g_legend(figs.thornley[[1]])

combined_graph <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90, gp=gpar(fontsize=12))),
                               arrangeGrob(arrangeGrob(figs.thornley[[1]] + theme(legend.position="none"),
                                                       figs.thornley[[2]] + theme(legend.position="none"),
                                                       figs.thornley[[3]] + theme(legend.position="none"),
                                                       ncol = 3),
                                           ncol = 1),
                               common_legend, 
                               ncol=3, widths=c(0.3, 5, 1.1))


