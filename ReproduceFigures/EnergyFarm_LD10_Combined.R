library(BioCro)
library(ggplot2)
library(grid)
library(gridExtra)

# Clear workspace
rm(list=ls())
# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Source 2022 files
source('2022_LD10_graphs.R')
# Load Data
load('organ_biomass_plot_2022_ld10.RData')
load('organ_biomass_plot_2022_ld10_partitioning.RData')
years <- c('2022')
figs.thornley <- list()
figs.partitioning <- list()
plot_biomass <- function(year, biocro_organ_biomass_tall, field_organ_biomass_tall, model = "partitioning"){
  # Reorganize data
  biocro_organ_biomass_tall$Source = "Simulated"
  field_organ_biomass_tall$Source = "Observed"
  r.all <- rbind(biocro_organ_biomass_tall, field_organ_biomass_tall)
  
  size.title <- 12
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
  f <- f + labs(title=element_blank(), x=paste0('Day of Year (',year,')'),y=NULL)
  f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
                 axis.text=element_text(size=size.axis),
                 axis.title=element_text(size=size.axislabel),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
                 plot.background = element_rect(fill = "transparent", colour = NA))
  if (year == '2021') {
    f <- f + scale_y_continuous(limits = c(0, 8), breaks = seq(0,10,2))
  }else{
    f <- f + scale_y_continuous(limits = c(0, 10), breaks = seq(0,10,2))
  }
  if (model == "thornley"){
    f <- f + theme(axis.text.y=element_blank())
  }
  f <- f + scale_x_continuous(breaks = seq(150,280,30))
  

  return (f)
}

figs.partitioning[[1]] <- plot_biomass(years[[1]], 
                                       biocro_organ_biomass_tall_2022_ld10_partitioning, 
                                       field_organ_biomass_tall_2022_ld10_partitioning)

figs.thornley[[1]] <- plot_biomass(years[[1]], 
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
                               arrangeGrob(arrangeGrob(figs.partitioning[[1]] + theme(legend.position="none"),
                                                       nrow = 1, top = 'Partitioning'),
                                           arrangeGrob(figs.thornley[[1]] + theme(legend.position="none"),
                                                       nrow = 1, top = 'Thornley'),
                                           ncol = 2),
                               common_legend, 
                               ncol=3, widths=c(0.3, 5, 1))


