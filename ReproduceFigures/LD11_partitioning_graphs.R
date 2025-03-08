# Clear the workspace
rm(list=ls())

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source files
source('../ParameterOptimization/soybean_parameter_expansion.R')
co2_opt = '_tbd_'
years <- c('2021', '2022', '2023')
Catms <- c(414.7, 417.2, 419.3) # from NOAA
source('../Data/Soybean-BioCro_Parameters/soybean_parameters.R')
source('../Data/Soybean-BioCro_Parameters/soybean_modules.R')
# Load packages
library(BioCro)
library(lattice)
library(reshape2)
library(ggplot2)

# Load shared data
load('../Data/Soybean-BioCro_Parameters/full_soybean_ld11.RData')

# Set up basic properties for the solver
full_soybean_ld11$ode_solver <- solver

ExpBiomass <- list()
weather.afteremergence <- list()
results <- list()
figs <- list()
lai.figs <- list()
allocation.figs <- list()
partitioning.figs <- list()

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

plot_partitioning <- function(res, year) {
  
  # Colorblind friendly color palette (https://personal.sron.nl/~pault/)
  col.palette.muted <- c("#332288", "#117733", "#999933", "#882255")
  
  size.title <- 12
  size.axislabel <-10
  size.axis <- 10
  size.legend <- 12
  
  r.partcoeffs.doy <- reshape2::melt(res[,c("DVI","kRoot","kLeaf","kStem","kGrain")],id.vars="DVI")
  r.partcoeffs.doy$value<-100*r.partcoeffs.doy$value
  f <- ggplot() + theme_classic()
  f <- f + geom_point(data=r.partcoeffs.doy, aes(x=DVI, y=value, colour=variable),show.legend = TRUE,size=0.25)
  f <- f + labs(title=element_blank(), x='DVI',y='% Allocated')
  f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
                 axis.text=element_text(size=size.axis),
                 axis.title=element_text(size=size.axislabel),
                 legend.position = c(.15,.45), legend.title = element_blank(),
                 legend.text=element_text(size=size.legend),
                 legend.background = element_rect(fill = "transparent",colour = NA),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
                 plot.background = element_rect(fill = "transparent", colour = NA))
  f <- f + scale_x_continuous(breaks = seq(-1, 2, 1))
  f <- f + guides(colour = guide_legend(override.aes = list(size=2)))
  f <- f + scale_colour_manual(values = col.palette.muted,labels=c("Root","Leaf","Stem","Pod"))
  
  return(f)
}

for (i in 1:length(years)){
  ExpBiomass[[i]] <- loadRData(paste0('../../energy-farm-biocro/soybean_ld11_biomass_', years[i],'/soybean_ld11_biomass_', years[i], '.RData'))
  weather <- loadRData(paste0('../../energy-farm-biocro/weather_', years[i], '/weather', years[i], '_hourly.RData'))
  full_soybean_ld11$parameters$Catm <- Catms[i]
  
  first_data_time <- ExpBiomass[[i]]$time[1]
  weather.aftersowing <- weather[weather$time >= first_data_time, ]
  
  soybean$parameters$b0 <- full_soybean_ld11$parameters$b0
  soybean$parameters$b1 <- full_soybean_ld11$parameters$b1
  
  soybean$parameters$vmax1 <- full_soybean_ld11$parameters$vmax1
  soybean$parameters$tpu_rate_max <-full_soybean_ld11$parameters$tpu_rate_max
  soybean$parameters$jmax <- full_soybean_ld11$parameters$jmax
  soybean$parameters$Rd <- 1.28
  
  # Obtain DVI from the original Soybean-BioCro
  result <- run_biocro(soybean$initial_values,
                       soybean$parameters,
                       weather.aftersowing,
                       soybean$direct_modules,
                       soybean$differential_modules,
                       soybean$ode_solver)

  
  results[[i]] <- result
  
  # Save plot into the list
  names(result)[names(result) == "Grain"] <- "Pod"
  biocro_organ_biomass <- result[c('time', 'Leaf', 'Stem', 'Root', 'Pod')]
  biocro_organ_biomass_tall <- melt(biocro_organ_biomass, id.vars = 'time')
  names(biocro_organ_biomass_tall) <- c('time','Organ', 'biomass')
  
  field_organ_biomass <- ExpBiomass[[i]][c('time', 'leaf', 'stem', 'root', 'pod')]
  names(field_organ_biomass)[names(field_organ_biomass) %in% c('leaf', 'stem', 'root', 'pod')] <- c('Leaf', 'Stem', 'Root', 'Pod')
  field_organ_biomass_tall <- melt(field_organ_biomass, id.vars = 'time')
  names(field_organ_biomass_tall) <- c('time','Organ', 'biomass')
  
  size.title <- 12
  size.axislabel <-12
  size.axis <- 12
  size.legend <- 8
  
  col.palette.muted <- c( "#117733", "#999933", "#332288", "#882255")
  
  figs[[i]] <- ggplot() + theme_classic() +
    geom_line(data = biocro_organ_biomass_tall, aes(x = time, y = biomass, color = Organ), linewidth = 1) +
    geom_point(data = field_organ_biomass_tall, aes(x = time, y = biomass, color = Organ), shape = 15, size = 3)+
    theme(plot.title=element_text(size=size.title, hjust=0.5),
          axis.text=element_text(size=size.axis),
          axis.title.x =element_text(size=size.axislabel),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA))+
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
    scale_x_continuous(breaks = seq(180,280,30))+
    labs(title=element_blank(), 
         x=paste0('Day of Year (', years[i], ')'), 
         y='Biomass (Mg/ha)')+
    scale_color_manual(values = col.palette.muted)
  
  # if (years[i]=='2022'){
  #   figs[[i]] <- figs[[i]] + theme(plot.background = element_rect(fill = "grey90", colour = NA))
  # }
  
  print(figs[[i]])
  save(biocro_organ_biomass_tall, file = paste0('organ_biomass_sim_', years[i],'_ld11_partitioning.RData'))
  
  # lai plots
  lai.figs[[i]] <- ggplot() + theme_classic() +
    geom_line(data = result, aes(x = time, y = lai), linewidth = 1) +
    geom_point(data = ExpBiomass[[i]], aes(x = time, y = LAI_from_LMA), shape = 15, size = 3)+
    theme(plot.title=element_text(size=size.title, hjust=0.5),
          axis.text=element_text(size=size.axis),
          axis.title.x =element_text(size=size.axislabel),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA))+
    # scale_y_continuous(limits = c(0, 7), breaks = seq(0, 7, 2)) +
    # scale_x_continuous(breaks = seq(180,280,30))+
    labs(title=element_blank(), 
         x=paste0('Day of Year (', years[i], ')'), 
         y='LAI')+
    scale_color_manual(values = col.palette.muted)
  
  # print(lai.figs[[i]])
  partitioning.figs[[i]] <- plot_partitioning(result, years[i])
}
print(partitioning.figs[[1]])
# library(grid)
# library(gridExtra)
# g_legend <-function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}
# 
# common_legend <- g_legend(figs[[1]])
# 
# combined_graph <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90, gp=gpar(fontsize=12))),
#                                arrangeGrob(arrangeGrob(figs[[1]] + theme(legend.position="none"), top = "Testing"),
#                                            arrangeGrob(figs[[2]] + theme(legend.position="none"), top = "Testing"),
#                                            arrangeGrob(figs[[3]] + theme(legend.position="none"), top = "Testing"),
#                                            ncol = 3, top = "LD11 at Energy Farm (Ambient CO2)"),
#                                common_legend, 
#                                ncol=3, widths=c(0.3, 5, 1.1))
# 
# common_legend <- g_legend(partitioning.figs[[1]])
# 
# 
# # combined_graph.allocation <- grid.arrange(arrangeGrob(textGrob('Remoblized C %', rot = 90, gp=gpar(fontsize=12))),
# #                                           arrangeGrob(arrangeGrob(allocation.figs[[1]] + theme(legend.position="none"),
# #                                                                   allocation.figs[[2]] + theme(legend.position="none"),
# #                                                                   allocation.figs[[3]] + theme(legend.position="none"),
# #                                                                   ncol = 3),
# #                                                       ncol = 1),
# #                                           ncol=2, widths=c(0.3, 5))
# 
# combined_graph.lai <- grid.arrange(arrangeGrob(textGrob('LAI', rot = 90, gp=gpar(fontsize=12))),
#                                           arrangeGrob(arrangeGrob(lai.figs[[1]] + theme(legend.position="none"),
#                                                                   lai.figs[[2]] + theme(legend.position="none"),
#                                                                   lai.figs[[3]] + theme(legend.position="none"),
#                                                                   ncol = 3),
#                                                       ncol = 1),
#                                           ncol=2, widths=c(0.3, 5))
# 
# # for(i in 1:3){
# #   print(xyplot(data=weather.afteremergence[[i]], temp~time, 
# #                scales = list(x = list(at = seq(150, 300, by = 20))),
# #                ylim = c(-5,40),
# #                main = years[i],
# #                pch = 16,
# #                cex = 0.5,
# #                xlab = 'DOY',
# #                ylab = 'Temperature (degree Celcius)',
# #                panel = function(...) {
# #                   panel.xyplot(...) # Default xyplot panel
# #                   panel.abline(h = 15, col = "orange", lty = 2)
# #                   panel.abline(h = 8, col = "red", lty = 2)
# #                }))
# # }
# # 
# r <- results[[2]]
# times <- c(186.5, 209.5, 236.5, 258.5)
# doys <- c(186, 209, 236, 258)
# layer_assim <- data.frame(DOY=numeric(),
#                           layer_number=numeric(),
#                           layer_assimilation=numeric())
# for (t in 1:length(times)){
#   time <- times[t]
#   doy <- doys[t]
#   idx <- which(r$time==time)
#   for (i in 1:9){
#     assim <- (r[idx, paste0('sunlit_Assim_layer_', i)]*
#               r[idx, paste0('sunlit_fraction_layer_', i)] +
#               r[idx, paste0('shaded_Assim_layer_', i)]*
#               r[idx, paste0('shaded_fraction_layer_', i)]) *
#       r[idx, 'lai'] /10 
#     new_row <- data.frame(
#       DOY = doy,
#       layer_number = i,
#       layer_assimilation = assim)
#     layer_assim <- rbind(layer_assim, new_row)
#   }
# }
# 
# plot <- ggplot(layer_assim, aes(x = layer_number, y = layer_assimilation)) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(~ DOY, nrow = 1, scales = "fixed") +
#   labs(# title = "Layer Assimilation by Layer Number on Different DOYs",
#        x = "Layer Number",
#        y = "Layer Assimilation (micromol / s)") +
#   theme_minimal() +
#   theme(
#     strip.background = element_rect(fill = "lightgrey"),
#     strip.text = element_text(face = "bold"),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.major.x = element_line(color = "grey90")
#   ) +
#   scale_x_continuous(breaks = unique(layer_assim$layer_number), 
#                      labels = as.integer(unique(layer_assim$layer_number)))
# 
# # Display the plot
# print(plot)
# 
# xyplot(data=ExpBiomass[[1]], LAI_from_LMA~time)
