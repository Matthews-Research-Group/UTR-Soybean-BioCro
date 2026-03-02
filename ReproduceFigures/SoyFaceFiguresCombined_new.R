library(BioCro)
library(UTRSoybeanBML)
library(BioCroWater)
library(ggplot2)
library(grid)
library(gridExtra)
library(lattice)

# Clear workspace
rm(list=ls())

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Parameters
source('../ParameterOptimization/soybean_parameter_expansion.R')
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')

# Source mse calculation and partitioning plotting functions
source('calculate_mse.R')
source('plot_partitioning.R')

# years, sowing dates, and harvesting dates of growing seasons being fit to
years <- c('2002', '2004', '2005', '2006')
sow.date <- c(152, 149, 148, 148)
harv.date <- c(288, 289, 270, 270)
co2_opt = '_ambient_'
# Model set up
initial_values <- soybean$initial_values
direct_modules <- soybean$direct_modules
differential_modules <- soybean$differential_modules
solver <- soybean$ode_solver
parameters <- soybean$parameters

# Set up BioCroWater
source('../Data/Soybean-BioCro_Parameters/set_up_BioCroWater.R')
initial_values       <- set_init_values(initial_values)
parameters           <- set_parameters(parameters)
parameters$kd        <- parameters$k_diffuse
direct_modules       <- set_direct_modules(direct_modules) 
differential_modules <- set_differential_modules(differential_modules) 

# Update UTR modules
source('../Data/Soybean-BioCro_Parameters/set_up_UTRSoybeanBML.R')
initial_values <- set_init_values(initial_values, ExpBiomass[[i]]) 
parameters           <- set_parameters(parameters)
parameters$Catm      <- 372
direct_modules       <- set_direct_modules(direct_modules) 
differential_modules <- set_differential_modules(differential_modules) 

# Update UTR parameters
fitted.utr.params <- optim_params_conversion(optim_params_short_SoyFACE)
names(fitted.utr.params) <- arg_names
parameters <-c(parameters, fitted.utr.params)[!duplicated(c(names(parameters), 
                                                            names(fitted.utr.params)), 
                                                          fromLast = TRUE)]

# Minor adjustments
parameters$time_zone_offset <- -6
initial_values$DVI <- 0

# Initialize lists
results <- list()
results.elevCO2 <- list()
weather.growingseason <- list()
soybean_optsolver <- list()
ExpBiomass <- list()
ExpBiomass.elevCO2 <- list()
LAI <- list()
LAI.elevCO2 <- list()
ExpBiomass.std <- list()
RootVals <- list()
numrows <- vector()
figs <- list()
lai.figs <- list()
pod_sensitivity <- data.frame()
allocation.figs <- list()
allocation.elevCO2.figs <- list()

# Define functions to create plots
plot_amb_elev_lai <- function(res, elev_res, year, lai, elev_lai) {
  
  # Colorblind friendly color palette (https://personal.sron.nl/~pault/)
  col.palette.muted <- c("#332288", "#117733", "#999933", "#882255")
  
  size.title <- 12
  size.axislabel <-10
  size.axis <- 10
  size.legend <- 12
  
  res$lai_tot <- (res$Leaf_substrate_carbon+res$Leaf_structural_carbon)* 0.3 * parameters$iSp
  elev_res$lai_tot <- (elev_res$Leaf_substrate_carbon+elev_res$Leaf_structural_carbon)* 0.3 * parameters$iSp
  s.lai <- cbind(res[,c("fractional_doy","lai_tot")],elev_res[,"lai_tot"])
  colnames(s.lai) <- c("fractional_doy","Amb","Elev")
  r.lai <- reshape2::melt(s.lai, id.vars = "fractional_doy")
  
  s.exp.lai <- cbind(lai[,c("DOY","LAI_mean")],elev_lai[,"LAI_mean"])
  colnames(s.exp.lai) <- c("DOY","Amb","Elev")
  r.exp.lai <- reshape2::melt(s.exp.lai, id.vars = "DOY")
  
  s.exp.std.lai <- cbind(lai[,c("DOY","LAI_std")],elev_lai[,"LAI_std"])
  colnames(s.exp.std.lai) <- c("DOY","Amb","Elev")
  r.exp.std.lai <- reshape2::melt(s.exp.std.lai, id.vars = "DOY")
  r.exp.std.lai$ymin <- r.exp.lai$value - r.exp.std.lai$value
  r.exp.std.lai$ymax <- r.exp.lai$value + r.exp.std.lai$value
  
  f <- ggplot() + theme_classic()
  f <- f + geom_point(data=r.lai, aes(x=fractional_doy, y=value, colour=variable),show.legend = FALSE,size=0.25)
  f <- f + geom_errorbar(data=r.exp.std.lai, aes(x=DOY, ymin=ymin, ymax=ymax), width=3.5, size=0.25, show.legend = FALSE)
  f <- f + geom_point(data=r.exp.lai, aes(x=DOY, y=value, fill=variable), shape=22, size=2, show.legend = FALSE, stroke=.5)
  f <- f + coord_cartesian(ylim = c(0,10)) + scale_x_continuous(breaks = seq(150,275,30))
  f <- f + labs(x='Day of Year' ,y=bquote("LAI"~(m^2~"/"~m^2)))
  f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
                 axis.text=element_text(size=size.axis),
                 axis.title=element_text(size=size.axislabel),
                 axis.title.y = element_blank(),
                 legend.position = c(.25,.85), legend.title = element_blank(),
                 legend.text=element_text(size=size.legend),
                 legend.background = element_rect(fill = "transparent",colour = NA),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
                 plot.background = element_rect(fill = "transparent", colour = NA))
  f <- f + guides(colour = guide_legend(override.aes = list(size=2)))
  f <- f + scale_fill_manual(values = col.palette.muted[2:3], guide = "none")
  f <- f + scale_colour_manual(values = col.palette.muted[2:3],labels=c('Ambient',bquote(Elevated~CO[2])))
  
  return(f)
}

for (i in 1:length(years)) {   
  yr <- years[i]
  weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength_wDVI.csv'))
  em.ind <- which(weather$DVI > 0)[1] # emergence doy
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather.growingseason[[i]] <- weather[em.ind:hd.ind,]
  weather.growingseason[[i]]$DVI <- NULL
  
  ExpBiomass[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass.csv'))
  colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  ExpBiomass.std[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[3], "Root"=0.17*sum(ExpBiomass[[i]][5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  numrows[i] <- nrow(weather.growingseason[[i]])
  invwts <- ExpBiomass.std[[i]]
  
  # Load LAI data
  if(i>1){
    LAI[[i]]<-read.csv(file=paste0('../Data/SoyFACE_data/',yr,'_ambient_lai.csv'))
    LAI.elevCO2[[i]]<-read.csv(file=paste0('../Data/SoyFACE_data/',yr,'_elevated_lai.csv'))
  }
  
  soybean_optsolver[[i]] <- partial_run_biocro(initial_values,
                                               parameters,
                                               weather.growingseason[[i]],
                                               direct_modules,
                                               differential_modules,
                                               solver,
                                               arg_names,
                                               verbose = FALSE)
  
  result <- soybean_optsolver[[i]](optim_params_conversion(optim_params_short_SoyFACE)) 
  results[[i]] <- result
  
  # organize simulated data
  r.lsrp.doy <- reshape2::melt(result[,c("fractional_doy","Root","Leaf","Stem","Pod")],id.vars="fractional_doy")
  
  # Leaf
  # organize the experimental data (mean and std)
  s.exp.leaf <- cbind(ExpBiomass[[i]][,c("DOY","Leaf")])
  colnames(s.exp.leaf) <- c("fractional_doy","Leaf") # DOY renamed as fractional_doy
  r.exp.leaf <- reshape2::melt(s.exp.leaf, id.vars = "fractional_doy") # DOY renamed as fractional_doy
  
  s.exp.std.leaf <- cbind(ExpBiomass.std[[i]][,c("DOY","Leaf")])
  colnames(s.exp.std.leaf) <- c("fractional_doy","Leaf") # DOY renamed as fractional_doy
  r.exp.std.leaf <- reshape2::melt(s.exp.std.leaf, id.vars = "fractional_doy") # DOY renamed as fractional_doy
  r.exp.std.leaf$ymin <- r.exp.leaf$value - r.exp.std.leaf$value
  r.exp.std.leaf$ymax <- r.exp.leaf$value + r.exp.std.leaf$value
  
  # Stem
  s.exp.stem <- cbind(ExpBiomass[[i]][,c("DOY","Stem")]) # DOY renamed as fractional_doy
  colnames(s.exp.stem) <- c("fractional_doy","Stem") # DOY renamed as fractional_doy
  r.exp.stem <- reshape2::melt(s.exp.stem, id.vars = "fractional_doy") # DOY renamed as fractional_doy
  
  s.exp.std.stem <- cbind(ExpBiomass.std[[i]][,c("DOY","Stem")]) # DOY renamed as fractional_doy
  colnames(s.exp.std.stem) <- c("fractional_doy","Stem") # DOY renamed as fractional_doy
  r.exp.std.stem <- reshape2::melt(s.exp.std.stem, id.vars = "fractional_doy") # DOY renamed as fractional_doy
  r.exp.std.stem$ymin <- r.exp.stem$value - r.exp.std.stem$value
  r.exp.std.stem$ymax <- r.exp.stem$value + r.exp.std.stem$value
  
  # Pod
  s.exp.pod <- cbind(ExpBiomass[[i]][,c("DOY","Pod")]) # DOY renamed as fractional_doy
  colnames(s.exp.pod) <- c("fractional_doy","Pod") # DOY renamed as fractional_doy
  r.exp.pod <- reshape2::melt(s.exp.pod, id.vars = "fractional_doy") # DOY renamed as fractional_doy
  
  s.exp.std.pod <- cbind(ExpBiomass.std[[i]][,c("DOY","Pod")]) # DOY renamed as fractional_doy
  colnames(s.exp.std.pod) <- c("fractional_doy","Pod") # DOY renamed as fractional_doy
  r.exp.std.pod <- reshape2::melt(s.exp.std.pod, id.vars = "fractional_doy") # DOY renamed as fractional_doy
  r.exp.std.pod$ymin <- r.exp.pod$value - r.exp.std.pod$value
  r.exp.std.pod$ymax <- r.exp.pod$value + r.exp.std.pod$value
  
  # Combine
  r.exp.ls <- rbind(r.exp.leaf, r.exp.stem, r.exp.pod)
  r.exp.ls$Source = "Observed"
  r.lsrp.doy$Source = "Simulated"
  r.all <- rbind(r.lsrp.doy, r.exp.ls)
  # Reverse the order as follow
  r.all$Organ <- factor(r.all$variable, levels = rev(levels(r.all$variable)))
  # Colorblind friendly color palette (https://personal.sron.nl/~pault/)
  col.palette.muted.organs <- c("Leaf"="#117733", 
                                "Stem"= "#999933", 
                                "Root"="#332288", 
                                "Pod"= "#882255")
  size.title <- 12
  size.axislabel <- 10
  size.axis <- 10
  size.legend <- 12
  
  f <- ggplot() + theme_classic()
  f <- f +
    geom_line(data = subset(r.all, Source == "Simulated"),  
                aes(x=fractional_doy,y=value, color=Organ), size=0.8, alpha = 0.8) +
    geom_point(data = subset(r.all, Source == "Observed"), 
                 aes(x=fractional_doy, y=value, color=Organ), shape=15, size=2, stroke=.5) +
  
    scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, 2)) +
    scale_color_manual(values = col.palette.muted.organs)

  
  # for leaf
  f <- f + geom_errorbar(data=r.exp.std.leaf, aes(x=fractional_doy, ymin=ymin, ymax=ymax),  # DOY renamed as fractional_doy
                         width=3.5, size=0.25, show.legend = FALSE)
  # for stem
  f <- f + geom_errorbar(data=r.exp.std.stem, aes(x=fractional_doy, ymin=ymin, ymax=ymax),   # DOY renamed as fractional_doy
                         width=3.5, size=0.25, show.legend = FALSE)
  
  # for pod
  f <- f + geom_errorbar(data=r.exp.std.pod, aes(x=fractional_doy, ymin=ymin, ymax=ymax),   # DOY renamed as fractional_doy
                         width=3.5, size=0.25, show.legend = FALSE)
  
  # change the plot labels and theme
  f <- f + labs(title=years[i], x='Day of Year', y=NULL)
  if (years[i] == '2002' || years[i] == '2005'){
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
  figs[[i]] <- f
  
  names(r.lsrp.doy) <- c('time', 'Organ', 'biomass', 'Source')
  names(r.exp.ls) <- c('time', 'Organ', 'biomass', 'Source')
  calculate_mse_soyFACE(years[[i]],r.lsrp.doy, r.exp.ls)
  allocation.figs[[i]] <- plot_partitioning(result, years[i])
}


# initialize lists for figures
figs.elevCO2 <- list()
co2_opt = '_co2_'
parameters$Catm      <- 550

# Initialize lists
weather.growingseason <- list()
soybean_optsolver <- list()
ExpBiomass.std <- list()
RootVals <- list()
numrows <- vector()

for (i in 1:length(years)) {  
  yr <- years[i]
  weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength_wDVI.csv'))
  # sd.ind <- which(weather$doy == sow.date[i])[1]
  em.ind <- which(weather$DVI > 0)[1] # emergence doy
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather.growingseason[[i]] <- weather[em.ind:hd.ind,]
  weather.growingseason[[i]]$DVI <- NULL
  
  ExpBiomass.elevCO2[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass.csv'))
  colnames(ExpBiomass.elevCO2[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  ExpBiomass.std[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  RootVals[[i]] <- data.frame("DOY"=ExpBiomass.elevCO2[[i]]$DOY[3], "Root"=0.17*sum(ExpBiomass.elevCO2[[i]][5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  numrows[i] <- nrow(weather.growingseason[[i]])
  invwts <- ExpBiomass.std[[i]]

  soybean_optsolver[[i]] <- partial_run_biocro(initial_values,
                                               parameters,
                                               weather.growingseason[[i]],
                                               direct_modules,
                                               differential_modules,
                                               solver,
                                               arg_names,
                                               verbose = FALSE)
  
  result <- soybean_optsolver[[i]](optim_params_conversion(optim_params_short_SoyFACE)) 
  results.elevCO2[[i]] <- result
  # check when the simulation stops if not running till the end
  if (dim(result)[1] < dim(weather.growingseason[[i]])[1]){
    print(max(result$DVI))
    print(result$Stem_substrate_carbon[which.max(result$DVI)])
    print(result$Pod_substrate_carbon[which.max(result$DVI)])
    
    print(result$Pod_utilization_rate[which.max(result$DVI)])
    print(result$substrate_transport_Stem_to_Pod[which.max(result$DVI)])
    
    print(result$Stem_substrate_carbon[which.max(result$DVI)]/result$Stem[which.max(result$DVI)])
    print(result$Pod_substrate_carbon[which.max(result$DVI)]/result$Pod[which.max(result$DVI)])
    
    print(result$Pod[which.max(result$DVI)])
  }
  # organize simulated data
  r.lsrp.doy <- reshape2::melt(result[,c("fractional_doy","Root","Leaf","Stem","Pod")],id.vars="fractional_doy")
  r.lsrp.doy$value<-r.lsrp.doy$value
  
  # Leaf
  # organize the expirimental data (mean and std)
  s.exp.leaf <- cbind(ExpBiomass.elevCO2[[i]][,c("DOY","Leaf")])
  colnames(s.exp.leaf) <- c("fractional_doy","Leaf") # DOY renamed as fractional_doy
  r.exp.leaf <- reshape2::melt(s.exp.leaf, id.vars = "fractional_doy") # DOY renamed as fractional_doy
  
  s.exp.std.leaf <- cbind(ExpBiomass.std[[i]][,c("DOY","Leaf")])
  colnames(s.exp.std.leaf) <- c("fractional_doy","Leaf") # DOY renamed as fractional_doy
  r.exp.std.leaf <- reshape2::melt(s.exp.std.leaf, id.vars = "fractional_doy") # DOY renamed as fractional_doy
  r.exp.std.leaf$ymin <- r.exp.leaf$value - r.exp.std.leaf$value
  r.exp.std.leaf$ymax <- r.exp.leaf$value + r.exp.std.leaf$value
  
  # Stem
  s.exp.stem <- cbind(ExpBiomass.elevCO2[[i]][,c("DOY","Stem")]) # DOY renamed as fractional_doy
  colnames(s.exp.stem) <- c("fractional_doy","Stem") # DOY renamed as fractional_doy
  r.exp.stem <- reshape2::melt(s.exp.stem, id.vars = "fractional_doy") # DOY renamed as fractional_doy
  
  s.exp.std.stem <- cbind(ExpBiomass.std[[i]][,c("DOY","Stem")]) # DOY renamed as fractional_doy
  colnames(s.exp.std.stem) <- c("fractional_doy","Stem") # DOY renamed as fractional_doy
  r.exp.std.stem <- reshape2::melt(s.exp.std.stem, id.vars = "fractional_doy") # DOY renamed as fractional_doy
  r.exp.std.stem$ymin <- r.exp.stem$value - r.exp.std.stem$value
  r.exp.std.stem$ymax <- r.exp.stem$value + r.exp.std.stem$value
  
  # Pod
  s.exp.pod <- cbind(ExpBiomass.elevCO2[[i]][,c("DOY","Pod")]) # DOY renamed as fractional_doy
  colnames(s.exp.pod) <- c("fractional_doy","Pod") # DOY renamed as fractional_doy
  r.exp.pod <- reshape2::melt(s.exp.pod, id.vars = "fractional_doy") # DOY renamed as fractional_doy
  
  s.exp.std.pod <- cbind(ExpBiomass.std[[i]][,c("DOY","Pod")]) # DOY renamed as fractional_doy
  colnames(s.exp.std.pod) <- c("fractional_doy","Pod") # DOY renamed as fractional_doy
  r.exp.std.pod <- reshape2::melt(s.exp.std.pod, id.vars = "fractional_doy") # DOY renamed as fractional_doy
  r.exp.std.pod$ymin <- r.exp.pod$value - r.exp.std.pod$value
  r.exp.std.pod$ymax <- r.exp.pod$value + r.exp.std.pod$value
  
  # Combine
  r.exp.ls <- rbind(r.exp.leaf, r.exp.stem, r.exp.pod)
  r.exp.ls$Source = "Observed"
  r.lsrp.doy$Source = "Simulated"
  r.all <- rbind(r.lsrp.doy, r.exp.ls)
  # Reverse the order as follow
  r.all$Organ <- factor(r.all$variable, levels = rev(levels(r.all$variable)))
  
  # Colorblind friendly color palette (https://personal.sron.nl/~pault/)
  col.palette.muted.organs <- c("Leaf"="#117733", 
                                "Stem"= "#999933", 
                                "Root"="#332288", 
                                "Pod"= "#882255")
  size.title <- 12
  size.axislabel <- 10
  size.axis <- 10
  size.legend <- 12
  
  f <- ggplot() + theme_classic()
  f <- f +
    geom_line(data = subset(r.all, Source == "Simulated"),  
              aes(x=fractional_doy,y=value, color=Organ), size=0.8, alpha = 0.8) +
    geom_point(data = subset(r.all, Source == "Observed"), 
               aes(x=fractional_doy, y=value, color=Organ), shape=15, size=2, stroke=.5) +
    
    scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, 2)) +
    scale_color_manual(values = col.palette.muted.organs)
  
  # for leaf
  f <- f + geom_errorbar(data=r.exp.std.leaf, aes(x=fractional_doy, ymin=ymin, ymax=ymax),  # DOY renamed as fractional_doy
                         width=3.5, size=0.25, show.legend = FALSE)
  # for stem
  f <- f + geom_errorbar(data=r.exp.std.stem, aes(x=fractional_doy, ymin=ymin, ymax=ymax),   # DOY renamed as fractional_doy
                         width=3.5, size=0.25, show.legend = FALSE)
  # for pod
  f <- f + geom_errorbar(data=r.exp.std.pod, aes(x=fractional_doy, ymin=ymin, ymax=ymax),   # DOY renamed as fractional_doy
                         width=3.5, size=0.25, show.legend = FALSE)
  
  # change the plot labels and theme
  f <- f + labs(title=years[i], x='Day of Year' ,y=NULL)
  f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
                 axis.text=element_text(size=size.axis),
                 axis.title=element_text(size=size.axislabel),
                 # legend.position = c(.15,.8), legend.title = element_blank(),
                 legend.text=element_text(size=size.legend),
                 # legend.background = element_rect(fill = "transparent",colour = NA),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
                 plot.background = element_rect(fill = "transparent", colour = NA))
  f <- f + scale_x_continuous(breaks = seq(150,280,30))
  figs.elevCO2[[i]] <- f
  
  names(r.lsrp.doy) <- c('time', 'Organ', 'biomass', 'Source')
  names(r.exp.ls) <- c('time', 'Organ', 'biomass', 'Source')
  calculate_mse_soyFACE(years[[i]],r.lsrp.doy, r.exp.ls)
  allocation.elevCO2.figs[[i]] <- plot_partitioning(result, years[i])
}

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

common_legend <- g_legend(figs[[1]])
heights <- c(0.8, 1)
combined_graph <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
                               arrangeGrob(arrangeGrob(figs[[1]] + theme(axis.title.x = element_blank(),
                                                                         axis.text.x = element_blank(),
                                                                         legend.position="none"),
                                                        figs[[2]] + theme(legend.position="none"),
                                                       nrow = 2, heights = heights),
                                            arrangeGrob(figs[[3]] + theme(axis.title.x = element_blank(),
                                                                          axis.text.x = element_blank(),
                                                                          legend.position="none"),
                                                        figs[[4]] + theme(legend.position="none"),
                                                        nrow = 2, 
                                                        heights = heights), # top = 'Testing'),
                                            ncol = 2, top = 'Ambient CO2'),
                                arrangeGrob(arrangeGrob(figs.elevCO2[[1]] + theme(axis.title.x = element_blank(),
                                                                                  axis.text.x = element_blank(),
                                                                                  legend.position="none"),
                                                        figs.elevCO2[[2]] + theme(legend.position="none"),
                                                        nrow = 2, 
                                                        heights = heights), # top = 'Testing'),
                                            arrangeGrob(figs.elevCO2[[3]] + theme(axis.title.x = element_blank(),
                                                                                  axis.text.x = element_blank(),
                                                                                  legend.position="none"),
                                                        figs.elevCO2[[4]] + theme(legend.position="none"),
                                                        nrow = 2,
                                                        heights = heights), # top = 'Testing'),
                                            ncol = 2, top = 'Elevated CO2'),
                                common_legend,
                                ncol=4, widths=c(0.3, 5,5,1.5))

ggsave('SoyFACE_UTR_model_biomass_graph.png', 
       plot = combined_graph, 
       width = 8,
       height = 4,
       units = "in",
       dpi = 300
)


combined_graph_v2 <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
                                  arrangeGrob(
                                    arrangeGrob(arrangeGrob(figs[[1]] + theme(axis.title.x = element_blank(),
                                                                              axis.text.x = element_blank(),
                                                                              legend.position="none"),
                                                            figs[[3]] + theme(legend.position="none"),
                                                            ncol = 2, left = '(A)'),
                                                top = 'Training'),
                                    
                                    arrangeGrob(arrangeGrob(figs[[2]] + theme(axis.title.x = element_blank(),
                                                                              axis.text.x = element_blank(),
                                                                              legend.position="none"),
                                                            figs[[4]] + theme(legend.position="none"),
                                                            nrow = 2, 
                                                            heights = heights, top = 'Ambient CO2'),
                                    
                                                arrangeGrob(figs.elevCO2[[1]] + theme(axis.title.x = element_blank(),
                                                                                      axis.text.x = element_blank(),
                                                                                      legend.position="none"),
                                                            figs.elevCO2[[2]] + theme(legend.position="none"),
                                                            nrow = 2, 
                                                            heights = heights, top = 'Elevated CO2'),
                                                arrangeGrob(figs.elevCO2[[3]] + theme(axis.title.x = element_blank(),
                                                                                      axis.text.x = element_blank(),
                                                                                      legend.position="none"),
                                                            figs.elevCO2[[4]] + theme(legend.position="none"),
                                                            nrow = 2,
                                                            heights = heights, top = 'Elevated CO2'),
                                                ncol = 3, left = '(B)'), 
                                    nrow = 2, heights = c(1.3,2)),
                                  common_legend,
                                  ncol=3, widths=c(0.3, 10, 1.5))

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
save(results, results.elevCO2, ExpBiomass, ExpBiomass.elevCO2, file = 'SoyFACE_results_and_measurements.RData')


for (i in 2:4){
  lai.figs[[i-1]] <- plot_amb_elev_lai(results[[i]], results.elevCO2[[i]], years[[i]], LAI[[i]], LAI.elevCO2[[i]])
  # print(lai.figs[[i-1]])
}
library(grid)
combined_graph.lai <- grid.arrange(arrangeGrob(textGrob(bquote("LAI"~(m^2~"/"~m^2)), rot = 90, gp=gpar(fontsize=12))),
                                   arrangeGrob(arrangeGrob(lai.figs[[1]] + theme(legend.position="none"),
                                                           lai.figs[[2]] + theme(legend.position="none"),
                                                           lai.figs[[3]] + theme(legend.position="none"),
                                                           ncol = 3),
                                               ncol = 1),
                                   ncol=2, widths=c(0.3, 5))

# for (i in 1:4){
#   percentage_change = (max(results.elevCO2[[i]]$Pod)-max(results[[i]]$Pod))/max(results[[i]]$Pod) * 100
#   print(percentage_change)
# }

for (i in 1:4){
  shoot_root_ratio = data.frame(
    fractional_doy = results[[i]]$fractional_doy,
    ambient_ratio = with(results[[i]], (Leaf+Stem+Pod)/Root),
    elevated_ratio = with(results.elevCO2[[i]], (Leaf+Stem+Pod)/Root))
  print(years[i])
  print(xyplot(data=shoot_root_ratio, ambient_ratio + elevated_ratio ~ fractional_doy, auto.key = TRUE))
}





# Plot LD11
years <- c('2021', '2022', '2023', '2024')
Catms <- c(414.7, 417.2, 419.3, 422.8) # from NOAA

updated_parameters <- read.csv('../Data/Soybean-BioCro_Parameters/updated-ld11-parameters.csv')
updated_idx <- match(updated_parameters$New.variable.name, names(parameters))
parameters[updated_idx] <- updated_parameters$LD11

ExpBiomass <- list()
weather.afteremergence <- list()
results <- list()
results_wBW <- list()
ld11.figs <- list()
lai.ld11.figs <- list()
allocation.ld11.figs <- list()

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


for (i in 1:length(years)){
  # for (i in 1:1){
  ExpBiomass[[i]] <- loadRData(paste0('../../energy-farm-biocro/soybean_ld11_biomass_', years[i],'/soybean_ld11_biomass_', years[i], '.RData'))
  harv.doy <- max(ExpBiomass[[i]]$doy)
  weather <- loadRData(paste0('../../energy-farm-biocro/weather_', years[i], '/weather', years[i], '_hourly.RData'))
  parameters$Catm <- Catms[i]
  # update initial values
  sub_frac <- 0.1           # substrate_fraction
  str_frac <- 1 - sub_frac  # structural_fraction
  seed_mass <- ExpBiomass[[i]]$initial_seed[1]
  j <- 2
  mass_t <- sum(ExpBiomass[[i]][j, c('leaf', 'stem', 'root')])
  leaf_frac <- ExpBiomass[[i]]$leaf[j]/mass_t
  stem_frac <- ExpBiomass[[i]]$stem[j]/mass_t
  root_frac <- ExpBiomass[[i]]$root[j]/mass_t
  
  cf <- 0.3 # optim_params_short_SoyFACE[1]
  
  updated_utr_initial_values <- list(
    Leaf_substrate_carbon = sub_frac * seed_mass * leaf_frac / cf,
    Leaf_structural_carbon = str_frac * seed_mass * leaf_frac / cf,
    Stem_substrate_carbon = sub_frac * seed_mass * stem_frac / cf,
    Stem_structural_carbon = str_frac * seed_mass * stem_frac / cf,
    Root_substrate_carbon =  sub_frac * seed_mass * root_frac / cf,
    Root_structural_carbon = str_frac * seed_mass * root_frac / cf)
  
  initial_values[names(updated_utr_initial_values)] <- updated_utr_initial_values
  
  first_data_time <- ExpBiomass[[i]]$time[1]
  weather.aftersowing <- weather[(((weather$doy-1)*24 + weather$hour) >= first_data_time) & weather$doy <= harv.doy, ]
  
  # Obtain DVI from the original Soybean-BioCro
  soybean_biocro_result <- run_biocro(soybean$initial_values,
                                      soybean$parameters,
                                      weather.aftersowing,
                                      soybean$direct_modules,
                                      soybean$differential_modules,
                                      soybean$ode_solver)
  
  weather.aftersowing$DVI <- soybean_biocro_result$DVI
  # start from emergence time
  weather.afteremergence[[i]] <- weather.aftersowing[-(1:which.min(abs(weather.aftersowing$DVI))),]
  weather.afteremergence[[i]]$DVI <- NULL
  weather.afteremergence[[i]]$time_zone_offset <- NULL
  parameters$RL_at_25 <- 1.28
  # parameters$par_energy_content <- 0.235
  result <- run_biocro(
    initial_values,
    parameters,
    weather.afteremergence[[i]],
    direct_modules,#[-1],
    differential_modules,
    solver)
  
  results[[i]] <- result
  
  # Save plot into the list
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
  
  ld11.figs[[i]] <- ggplot() + theme_classic() +
    geom_line(data = biocro_organ_biomass_tall, aes(x = time/24 + 1, y = biomass, color = Organ), linewidth = 1, alpha=0.8) +
    geom_point(data = field_organ_biomass_tall, aes(x = time/24 + 1, y = biomass, color = Organ), shape = 15, size = 3)+
    theme(plot.title=element_text(size=size.title, hjust=0.5),
          axis.text=element_text(size=size.axis),
          axis.title.x =element_text(size=size.axislabel),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA))+
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
    # scale_x_continuous(breaks = seq(180,280,30))+
    labs(title=years[i], 
         x=paste0('Day of Year'), 
         y='Biomass (Mg/ha)')+
    scale_color_manual(values = col.palette.muted)
  
  # print(ld11.figs[[i]])
  save(biocro_organ_biomass_tall, file = paste0('organ_biomass_sim_', years[i],'_ld11.RData'))
  save(field_organ_biomass_tall, file = paste0('organ_biomass_mea_', years[i],'_ld11.RData'))
  
  # lai plots
  lai.ld11.figs[[i]] <- ggplot() + theme_classic() +
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

  # print(lai.ld11.figs[[i]])
  
  allocation.ld11.figs[[i]] <- plot_partitioning(result, years[i])
  
  calculate_mse(years[[i]],biocro_organ_biomass_tall, field_organ_biomass_tall)
}

# extract the common legend
g_legend <-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

common_legend <- g_legend(ld11.figs[[1]])

combined_graph <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90, gp=gpar(fontsize=12))),
                               arrangeGrob(arrangeGrob(ld11.figs[[1]] + theme(legend.position="none"), top = "Testing"),
                                           arrangeGrob(ld11.figs[[2]] + theme(legend.position="none"), top = "Testing"),
                                           arrangeGrob(ld11.figs[[3]] + theme(legend.position="none"), top = "Testing"),
                                           arrangeGrob(ld11.figs[[4]] + theme(legend.position="none"), top = "Testing"),
                                           ncol = 4, top = "LD11 at Energy Farm (Ambient CO2)"),
                               common_legend, 
                               ncol=3, widths=c(0.3, 5, 1.1))

common_legend <- g_legend(allocation.ld11.figs[[1]])
combined_graph.allocation <- grid.arrange(arrangeGrob(textGrob('Allocation %', rot = 90, gp=gpar(fontsize=12))),
                                          arrangeGrob(arrangeGrob(allocation.ld11.figs[[1]] + theme(legend.position="none"),
                                                                  allocation.ld11.figs[[2]] + theme(legend.position="none"),
                                                                  allocation.ld11.figs[[3]] + theme(legend.position="none"),
                                                                  allocation.ld11.figs[[4]] + theme(legend.position="none"),
                                                                  ncol = 4),
                                                      ncol = 1),
                                          common_legend,
                                          ncol=3, widths=c(0.3, 5, 1.1))

combined_graph.allocation_v2 <- grid.arrange(arrangeGrob(textGrob('Allocation %', rot = 90, gp=gpar(fontsize=12))),
                                            arrangeGrob(arrangeGrob(allocation.figs[[1]] + theme(legend.position="none"),
                                                                    allocation.figs[[2]] + theme(legend.position="none"),
                                                                    allocation.figs[[3]] + theme(legend.position="none"),
                                                                    allocation.figs[[4]] + theme(legend.position="none"),
                                                                    ncol = 4, top = 'Pioneer 93B15 at ambient CO2'),
                                                        arrangeGrob(allocation.elevCO2.figs[[1]] + theme(legend.position="none"),
                                                                    allocation.elevCO2.figs[[2]] + theme(legend.position="none"),
                                                                    allocation.elevCO2.figs[[3]] + theme(legend.position="none"),
                                                                    allocation.elevCO2.figs[[4]] + theme(legend.position="none"),
                                                                    ncol = 4, top = 'Pioneer 93B15 at elevated CO2'),
                                                        arrangeGrob(allocation.ld11.figs[[1]] + theme(legend.position="none"),
                                                                    allocation.ld11.figs[[2]] + theme(legend.position="none"),
                                                                    allocation.ld11.figs[[3]] + theme(legend.position="none"),
                                                                    allocation.ld11.figs[[4]] + theme(legend.position="none"),
                                                                    ncol = 4, top = 'LD11-2170 at ambient CO2'),
                                                        nrow = 3),
                                            common_legend,
                                            ncol=3, widths=c(0.2, 5, 0.8))

ggsave('Fig3-allocation.png', 
       plot = combined_graph.allocation_v2, 
       width = 10,
       height = 7,
       units = "in",
       dpi = 600
)

combined_graph_v3 <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
                                  arrangeGrob(
                                    arrangeGrob(arrangeGrob(figs[[1]] + theme(legend.position="none"),
                                                            figs[[3]] + theme(legend.position="none"),
                                                            ncol = 2),
                                                top = textGrob("(A) Parameterization - Pioneer 93B15 at Ambient CO2", 
                                                               gp=gpar(fontface="bold", fontsize=12))),
                                    
                                    arrangeGrob(arrangeGrob(figs[[2]] + theme(axis.title.x = element_blank(),
                                                                              axis.text.x = element_blank(),
                                                                              legend.position="none"),
                                                            figs[[4]] + theme(legend.position="none"),
                                                            nrow = 2, 
                                                            heights = heights, top = 'Ambient CO2'),
                                                
                                                arrangeGrob(figs.elevCO2[[1]] + theme(axis.title.x = element_blank(),
                                                                                      axis.text.x = element_blank(),
                                                                                      legend.position="none"),
                                                            figs.elevCO2[[2]] + theme(legend.position="none"),
                                                            nrow = 2, 
                                                            heights = heights, top = 'Elevated CO2'),
                                                arrangeGrob(figs.elevCO2[[3]] + theme(axis.title.x = element_blank(),
                                                                                      axis.text.x = element_blank(),
                                                                                      legend.position="none"),
                                                            figs.elevCO2[[4]] + theme(legend.position="none"),
                                                            nrow = 2,
                                                            heights = heights, top = 'Elevated CO2'),
                                                ncol = 3, 
                                                top = textGrob("(B) Pioneer 93B15 at SoyFACE", 
                                                               gp=gpar(fontface="bold", fontsize=12))), 
                                    
                                    arrangeGrob(arrangeGrob(ld11.figs[[1]] + theme(axis.title.x = element_blank(),
                                                                                      axis.text.x = element_blank(),
                                                                                      legend.position="none"),
                                                            ld11.figs[[3]] + theme(legend.position="none"),
                                                            nrow = 2, 
                                                            heights = heights),
                                                arrangeGrob(ld11.figs[[2]] + theme(axis.title.x = element_blank(),
                                                                                      axis.text.x = element_blank(),
                                                                                      legend.position="none"),
                                                            ld11.figs[[4]] + theme(legend.position="none"),
                                                            nrow = 2,
                                                            heights = heights),
                                                ncol = 2, 
                                                top = textGrob("(C) LD11-2170 at Ambient CO2", 
                                                               gp=gpar(fontface="bold", fontsize=12))),
                                    
                                    nrow = 3, heights = c(1.3, 2, 2)),
                                  common_legend,
                                  ncol=3, widths=c(0.3, 10, 1.1))

ggsave('Fig2-biomass.png', 
       plot = combined_graph_v3, 
       width = 7.5,
       height = 12,
       units = "in",
       dpi = 600
)

library(dplyr)
total_precip <- sapply(1:4, function(i) {
  results[[i]] %>%
    filter(DVI >= 0 & DVI <= 1) %>%
    summarise(total = sum(precip, na.rm = TRUE)) %>%
    pull(total)
})

names(total_precip) <- years
print(total_precip)


# Plot substrate C concentration
TNC.data <- read.csv('../Data/2022_Carb_data/2022_LD11_TNC_new.csv')
leaf.tnc.all <- read.csv('../Data/2022_Carb_data/Leaf_all.csv')
stem.tnc.all <- read.csv('../Data/2022_Carb_data/Stem_all.csv')

# convert to long format
convert_data_long <- function(tnc.all){
  tnc.long <- data.frame()
  col.names <- colnames(tnc.all)[1:7]
  for (i in 1:nrow(tnc.all)){
    for (j in 1:4){
      if (tnc.all[i, paste0('Valid_',j)]==1){
        new_row <- c(tnc.all[i, col.names], TNC = tnc.all[i, paste0('Carb_tot_',j)])
        tnc.long <- rbind(tnc.long, new_row)
      }
    }
  }
  return (tnc.long)
}
leaf.tnc.long <- convert_data_long(leaf.tnc.all)
stem.tnc.long <- convert_data_long(stem.tnc.all)
# convert unit from nmol glucose / mg (mol glucose / Mg) to mol C / kg
leaf.tnc.long$TNC <- leaf.tnc.long$TNC * 6 / 1000
stem.tnc.long$TNC <- stem.tnc.long$TNC * 6 / 1000

# Simulated
r <- results[[2]]
sim_substrate_C_by_mass <- data.frame(
  time = r$fractional_doy,
  hour = r$hour,
  Leaf = 10 * r$Leaf_substrate_carbon / r$Leaf, # convert from (mol/m2) / (Mg/ha) = 10^4 * mol / Mg = 10 mol / kg
  Stem = 10 * r$Stem_substrate_carbon / r$Stem
)
data_long_per_mass <- tidyr::gather(sim_substrate_C_by_mass[1:2600,
                                                            c('time', 'Leaf', 'Stem')], 
                                    key="Type", value="Value", -time)


# Simulated + Measured
# Leaf
sim_leaf_tnc_by_mass <- sim_substrate_C_by_mass[, c('time','Leaf','hour')]
colnames(sim_leaf_tnc_by_mass)[which(colnames(sim_leaf_tnc_by_mass)=='Leaf')] <- 'TNC'
sim_leaf_tnc_by_mass$Source <- 'Simulated'

leaf.tnc.mean <- leaf.tnc.all[,c('time','Carb_total_mean', 'hour','DOY_date')]
colnames(leaf.tnc.mean)[which(colnames(leaf.tnc.mean)=='Carb_total_mean')] <- 'TNC'
leaf.tnc.mean$Source <- 'Measured (mean)'
# convert unit from nmol glucose / mg (mol glucose / Mg) to mol C / kg
leaf.tnc.mean$TNC <- leaf.tnc.mean$TNC * 6 / 1000

leaf.tnc.long$Source <- 'Measured (individual)'

Leaf.carb.data <- rbind(sim_leaf_tnc_by_mass[,c('time','TNC','Source')], 
                        leaf.tnc.mean[, c('time','TNC','Source')],
                        leaf.tnc.long[, c('time', 'TNC', 'Source')] )


leaf_tnc_plot <- ggplot(Leaf.carb.data, aes(time, TNC, group = Source)) +
  geom_point(aes(shape=Source, color=Source, size=Source, alpha = Source))+
  scale_shape_manual(values=c(8, 18, 16)) +
  scale_size_manual(values=c(1.2, 4, 1)) +
  scale_alpha_manual(values=c(0.6, 1, 0.4)) +
  scale_color_manual(values=c('#D81B60','#117733','grey')) +
  theme_classic() +
  theme(plot.title=element_text(size=size.title, hjust=0.5),
        axis.text=element_text(size=size.axis),
        axis.title=element_text(size=size.axislabel),
        legend.position = c(0.85, 0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(),
       x='Day of Year (2022)',
       y = NULL)
# y='Leaf Substrate C (mol glucose eq./Mg)')

# Stem
sim_stem_tnc_by_mass <- sim_substrate_C_by_mass[, c('time','Stem','hour')]
colnames(sim_stem_tnc_by_mass)[which(colnames(sim_stem_tnc_by_mass)=='Stem')] <- 'TNC'
sim_stem_tnc_by_mass$Source <- 'Simulated'

stem.tnc.mean <- stem.tnc.all[,c('time','Carb_total_mean', 'hour','DOY_date')]
colnames(stem.tnc.mean)[which(colnames(stem.tnc.mean)=='Carb_total_mean')] <- 'TNC'
stem.tnc.mean$Source <- 'Measured (mean)'
# convert unit from nmol glucose / mg (mol glucose / Mg) to mol C / kg
stem.tnc.mean$TNC <- stem.tnc.mean$TNC * 6 / 1000
stem.tnc.long$Source <- 'Measured (individual)'

Stem.carb.data <- rbind(sim_stem_tnc_by_mass[,c('time','TNC','Source')], 
                        stem.tnc.mean[, c('time','TNC','Source')],
                        stem.tnc.long[, c('time', 'TNC', 'Source')] )

stem_tnc_plot <- ggplot(Stem.carb.data, aes(time, TNC, group = Source)) +
  geom_point(aes(shape=Source, color=Source, size=Source, alpha = Source))+
  scale_shape_manual(values=c(8, 18, 16)) +
  scale_size_manual(values=c(1, 4, 1)) +
  scale_alpha_manual(values=c(0.6, 1, 0.4)) +
  scale_color_manual(values=c('#D81B60','#999933','grey')) +
  scale_y_continuous(limits=c(0, 10), n.breaks=6)+
  theme_classic() +
  theme(plot.title=element_text(size=size.title, hjust=0.5),
        axis.text=element_text(size=size.axis),
        axis.title=element_text(size=size.axislabel),
        legend.position = c(0.85, 0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(),
       x='Day of Year (2022)',
       y = NULL)
# y='Stem Substrate C (mol glucose eq./Mg)')

tnc_growingseasion_plot <- grid.arrange(arrangeGrob(textGrob('Substrate C (mol C / kg)', rot = 90, gp=gpar(fontsize=12))),
                                             arrangeGrob(arrangeGrob(leaf_tnc_plot, top = '(A) Leaf'),
                                                         arrangeGrob(stem_tnc_plot, top = '(B) Stem'),
                                                                     ncol = 2, top = 'Pioneer 93B15 (2022) Leaf and Stem Substrate C Concentrations'),
                                             ncol=2, widths=c(0.2, 10))

ggsave('Fig4-substrate-annual.png', 
       plot = tnc_growingseasion_plot, 
       width = 9,
       height = 3.5,
       units = "in",
       dpi = 600
)

# Diurnal changes of substrate C
# take out the last TNC data because in simulation the crop has stopped growing
sim_leaf_tnc_by_mass$DOY <- as.integer(sim_leaf_tnc_by_mass$time)
leaf.tnc.mean$DOY <- as.integer(leaf.tnc.mean$time)
ind.sampling.days <- which(sim_leaf_tnc_by_mass$DOY %in% leaf.tnc.all$DOY)
sim_leaf_tnc_sampling_days <- sim_leaf_tnc_by_mass[ind.sampling.days,]

date.list <- data.frame(
  DOY = c(186, 187, 209, 236, 258),
  DOY_date = c('DOY 186: 07/05', 'DOY 187: 07/06', 'DOY 209: 07/28', 'DOY 236: 08/24', 'DOY 258: 09/15'))

sim_leaf_tnc_sampling_days <- left_join(sim_leaf_tnc_sampling_days, date.list, by = 'DOY')

col.names <- c('DOY_date','hour','TNC', 'Source')
leaf.tnc.sampling.days <- rbind(sim_leaf_tnc_sampling_days[,col.names],
                                leaf.tnc.long[,col.names],
                                leaf.tnc.mean[,col.names])

ggplot(leaf.tnc.sampling.days, aes(hour, TNC, group = Source)) + 
  geom_point(aes(shape=Source, color=Source, size=Source))+
  facet_wrap("DOY_date") +
  scale_shape_manual(values=c(17, 18, 16)) +
  scale_size_manual(values=c(1, 2.5, 0.5)) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  scale_x_continuous(breaks = seq(0,24,6))+
  labs(title=element_blank(), 
       x='Hour of the Day',
       y='Leaf Substrate C (mol glucose eq./Mg)')

sim_stem_tnc_by_mass$DOY <- as.integer(sim_stem_tnc_by_mass$time)
stem.tnc.mean$DOY <- as.integer(stem.tnc.mean$time)
ind.sampling.days <- which(sim_stem_tnc_by_mass$DOY %in% stem.tnc.all$DOY)
sim_stem_tnc_sampling_days <- sim_stem_tnc_by_mass[ind.sampling.days,]

date.list <- data.frame(
  DOY = c(186, 187, 209, 236, 258, 278),
  DOY_date = c('DOY 186: 07/05', 'DOY 187: 07/06', 
               'DOY 209: 07/28', 'DOY 236: 08/24', 
               'DOY 258: 09/15', 'DOY 278: 10/05'))

sim_stem_tnc_sampling_days <- left_join(sim_stem_tnc_sampling_days, date.list, by = 'DOY')

col.names <- c('DOY_date','hour','TNC', 'Source')
stem.tnc.sampling.days <- rbind(sim_stem_tnc_sampling_days[,col.names],
                                stem.tnc.long[,col.names],
                                stem.tnc.mean[,col.names])

ggplot(stem.tnc.sampling.days, aes(hour, TNC, group = Source)) + 
  geom_point(aes(shape=Source, color=Source, size=Source))+
  facet_wrap("DOY_date") +
  scale_shape_manual(values=c(17, 18, 16)) +
  scale_size_manual(values=c(1, 2.5, 0.5)) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  scale_x_continuous(breaks = seq(0,24,6))+
  labs(title=element_blank(), 
       x='Hour',
       y='Stem Substrate C (mol glucose eq./Mg)')

leaf.tnc.sampling.days$Organ <- 'Leaf'
stem.tnc.sampling.days$Organ <- 'Stem'

TNC.sampling.days <- rbind(leaf.tnc.sampling.days, stem.tnc.sampling.days)

TNC.sampling.days <- TNC.sampling.days[-which(TNC.sampling.days$DOY_date=='DOY 278: 10/05'),]
tnc_diurnal_plot <- ggplot(TNC.sampling.days, aes(hour, TNC, group = Source)) + 
  geom_point(data=subset(TNC.sampling.days, Source != 'Simulated'), aes(shape=Source, color=Organ, size=Source, alpha = Source))+
  geom_line(data=subset(TNC.sampling.days, Source == 'Simulated' & Organ == 'Leaf'), aes(color=Organ))+
  geom_line(data=subset(TNC.sampling.days, Source == 'Simulated' & Organ == 'Stem'), aes(color=Organ))+
  facet_wrap("DOY_date") +
  scale_shape_manual(values=c(8, 18)) +
  scale_size_manual(values=c(2, 4)) +
  scale_color_manual(values = c('#117733', '#999933'))+
  scale_alpha_manual(values=c(0.5, 1)) +
  theme_classic() +
  theme(legend.position = c(0.84, 0.18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  scale_x_continuous(breaks = seq(0,24,6))+
  labs(title=element_blank(), 
       x='Hour',
       y='Substrate C (mol C / kg)')
# y='Substrate C (mol glucose eq./Mg)')
ggsave('Fig5-substrate-diurnal.png', 
       plot = tnc_diurnal_plot, 
       width = 9,
       height = 3.5,
       units = "in",
       dpi = 600
)
