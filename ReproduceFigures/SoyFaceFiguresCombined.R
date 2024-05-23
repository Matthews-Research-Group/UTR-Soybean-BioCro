library(BioCro)
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

# years, sowing dates, and harvesting dates of growing seasons being fit to
years <- c('2002', '2004', '2005', '2006')
sow.date <- c(152, 149, 148, 148)
harv.date <- c(288, 289, 270, 270)

co2_opt = '_ambient_'
# load parameter files
setwd('../Data/Soybean-BioCro_Parameters')

param_files <- list.files(pattern = "[.]R$", recursive = TRUE)

param_files <- param_files[-c(which(param_files == "soybean_initial_values.R"))]
lapply(param_files, source)

# Initialize lists
results <- list()
results.elevCO2 <- list()
weather.growingseason <- list()
soybean_optsolver <- list()
ExpBiomass <- list()
ExpBiomass.elevCO2 <- list()
ExpBiomass.std <- list()
RootVals <- list()
numrows <- vector()
figs <- list()

for (i in 1:length(years)) {   
  yr <- years[i]
  weather <- read.csv(file = paste0('../Weather_data/', yr,'_Bondville_IL_daylength_wDVI.csv'))
  # sd.ind <- which(weather$doy == sow.date[i])[1]
  em.ind <- which(weather$DVI > 0)[1] # emergence doy
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather.growingseason[[i]] <- weather[em.ind:hd.ind,]
  
  ExpBiomass[[i]] <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass.csv'))
  colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  ExpBiomass.std[[i]] <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[3], "Root"=0.17*sum(ExpBiomass[[i]][5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  numrows[i] <- nrow(weather.growingseason[[i]])
  invwts <- ExpBiomass.std[[i]]
  # load parameter files
  source('soybean_initial_values.R')
  initial_state$DVI <- NULL
  soybean_optsolver[[i]] <- partial_run_biocro(initial_state,
                                               parameters,
                                               weather.growingseason[[i]],
                                               steady_state_module_names,
                                               derivative_module_names,
                                               solver,
                                               arg_names,
                                               verbose = FALSE)
  
  result <- soybean_optsolver[[i]](optim_params_conversion(optim_params_short_SoyFACE)) 
  results[[i]] <- result
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
  r.lsrp.doy <- reshape2::melt(result[,c("time","Root","Leaf","Stem","Pod")],id.vars="time")
  r.lsrp.doy$value<-r.lsrp.doy$value
  
  # Leaf
  # organize the expirimental data (mean and std)
  s.exp.leaf <- cbind(ExpBiomass[[i]][,c("DOY","Leaf")])
  colnames(s.exp.leaf) <- c("time","Leaf") # DOY renamed as time
  r.exp.leaf <- reshape2::melt(s.exp.leaf, id.vars = "time") # DOY renamed as time
  
  s.exp.std.leaf <- cbind(ExpBiomass.std[[i]][,c("DOY","Leaf")])
  colnames(s.exp.std.leaf) <- c("time","Leaf") # DOY renamed as time
  r.exp.std.leaf <- reshape2::melt(s.exp.std.leaf, id.vars = "time") # DOY renamed as time
  r.exp.std.leaf$ymin <- r.exp.leaf$value - r.exp.std.leaf$value
  r.exp.std.leaf$ymax <- r.exp.leaf$value + r.exp.std.leaf$value
  
  # Stem
  s.exp.stem <- cbind(ExpBiomass[[i]][,c("DOY","Stem")]) # DOY renamed as time
  colnames(s.exp.stem) <- c("time","Stem") # DOY renamed as time
  r.exp.stem <- reshape2::melt(s.exp.stem, id.vars = "time") # DOY renamed as time
  
  s.exp.std.stem <- cbind(ExpBiomass.std[[i]][,c("DOY","Stem")]) # DOY renamed as time
  colnames(s.exp.std.stem) <- c("time","Stem") # DOY renamed as time
  r.exp.std.stem <- reshape2::melt(s.exp.std.stem, id.vars = "time") # DOY renamed as time
  r.exp.std.stem$ymin <- r.exp.stem$value - r.exp.std.stem$value
  r.exp.std.stem$ymax <- r.exp.stem$value + r.exp.std.stem$value
  
  # Pod
  s.exp.pod <- cbind(ExpBiomass[[i]][,c("DOY","Pod")]) # DOY renamed as time
  colnames(s.exp.pod) <- c("time","Pod") # DOY renamed as time
  r.exp.pod <- reshape2::melt(s.exp.pod, id.vars = "time") # DOY renamed as time
  
  s.exp.std.pod <- cbind(ExpBiomass.std[[i]][,c("DOY","Pod")]) # DOY renamed as time
  colnames(s.exp.std.pod) <- c("time","Pod") # DOY renamed as time
  r.exp.std.pod <- reshape2::melt(s.exp.std.pod, id.vars = "time") # DOY renamed as time
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
  col.palette.muted <- c("#332288", "#999933", "#117733", "#882255")
  
  size.title <- 12
  size.axislabel <- 10
  size.axis <- 10
  size.legend <- 12
  
  f <- ggplot() + theme_classic()
  
  f <- f + geom_point(data=r.all, aes(x=time, y=value,
                                      color=Organ,
                                      size = Source, shape = Source),
                      show.legend = TRUE, stroke=0.5) +
    scale_shape_manual(values = c(15, 16)) +
    scale_size_manual(values = c(2, 0.25)) +
    scale_color_manual(values = col.palette.muted)
  
  # for leaf
  f <- f + geom_errorbar(data=r.exp.std.leaf, aes(x=time, ymin=ymin, ymax=ymax),  # DOY renamed as time
                         width=3.5, size=0.25, show.legend = FALSE)
  # for stem
  f <- f + geom_errorbar(data=r.exp.std.stem, aes(x=time, ymin=ymin, ymax=ymax),   # DOY renamed as time
                         width=3.5, size=0.25, show.legend = FALSE)
  
  # for pod
  f <- f + geom_errorbar(data=r.exp.std.pod, aes(x=time, ymin=ymin, ymax=ymax),   # DOY renamed as time
                         width=3.5, size=0.25, show.legend = FALSE)
  
  # change the plot labels and theme
  print(years[i])
  f <- f + labs(title=element_blank(), x=paste0('Day of Year (',years[i],')'),y=NULL)
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
}

canopy_assim_daily <- aggregate(result$canopy_assimilation_rate,list(result$doy), FUN=sum) * 10 / 3
leaf_export_daily <- aggregate(result$substrate_transport_Leaf_to_Stem,list(result$doy), FUN=sum)
leaf_utilization_daily <- aggregate(result$Leaf_utilization_rate,list(result$doy), FUN=sum)
stem_utilization_daily <- aggregate(result$Stem_utilization_rate,list(result$doy), FUN=sum)
df <- data.frame(DOY = leaf_export_daily$Group.1, 
                 Canopy_Assimilation_mol_per_day = canopy_assim_daily$x,
                 Leaf_Export_mol_per_day = leaf_export_daily$x,
                 Leaf_Utilization_mol_per_day = leaf_utilization_daily$x,
                 Stem_Utilization_mol_per_day = stem_utilization_daily$x)
xyplot(data = df,
       Canopy_Assimilation_mol_per_day+
         Leaf_Export_mol_per_day+
         Leaf_Utilization_mol_per_day+
         Stem_Utilization_mol_per_day~
         DOY, type = c('p','l'), auto=TRUE)
sum(result$substrate_transport_Leaf_to_Stem)
sum(result$Stem_utilization_rate)

# initialize lists for figures
figs.elevCO2 <- list()
co2_opt = '_co2_'

# Initialize lists
weather.growingseason <- list()
soybean_optsolver <- list()
ExpBiomass.std <- list()
RootVals <- list()
numrows <- vector()

# load parameter files
param_files <- list.files(pattern = "[.]R$", recursive = TRUE)

param_files <- param_files[-c(which(param_files == "soybean_initial_values.R"))]
lapply(param_files, source)

for (i in 1:length(years)) {  
# for (i in c(1,3)) { 
  yr <- years[i]
  weather <- read.csv(file = paste0('../Weather_data/', yr,'_Bondville_IL_daylength_wDVI.csv'))
  # sd.ind <- which(weather$doy == sow.date[i])[1]
  em.ind <- which(weather$DVI > 0)[1] # emergence doy
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather.growingseason[[i]] <- weather[em.ind:hd.ind,]
  
  ExpBiomass.elevCO2[[i]] <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass.csv'))
  colnames(ExpBiomass.elevCO2[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  ExpBiomass.std[[i]] <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  RootVals[[i]] <- data.frame("DOY"=ExpBiomass.elevCO2[[i]]$DOY[3], "Root"=0.17*sum(ExpBiomass.elevCO2[[i]][5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  numrows[i] <- nrow(weather.growingseason[[i]])
  invwts <- ExpBiomass.std[[i]]
  # load parameter files
  source('soybean_initial_values.R')
  initial_state$DVI <- NULL
  soybean_optsolver[[i]] <- partial_run_biocro(initial_state,
                                               parameters,
                                               weather.growingseason[[i]],
                                               steady_state_module_names,
                                               derivative_module_names,
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
  r.lsrp.doy <- reshape2::melt(result[,c("time","Root","Leaf","Stem","Pod")],id.vars="time")
  r.lsrp.doy$value<-r.lsrp.doy$value
  
  # Leaf
  # organize the expirimental data (mean and std)
  s.exp.leaf <- cbind(ExpBiomass.elevCO2[[i]][,c("DOY","Leaf")])
  colnames(s.exp.leaf) <- c("time","Leaf") # DOY renamed as time
  r.exp.leaf <- reshape2::melt(s.exp.leaf, id.vars = "time") # DOY renamed as time
  
  s.exp.std.leaf <- cbind(ExpBiomass.std[[i]][,c("DOY","Leaf")])
  colnames(s.exp.std.leaf) <- c("time","Leaf") # DOY renamed as time
  r.exp.std.leaf <- reshape2::melt(s.exp.std.leaf, id.vars = "time") # DOY renamed as time
  r.exp.std.leaf$ymin <- r.exp.leaf$value - r.exp.std.leaf$value
  r.exp.std.leaf$ymax <- r.exp.leaf$value + r.exp.std.leaf$value
  
  # Stem
  s.exp.stem <- cbind(ExpBiomass.elevCO2[[i]][,c("DOY","Stem")]) # DOY renamed as time
  colnames(s.exp.stem) <- c("time","Stem") # DOY renamed as time
  r.exp.stem <- reshape2::melt(s.exp.stem, id.vars = "time") # DOY renamed as time
  
  s.exp.std.stem <- cbind(ExpBiomass.std[[i]][,c("DOY","Stem")]) # DOY renamed as time
  colnames(s.exp.std.stem) <- c("time","Stem") # DOY renamed as time
  r.exp.std.stem <- reshape2::melt(s.exp.std.stem, id.vars = "time") # DOY renamed as time
  r.exp.std.stem$ymin <- r.exp.stem$value - r.exp.std.stem$value
  r.exp.std.stem$ymax <- r.exp.stem$value + r.exp.std.stem$value
  
  # Pod
  s.exp.pod <- cbind(ExpBiomass.elevCO2[[i]][,c("DOY","Pod")]) # DOY renamed as time
  colnames(s.exp.pod) <- c("time","Pod") # DOY renamed as time
  r.exp.pod <- reshape2::melt(s.exp.pod, id.vars = "time") # DOY renamed as time
  
  s.exp.std.pod <- cbind(ExpBiomass.std[[i]][,c("DOY","Pod")]) # DOY renamed as time
  colnames(s.exp.std.pod) <- c("time","Pod") # DOY renamed as time
  r.exp.std.pod <- reshape2::melt(s.exp.std.pod, id.vars = "time") # DOY renamed as time
  r.exp.std.pod$ymin <- r.exp.pod$value - r.exp.std.pod$value
  r.exp.std.pod$ymax <- r.exp.pod$value + r.exp.std.pod$value
  
  # Combine
  r.exp.ls <- rbind(r.exp.leaf, r.exp.stem, r.exp.pod)
  r.exp.ls$Source = "Observed"
  r.lsrp.doy$Source = "Simulated"
  r.all <- rbind(r.lsrp.doy, r.exp.ls)
  # Reverse the order as follow
  r.all$Organ <- factor(r.all$variable, levels = rev(levels(r.all$variable)))
  
  # combine the simulated and experimental data
  f <- ggplot() + theme_classic()
  f <- f + geom_point(data=r.all, aes(x=time, y=value,
                                      color=Organ,
                                      size = Source, shape = Source),
                      show.legend = TRUE, stroke=0.5) +
    scale_shape_manual(values = c(15, 16)) +
    scale_size_manual(values = c(2, 0.25)) +
    scale_color_manual(values = col.palette.muted)
  # for leaf
  f <- f + geom_errorbar(data=r.exp.std.leaf, aes(x=time, ymin=ymin, ymax=ymax),  # DOY renamed as time
                         width=3.5, size=0.25, show.legend = FALSE)
  # for stem
  f <- f + geom_errorbar(data=r.exp.std.stem, aes(x=time, ymin=ymin, ymax=ymax),   # DOY renamed as time
                         width=3.5, size=0.25, show.legend = FALSE)
  # for pod
  f <- f + geom_errorbar(data=r.exp.std.pod, aes(x=time, ymin=ymin, ymax=ymax),   # DOY renamed as time
                         width=3.5, size=0.25, show.legend = FALSE)
  
  # change the plot labels and theme
  print(years[i])
  f <- f + labs(title=element_blank(), x=paste0('Day of Year (',years[i],')'),y=NULL)
  f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
                 axis.text=element_text(size=size.axis),
                 axis.title=element_text(size=size.axislabel),
                 # legend.position = c(.15,.8), legend.title = element_blank(),
                 # legend.text=element_text(size=size.legend),
                 # legend.background = element_rect(fill = "transparent",colour = NA),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
                 plot.background = element_rect(fill = "transparent", colour = NA))
  f <- f + scale_x_continuous(breaks = seq(150,280,30))
  figs.elevCO2[[i]] <- f
}

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

common_legend <- g_legend(figs[[1]])

combined_graph <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
                               arrangeGrob(arrangeGrob(figs[[1]] + theme(legend.position="none"),
                                                        figs[[3]] + theme(legend.position="none"),
                                                        nrow = 2, top = 'Training'),
                                            arrangeGrob(figs[[2]] + theme(legend.position="none"),
                                                        figs[[4]] + theme(legend.position="none"),
                                                        nrow = 2, top = 'Testing'),
                                            ncol = 2, top = 'Ambient CO2'),
                                arrangeGrob(arrangeGrob(figs.elevCO2[[1]] + theme(legend.position="none"),
                                                        figs.elevCO2[[3]] + theme(legend.position="none"),
                                                        nrow = 2, top = 'Testing'),
                                            arrangeGrob(figs.elevCO2[[2]] + theme(legend.position="none"),
                                                        figs.elevCO2[[4]] + theme(legend.position="none"),
                                                        nrow = 2, top = 'Testing'),
                                            ncol = 2, top = 'Elevated CO2'),
                                common_legend,
                                ncol=4, widths=c(0.3, 5,5,1.2))


# day_result <- result[1210:1233,c('hour',
#                                  'Leaf_substrate_carbon',
#                                  'Stem_substrate_carbon', 
#                                  'Root_substrate_carbon')]
# colnames(day_result) = c('hour', 'Leaf', 'Stem', 'Root')
# 
# day_result <- reshape2::melt(day_result, 
#                             id.vars = "hour")
# colnames(day_result) = c('Hour', 'Organ', 'Value')
# ggplot(data = day_result, aes(x=Hour, y=Value)) + 
#   theme_classic() +
#   ylab('Substrate C (mol / m^2 / hr)') +
#   geom_point(aes(col=Organ)) +
#   theme(legend.position="top")

# canopy_assim_daily <- aggregate(result$canopy_assimilation_rate,list(result$doy), FUN=sum) * 10 / 3
# leaf_export_daily <- aggregate(result$substrate_transport_Leaf_to_Stem,list(result$doy), FUN=sum)
# leaf_utilization_daily <- aggregate(result$Leaf_utilization_rate,list(result$doy), FUN=sum)
# stem_utilization_daily <- aggregate(result$Stem_utilization_rate,list(result$doy), FUN=sum)
# df <- data.frame(DOY = leaf_export_daily$Group.1, 
#                  Canopy_Assimilation_mol_per_day = canopy_assim_daily$x,
#                  Leaf_Export_mol_per_day = leaf_export_daily$x,
#                  Leaf_Utilization_mol_per_day = leaf_utilization_daily$x,
#                  Stem_Utilization_mol_per_day = stem_utilization_daily$x)
# xyplot(data = df,
#        Canopy_Assimilation_mol_per_day+
#          Leaf_Export_mol_per_day+
#          Leaf_Utilization_mol_per_day+
#          Stem_Utilization_mol_per_day~
#          DOY, type = c('p','l'), auto=TRUE)
# sum(result$substrate_transport_Leaf_to_Stem)
# sum(result$Stem_utilization_rate)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
save(results, results.elevCO2, ExpBiomass, ExpBiomass.elevCO2, file = 'SoyFACE_results_and_measurements.RData')

i <- 2
sDVI <- parameters$stop_growth_dvi
xyplot(results.elevCO2[[i]]$Leaf_substrate_carbon[which(results.elevCO2[[i]]$DVI<sDVI)]+
         results[[i]]$Leaf_substrate_carbon[which(results[[1]]$DVI<sDVI)]
       ~results[[i]]$time[which(results.elevCO2[[i]]$DVI<sDVI)])

for (i in 1:4){
  print(xyplot(data = results[[i]][1:which.min(abs(results[[i]]$DVI-2)),], 
         Leaf_substrate_carbon/Leaf_structural_carbon+
           Stem_substrate_carbon/Stem_structural_carbon~
           time,
         auto=TRUE,
         ylim = c(0,1),
         main = years[i]))
}
