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
LAI <- list()
LAI.elevCO2 <- list()
ExpBiomass.std <- list()
RootVals <- list()
numrows <- vector()
figs <- list()
lai.figs <- list()
pod_sensitivity <- data.frame()

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
  s.lai <- cbind(res[,c("time","lai_tot")],elev_res[,"lai_tot"])
  colnames(s.lai) <- c("time","Amb","Elev")
  r.lai <- reshape2::melt(s.lai, id.vars = "time")
  
  s.exp.lai <- cbind(lai[,c("DOY","LAI_mean")],elev_lai[,"LAI_mean"])
  colnames(s.exp.lai) <- c("DOY","Amb","Elev")
  r.exp.lai <- reshape2::melt(s.exp.lai, id.vars = "DOY")
  
  s.exp.std.lai <- cbind(lai[,c("DOY","LAI_std")],elev_lai[,"LAI_std"])
  colnames(s.exp.std.lai) <- c("DOY","Amb","Elev")
  r.exp.std.lai <- reshape2::melt(s.exp.std.lai, id.vars = "DOY")
  r.exp.std.lai$ymin <- r.exp.lai$value - r.exp.std.lai$value
  r.exp.std.lai$ymax <- r.exp.lai$value + r.exp.std.lai$value
  
  f <- ggplot() + theme_classic()
  f <- f + geom_point(data=r.lai, aes(x=time, y=value, colour=variable),show.legend = FALSE,size=0.25)
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
  
  # Load LAI data
  if(i>1){
    LAI[[i]]<-read.csv(file=paste0('../SoyFACE_data/',yr,'_ambient_lai.csv'))
    LAI.elevCO2[[i]]<-read.csv(file=paste0('../SoyFACE_data/',yr,'_elevated_lai.csv'))
  }
  
  # load parameter files
  source('soybean_initial_values.R')
  initial_state$DVI <- NULL
  # parameters$StomataWS <- 1
  soybean_optsolver[[i]] <- partial_run_biocro(initial_state,
                                               parameters,
                                               weather.growingseason[[i]],
                                               steady_state_module_names, # [-1],
                                               derivative_module_names,
                                               solver,
                                               arg_names,
                                               verbose = FALSE)
  
  result <- soybean_optsolver[[i]](optim_params_conversion(optim_params_short_SoyFACE)) 
  results[[i]] <- result
  
  # organize simulated data
  r.lsrp.doy <- reshape2::melt(result[,c("time","Root","Leaf","Stem","Pod")],id.vars="time")
  
  # Leaf
  # organize the experimental data (mean and std)
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
                aes(x=time,y=value, color=Organ), size=0.8, alpha = 0.8) +
    geom_point(data = subset(r.all, Source == "Observed"), 
                 aes(x=time, y=value, color=Organ), shape=15, size=2, stroke=.5) +
  
    scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, 2)) +
    scale_color_manual(values = col.palette.muted.organs)

  
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
  f <- f + labs(title=years[i], x='Day of Year', y=NULL)
  if (years[i] == '2002' || years[i] == '2004'){
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
  # parameters$StomataWS <- 1
  soybean_optsolver[[i]] <- partial_run_biocro(initial_state,
                                               parameters,
                                               weather.growingseason[[i]],
                                               steady_state_module_names,#[-1],
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
              aes(x=time,y=value, color=Organ), size=0.8, alpha = 0.8) +
    geom_point(data = subset(r.all, Source == "Observed"), 
               aes(x=time, y=value, color=Organ), shape=15, size=2, stroke=.5) +
    
    scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, 2)) +
    scale_color_manual(values = col.palette.muted.organs)
  
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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
save(results, results.elevCO2, ExpBiomass, ExpBiomass.elevCO2, file = 'SoyFACE_results_and_measurements.RData')

source('plot_partitioning.R')

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
    time = results[[i]]$time,
    ambient_ratio = with(results[[i]], (Leaf+Stem+Pod)/Root),
    elevated_ratio = with(results.elevCO2[[i]], (Leaf+Stem+Pod)/Root))
  print(years[i])
  print(xyplot(data=shoot_root_ratio, ambient_ratio + elevated_ratio ~ time, auto.key = TRUE))
}
