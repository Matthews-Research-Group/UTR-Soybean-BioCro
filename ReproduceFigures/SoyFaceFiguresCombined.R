library(BioCro)
library(ggplot2)
library(grid)
library(gridExtra)
library(lattice)

# Clear workspace
rm(list=ls())

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Cost function
source('../ParameterOptimization/soybean_parameter_expansion.R')

# years, sowing dates, and harvesting dates of growing seasons being fit to
years <- c('2002', '2004', '2005', '2006')
sow.date <- c(152, 149, 148, 148)
harv.date <- c(288, 289, 270, 270)

# initialize lists for figures
figs <- list()
co2_opt = '_ambient_'

# names of fitted parameters
arg_names <- c('Leaf_carbon_to_mass_factor', 'Stem_carbon_to_mass_factor', # 1, 2 
               'Root_carbon_to_mass_factor', 'Pod_carbon_to_mass_factor',  # 3, 4
               'Leaf_utilization_rate_constant', 'Stem_utilization_rate_constant', # 5, 6
               'Root_utilization_rate_constant', 'Pod_utilization_rate_constant',  # 7, 8
               'Leaf_utilization_km', 'Stem_utilization_km', 'Root_utilization_km', 'Pod_utilization_km', # 9,10,11,12
               'Leaf_respiration_factor', 'Stem_respiration_factor', # 13, 14
               'Root_respiration_factor', 'Pod_respiration_factor', # 15, 16
               'substrate_conductance_Leaf_to_Stem', 'substrate_conductance_Stem_to_Root', # 17, 18
               'substrate_conductance_Stem_to_Pod', # 'transportation_beta_exponent', # 19, 20
               'Leaf_senescence_rate_max','Stem_senescence_rate_max', # 21, 22
               'Leaf_senescence_alpha', 'Stem_senescence_alpha',# 23, 24
               'Leaf_senescence_beta', 'Stem_senescence_beta', # 25, 26
               'Pod_start_dvi', 'stop_growth_dvi') # 27, 28
arg_names <- c('Leaf_carbon_to_mass_factor', 'Stem_carbon_to_mass_factor', # 1, 2 
               'Root_carbon_to_mass_factor', 'Pod_carbon_to_mass_factor',  # 3, 4
               'Leaf_utilization_rate_constant', 'Stem_utilization_rate_constant', # 5, 6
               'Root_utilization_rate_constant', 'Pod_utilization_rate_constant',  # 7, 8
               'Leaf_utilization_km', 'Stem_utilization_km', 'Root_utilization_km', 'Pod_utilization_km', # 9,10,11,12
               'Leaf_respiration_factor', 'Stem_respiration_factor', # 13, 14
               'Root_respiration_factor', 'Pod_respiration_factor', # 15, 16
               'substrate_conductance_Leaf_to_Stem', 'substrate_conductance_Stem_to_Root', # 17, 18
               'substrate_conductance_Stem_to_Pod', #  'transportation_beta_exponent', # 19
               'Leaf_senescence_rate_max','Stem_senescence_rate_max', 'Root_senescence_rate_max', # 20, 21，22
               'Leaf_senescence_alpha', 'Stem_senescence_alpha', 'Root_senescence_alpha',# 23，24, 25,
               'Leaf_senescence_beta', 'Stem_senescence_beta', 'Root_senescence_beta')# 26，27, 28

# optim_params_short <-c(0.242485,    0.004276,    0.032644,    0.053587,    
#                        0.586933,    0.062188,    0.902778,    0.832786,    
#                        0.485468,    0.357419,    0.352924,    0.373693,   
#                        1.157740,    0.858463,    0.008441,    0.020550,   #  0.858463
#                        16.118925,   20.463494,   
#                        -9.096414,   -9.807364,    0.825496,    1.978349) #0.825496
# 
# optim_params_short <-c(0.233784,    0.004728,    0.034942,    0.022596,
#                       0.962753,    0.164531,    1.840606,    1.357177,
#                       3.092379,    0.293349,    0.317670,    0.779145,
#                       1.198915,    0.013672,    0.039876,
#                       10.891061,   14.012313,   -5.537149,   -5.798883,    
#                       0.883654,    1.995483)
optim_params_short <-c(0.280385,    0.013297,    0.025664,    0.030763,    
                       0.264599,    0.784713,    1.045214,    1.492283,    
                       0.559731,    0.124240,    0.960166,    0.935573,    
                       1.138520,    0.021318,    0.010909,    9.443322,    
                       8.157651,   -4.506269,   -3.057007,    1,    1.965775) #0.982365
optim_params_short <-c(0.273470,    0.028917,    0.067696,    0.048373,    
                       0.715883,    1.511239,    1.799835,    1.165064,    
                       0.564379,    0.130939,    0.548013,    0.746798,    
                       1.673044,    0.049328,    0.002033,   16.316789,   
                       17.809283,   -8.156955,   -9.435581,    0.993601,    1.895145)
optim_params_short <-c(0.261472,    0.023359,    0.031246,    0.065703,    
                       0.339670,    1.480430,    0.876905,    1.537500,    
                       0.119764,    0.155933,    0.536036,    0.411301,    
                       1.002809,    
                       0.080179,    0.026975,    0.036658,   
                       14.857419,   14.663911,    9.776203,   
                       -6.546716,   -5.308754,   -3.626154)
optim_params_short <-c(0.321679,    0.019472,    0.019775,    0.055453,    
                       0.460345,    1.437916,    0.665640,    1.092556,    
                       1.535427,    0.195530,    0.678282,    0.344428,    
                       1.424100,    0.098971,    0.098040,    0.072484,   
                       14.963209,   13.831948,    3.832358,   
                       -7.022856,   -4.641297,   -0.657350)

optim_params_short <-c(0.24428360,  0.02270617,  0.04000102,  0.03756000,  0.55270866,  1.31624494,  1.47867831,
                       1.43180810,  0.90805499,  0.08135216,  0.68588900,  0.51311881,  1.88708826,  0.06991942,
                       0.09900254,  0.01206099, 17.76503630, 24.41046211,  6.97516144, -9.15102399, -9.92872499,
                       -2.01417925)
# calculate the cost for references
optim_params_short <-c(0.244284,    0.022706,    0.040001,    0.037560,    
                       0.552709,    1.316245,    1.478678,    1.431808,    
                       0.908055,    0.081352,    0.685889,    0.513119,    
                       1.887088,    0.069919,    0.099003,    0.012061,   
                       17.765036,   24.410462,    6.975161,   
                       -9.151024,   -9.928725,   -2.014179)

setwd('../Data/Soybean-BioCro_Parameters')

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

# load parameter files
param_files <- list.files(pattern = "[.]R$", recursive = TRUE)

param_files <- param_files[-c(which(param_files == "soybean_initial_values.R"))]
lapply(param_files, source)

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
  soybean_optsolver[[i]] <- partial_run_biocro(initial_state,
                                               parameters,
                                               weather.growingseason[[i]],
                                               steady_state_module_names,
                                               derivative_module_names,
                                               solver,
                                               arg_names,
                                               verbose = FALSE)
  
  result <- soybean_optsolver[[i]](optim_params_conversion(optim_params_short))
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
  soybean_optsolver[[i]] <- partial_run_biocro(initial_state,
                                               parameters,
                                               weather.growingseason[[i]],
                                               steady_state_module_names,
                                               derivative_module_names,
                                               solver,
                                               arg_names,
                                               verbose = FALSE)
  
  result <- soybean_optsolver[[i]](optim_params_conversion(optim_params_short))
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
