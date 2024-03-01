# The work directory is hard coded since the previous one doesn't work well all the time. 
# Please change it according to your own working directory.
# Parameter loading and plotting are put into this R code to check the variables more easily.
# To be moved to different files once everything works out.
# Clear workspace
rm(list=ls())

library(BioCro)
library(UTRSoybeanBML)
library(DEoptim)
library(ggplot2)
library(lattice)

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Cost function
source('../ParameterOptimization/soybean_parameter_expansion.R')
source('../ParameterOptimization/multiyear_BioCro_optim_lsrg.R')

# set year and CO2 level
co2_opt = '_ambient_' # '_ambient_' or '_co2_'

# years, sowing dates, and harvesting dates of growing seasons being fit to
year <- c('2002', '2004', '2005', '2006')
sow.date <- c(152, 149, 148, 148)
harv.date <- c(288, 289, 270, 270)
weather.growingseason <- list()

# names of fitted parameters
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

# Initialize variables for the cost function
soybean_optsolver <- list()
ExpBiomass <- list()
ExpBiomass.std <- list()
RootVals <- list()
weights <- list()
numrows <- vector()

setwd('../Data/Soybean-BioCro_Parameters')
# load parameter files
param_files <- list.files(pattern = "[.]R$", recursive = TRUE)
param_files <- param_files[-c(which(param_files == "soybean_initial_values.R"))]
lapply(param_files, source)

for (i in 1:length(year)) { 
  yr <- year[i]
  weather <- read.csv(file = paste0('../Weather_data/', yr,'_Bondville_IL_daylength_wDVI.csv'))
  emergence.idx <- which(weather$DVI>0)[1]
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather.growingseason[[i]] <- weather[emergence.idx:hd.ind,]
  
  ExpBiomass[[i]] <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass.csv'))
  colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  ExpBiomass.std[[i]] <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  # RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[3], 
  #                             "Root"=0.17*sum(ExpBiomass[[i]][5,2:4])) 
  RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[-c(1,2)], 
                              "Root"=0.17*sum(ExpBiomass[[i]][5,2:4])) 
  # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  numrows[i] <- nrow(weather.growingseason[[i]])
  invwts <- ExpBiomass.std[[i]]
  weights[[i]] <- log(1/(invwts[,2:ncol(invwts)]+1e-5))
  
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
}

# Optimization
wts2 <- data.frame("Stem" = 1, "Leaf" = 1, "Pod" = 1, "Root" = 0.5, "CumLitter" = 0.5)

# cost function
cost_func <- function(x){
  multiyear_BioCro_optim(x, soybean_optsolver[c(1,3)], ExpBiomass[c(1,3)], 
                         numrows[c(1,3)], weights[c(1,3)], wts2, RootVals[c(1,3)])
}
# calculate the cost for references
optim_params_short <-c(0.244284,    0.022706,    0.040001,    0.037560,    
                       0.552709,    1.316245,    1.478678,    1.431808,    
                       0.908055,    0.081352,    0.685889,    0.513119,    
                       1.887088,    0.069919,    0.099003,    0.012061,   
                       17.765036,   24.410462,    6.975161,   
                       -9.151024,   -9.928725,   -2.014179)
optim_params_short <-c(0.24428360,  0.02270617,  0.04000102,  0.03756000,  0.55270866,  1.31624494,  1.47867831,
                       1.43180810,  0.90805499,  0.08135216,  0.68588900,  0.51311881,  1.88708826,  0.06991942,
                       0.09900254,  0.01206099, 17.76503630, 24.41046211,  6.97516144, -9.15102399, -9.92872499,
                       -2.01417925)
# testing
for (i in 1:length(year)) {
  print(paste0('Year ', year[i]))
  result <- match.fun(soybean_optsolver[[i]])(optim_params_conversion(optim_params_short))
  print(paste0('1st nan value occurrs at time: ', result$time[is.nan(result$Leaf)][1]))
  print(paste0('1st nan value occurrs at DVI: ', result$DVI[is.nan(result$Leaf)][1]))
  print(paste0('nrow(result):', nrow(result), ' nrow(weather):', nrow(weather.growingseason[[i]])))
  xyplot(data = result, Leaf+Stem+Root+Pod~time, auto.key = TRUE)
}

