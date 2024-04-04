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
co2_opt = '_co2_' # '_ambient_' or '_co2_'

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
optim_params_short <-c(0.3,    0.01,    0.02,    0.02,     0.8, # 1-5  # carbon to mass factor and utilization rate constants
                       0.2,    0.1,    0.05,    0.1,    0.5,     # 6-10 Km s and respiration factor
                       0.4,    0.2,    5,  # 11-13  substrate conductance
                       0.01,     0.01,    0.002, # 14-16    senescence max rates [/hr]
                       5.0,     8.0,     1.0,   # 17-22 alphas and betas
                       1.7,     1.9,     1.6)
optim_params_short <-c(0.395667,    0.004489,    0.018746,    0.090290,    0.588019,    
                       0.419076,    0.177943,    0.466320,    0.211607,    0.200006,    
                       0.025633,    0.408561,    1.105566,    
                       0.041240,    0.002908,    0.000173,    
                       8.783273,    9.921733,    1.711259,    
                       1.872338,    1.992031,    1.665296)

# testing
for (i in 1:length(year)) {
  print(paste0('Year ', year[i]))
  result <- match.fun(soybean_optsolver[[i]])(optim_params_conversion(optim_params_short))
  print(paste0('1st nan value occurrs at time: ', result$time[is.nan(result$Leaf)][1]))
  print(paste0('1st nan value occurrs at DVI: ', result$DVI[is.nan(result$Leaf)][1]))
  print(paste0('nrow(result):', nrow(result), ' nrow(weather):', nrow(weather.growingseason[[i]])))
  print(xyplot(data = result, Leaf+Stem+Root+Pod~time, auto.key = TRUE, main = year[i]))
  Pod_start_idx <- which.min(abs(result$DVI-1))+1
  print(paste0("The first transport to Pod was: ", result$substrate_transport_Stem_to_Pod[Pod_start_idx]))
  print(paste0("The Stem Substrate C concentration was: ", 
               result$Stem_substrate_carbon[Pod_start_idx]/
                 result$Stem_structural_carbon[Pod_start_idx]))
  
  print(paste0("The Pod Substrate C concentration was: ", 
               result$Pod_substrate_carbon[Pod_start_idx]/
                 result$Pod_structural_carbon[Pod_start_idx]))
}

result <- match.fun(soybean_optsolver[[1]])(optim_params_conversion(optim_params_short))
xyplot(data=result, 
       substrate_transport_Stem_to_Pod
       ~time, auto=TRUE)

