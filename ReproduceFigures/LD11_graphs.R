# Clear the workspace
rm(list=ls())

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source files
source('../ParameterOptimization/soybean_parameter_expansion.R')
co2_opt = '_tbd_'
years <- c('2021', '2022')
Catms <- c(414.7, 417.2) # from NOAA
source('../Data/Soybean-BioCro_Parameters/soybean_parameters.R')
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')
source('../Data/Soybean-BioCro_Parameters/soybean_modules.R')
source('plot_partitioning.R')
# Load packages
library(BioCro)
library(lattice)
library(reshape2)
library(ggplot2)

# Load shared data
load('../Data/Soybean-BioCro_Parameters/full_soybean_ld11.RData')

# Modify the modules
full_soybean_ld11$direct_modules <- steady_state_module_names
full_soybean_ld11$differential_modules <- derivative_module_names

# Set up basic properties for the solver
full_soybean_ld11$ode_solver <- solver

# Update UTR parameters
fitted.utr.params <- optim_params_conversion(optim_params_short)
names(fitted.utr.params) <- arg_names
parameters <-c(parameters, fitted.utr.params)[!duplicated(c(names(parameters), 
                                                            names(fitted.utr.params)), 
                                                          fromLast = TRUE)]

update_parameters = TRUE
if (update_parameters){
  updated.params <- full_soybean_ld11$parameters[names(full_soybean_ld11$parameters) %in% names(parameters)]
  full_soybean_ld11$parameters <-c(parameters, 
                                   updated.params)[!duplicated(c(names(parameters), 
                                                                 names(updated.params)), 
                                                               fromLast = TRUE)]
}else{
  full_soybean_ld11$parameters <- parameters
}

ExpBiomass <- list()
weather.afteremergence <- list()

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

for (i in 1:length(years)){
  ExpBiomass[[i]] <- loadRData(paste0('../../energy-farm-biocro/soybean_ld11_biomass_', years[i],'/soybean_ld11_biomass_', years[i], '.RData')) 
  weather <- loadRData(paste0('../../energy-farm-biocro/weather_', years[i], '/weather', years[i], '_hourly.RData'))
  weather.supplement <- loadRData(paste0('../Data/Weather_data/weather', years[i], 'supplement.RData'))
  full_soybean_ld11$parameters$Catm <- Catms[i]
  
  # update initial values
  sub_frac <- 0.1           # substrate_fraction
  str_frac <- 1 - sub_frac  # structural_fraction
  seed_mass <- ExpBiomass[[i]]$initial_seed[1]
  j <- 2
  mass_t <- sum(ExpBiomass[[i]][i, c('leaf', 'stem', 'root')])
  leaf_frac <- ExpBiomass[[i]]$leaf[j]/mass_t
  stem_frac <- ExpBiomass[[i]]$stem[j]/mass_t
  root_frac <- ExpBiomass[[i]]$root[j]/mass_t
  cf <- optim_params_short[1]
  
  initial_state <- list(
    Leaf_respiration_loss = 0.0,
    Stem_respiration_loss = 0.0,
    Root_respiration_loss = 0.0,
    Pod_respiration_loss = 0.0,
    # senescence related initial state
    Leaf_senescence_loss = 0.0,
    Stem_senescence_loss = 0.0,
    Root_senescence_loss = 0.0,
    Pod_senescence_loss = 0.0,
    # Other variables
    TTc = 0.0,
    soil_water_content = 0.32,
    cws1=                     0.32,          # dimensionless, current water status, soil layer 1
    cws2=                     0.32,          # dimensionless, current water status, soil layer 2
    Rhizome =                 0.0000001,     # Mg / ha
    RhizomeLitter =           0,               # Mg / ha
    Leaf_substrate_carbon = sub_frac * seed_mass * leaf_frac / cf,
    Leaf_structural_carbon = str_frac * seed_mass * leaf_frac / cf,
    Stem_substrate_carbon = sub_frac * seed_mass * stem_frac / cf, 
    Stem_structural_carbon = str_frac * seed_mass * stem_frac / cf,
    Root_substrate_carbon =  sub_frac * seed_mass * root_frac / cf,
    Root_structural_carbon = str_frac * seed_mass * root_frac / cf,
    Pod_substrate_carbon = 1e-4,  
    Pod_structural_carbon = 9e-4)
  
  update_init = TRUE
  if (update_init){
    updated.init <- full_soybean_ld11$initial_values[names(full_soybean_ld11$initial_values) %in% names(initial_state)]
    full_soybean_ld11$initial_values <-c(initial_state, 
                                         updated.init)[!duplicated(c(names(initial_state), 
                                                                     names(updated.init)), 
                                                                   fromLast = TRUE)]
  }else{
    full_soybean_ld11$initial_values <- initial_state
  }
  
  # Make some decisions about what to do
  MAKE_OPTIONAL_PLOTS <- TRUE
  SET_NEW_PARAMETER_VALUES <- TRUE
  VERBOSE_MODEL_VALIDATION <- FALSE
  SLA_AS_DRIVER <- TRUE
  
  first_data_time <- ExpBiomass[[i]]$time[1]
  weather.aftersowing <- weather[weather$time >= first_data_time, ]
  weather.aftersowing$DVI <- weather.supplement$DVI
  # start from emergence time
  weather.afteremergence[[i]] <- weather.aftersowing[-(1:which.min(abs(weather.aftersowing$DVI))),]
  
  if (SLA_AS_DRIVER) {
    # The experimental data indicates a non-monotonic dependence of SLA on
    # time. As a test, try including it as a driver instead of a parameter,
    # where the values are interpolated from the experimental ones.
    ExpBiomass[[i]]$SLA[1] <- ExpBiomass[[i]]$SLA[2]-
      (ExpBiomass[[i]]$time[2]-ExpBiomass[[i]]$time[1])* 
      (ExpBiomass[[i]]$SLA[3]-ExpBiomass[[i]]$SLA[2])/
      (ExpBiomass[[i]]$time[3]-ExpBiomass[[i]]$time[2])
    sla_func <- approxfun(
      ExpBiomass[[i]]$time,
      ExpBiomass[[i]]$SLA,
      rule = 2,
      na.rm = TRUE,
      method = 'linear'
    )
    weather.afteremergence[[i]]$iSp <- sla_func(weather.afteremergence[[i]]$time)
    
    full_soybean_ld11$parameters$iSp <- NULL
  }
  
  full_soybean_ld11$parameters$timestep <- 1
  full_soybean_ld11$parameters$time_zone_offset <- NULL
  full_soybean_ld11$parameters$Rd = 1.28
  
  result <- with(full_soybean_ld11, {run_biocro(
    initial_values,
    parameters,
    weather.afteremergence[[i]],
    direct_modules,
    differential_modules,
    ode_solver
  )})
  # Plot the biomass along with the measured values
  print(xyplot(
    Leaf + Stem + Root + Pod ~ time,
    data = result,
    type = 'l',
    auto = TRUE,
    grid = TRUE,
    xlab = paste0('Day of year (', years[i], ')'),
    ylab = 'Soybean biomass (Mg / ha)',
    ylim = c(-1, 6),
    panel = function(...) {
      panel.xyplot(...)
      panel.points(
        ExpBiomass[[i]]$leaf ~ ExpBiomass[[i]]$time,
        type = 'b',
        col = 'darkblue',
        pch = 16,
        lty = 2
      )
      panel.points(
        ExpBiomass[[i]]$root ~ ExpBiomass[[i]]$time,
        type = 'b',
        col = 'darkgreen',
        pch = 16,
        lty = 2
      )
      panel.points(
        ExpBiomass[[i]]$stem ~ ExpBiomass[[i]]$time,
        type = 'b',
        col = 'darkmagenta',
        pch = 16,
        lty = 2
      )
      panel.points(
        ExpBiomass[[i]]$pod ~ ExpBiomass[[i]]$time,
        type = 'b',
        col = 'darkred',
        pch = 16,
        lty = 2
      )
    }
  ))
  allocation_percentage_tall <- plot_partitioning(result, years[i])
}





