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
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')
source('../Data/Soybean-BioCro_Parameters/soybean_modules.R')
source('../ReproduceFigures/plot_partitioning.R')
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
fitted.utr.params <- optim_params_conversion(optim_params_short_SoyFACE)
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
results <- list()
figs <- list()
lai.figs <- list()
allocation.figs <- list()
pod_change <- data.frame()
pod_change_utilization_rate <- data.frame()

loadRData <- function(fileName){
  # loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
size.title <- 12
size.axislabel <-12
size.axis <- 12
size.legend <- 8
col.palette.muted <- c( "#117733", "#999933", "#332288", "#882255")

plot_biomass <- function(simulation_result, measurements, year){
  # Save plot into the list
  biocro_organ_biomass <- simulation_result[c('time', 'Leaf', 'Stem', 'Root', 'Pod')]
  biocro_organ_biomass_tall <- melt(biocro_organ_biomass, id.vars = 'time')
  names(biocro_organ_biomass_tall) <- c('time','Organ', 'biomass')
  
  field_organ_biomass <- measurements[c('time', 'leaf', 'stem', 'root', 'pod')]
  names(field_organ_biomass)[names(field_organ_biomass) %in% c('leaf', 'stem', 'root', 'pod')] <- c('Leaf', 'Stem', 'Root', 'Pod')
  field_organ_biomass_tall <- melt(field_organ_biomass, id.vars = 'time')
  names(field_organ_biomass_tall) <- c('time','Organ', 'biomass')
  
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
    scale_y_continuous(limits = c(0, 7), breaks = seq(0, 7, 2)) +
    scale_x_continuous(breaks = seq(180,280,30))+
    labs(title=element_blank(), 
         x=paste0('Day of Year (', year, ')'), 
         y='Biomass (Mg/ha)')+
    scale_color_manual(values = col.palette.muted)
}

for (i in 1:length(years)){
# for (i in 1:1){
  ExpBiomass[[i]] <- loadRData(paste0('../../energy-farm-biocro/soybean_ld11_biomass_', years[i],'/soybean_ld11_biomass_', years[i], '.RData'))
  weather <- loadRData(paste0('../../energy-farm-biocro/weather_', years[i], '/weather', years[i], '_hourly.RData'))
  full_soybean_ld11$parameters$Catm <- Catms[i]
  
  # update initial values
  sub_frac <- 0.1           # substrate_fraction
  str_frac <- 1 - sub_frac  # structural_fraction
  seed_mass <- ExpBiomass[[i]]$initial_seed[1]
  j <- 2
  mass_t <- sum(ExpBiomass[[i]][j, c('leaf', 'stem', 'root')])
  leaf_frac <- ExpBiomass[[i]]$leaf[j]/mass_t
  stem_frac <- ExpBiomass[[i]]$stem[j]/mass_t
  root_frac <- ExpBiomass[[i]]$root[j]/mass_t
  # print(paste0("fractions - leaf: ", leaf_frac, " stem: ", stem_frac, " root: ",root_frac))
  cf <- 0.3 # optim_params_short_SoyFACE[1]
  
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
  SLA_AS_DRIVER <- FALSE
  
  first_data_time <- ExpBiomass[[i]]$time[1]
  weather.aftersowing <- weather[weather$time >= first_data_time, ]
  
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
  # full_soybean_ld11$parameters$Stem_senescence_beta <- 2.1
  
  result <- with(full_soybean_ld11, {run_biocro(
    initial_values,
    parameters,
    weather.afteremergence[[i]],
    direct_modules,
    differential_modules,
    ode_solver
  )})
  results[[i]] <- result
  figs[[i]] <- plot_biomass(result, ExpBiomass[[i]], years[i])
  
  # # change the parameters by 10%
  # for (j in 1:length(arg_names)){
  #   lower_params <- full_soybean_ld11$parameters
  #   lower_params[[arg_names[j]]] <- full_soybean_ld11$parameters[[arg_names[j]]] * 0.9
  #   higher_params <- full_soybean_ld11$parameters
  #   higher_params[[arg_names[j]]] <- full_soybean_ld11$parameters[[arg_names[j]]] * 1.1
  #   
  #   result_lower_param <- with(full_soybean_ld11, {run_biocro(
  #     initial_values,
  #     lower_params,
  #     weather.afteremergence[[i]],
  #     direct_modules,
  #     differential_modules,
  #     ode_solver
  #   )})
  #   lower_param_pod_change <- (max(result_lower_param$Pod)-max(result$Pod))/max(result$Pod)
  #   
  #   result_higher_param <- with(full_soybean_ld11, {run_biocro(
  #     initial_values,
  #     higher_params,
  #     weather.afteremergence[[i]],
  #     direct_modules,
  #     differential_modules,
  #     ode_solver
  #   )})
  #   higher_param_pod_change <- (max(result_higher_param$Pod)-max(result$Pod))/max(result$Pod)
  #   pod_change[arg_names[j],paste0(years[i],' -10%')] <- paste0(round(lower_param_pod_change*100, 2), "%")
  #   pod_change[arg_names[j],paste0(years[i],' +10%')] <- paste0(round(higher_param_pod_change*100, 2), "%")
  # }
  # print(pod_change)
  # write.csv(pod_change, "LD11_pod_sensitivity_10percent.csv")
  # 
  for (j in c(1:4,12:14)){
    lower_params <- full_soybean_ld11$parameters
    lower_params[[arg_names[j]]] <- full_soybean_ld11$parameters[[arg_names[j]]] * 0.5
    higher_params <- full_soybean_ld11$parameters
    higher_params[[arg_names[j]]] <- full_soybean_ld11$parameters[[arg_names[j]]] * 2

    result_lower_param <- with(full_soybean_ld11, {run_biocro(
      initial_values,
      lower_params,
      weather.afteremergence[[i]],
      direct_modules,
      differential_modules,
      ode_solver
    )})
    lower_param_pod_change <- (max(result_lower_param$Pod)-max(result$Pod))/max(result$Pod)

    result_higher_param <- with(full_soybean_ld11, {run_biocro(
      initial_values,
      higher_params,
      weather.afteremergence[[i]],
      direct_modules,
      differential_modules,
      ode_solver
    )})
    higher_param_pod_change <- (max(result_higher_param$Pod)-max(result$Pod))/max(result$Pod)
    pod_change_utilization_rate[arg_names[j],paste0(years[i],' x0.01')] <- paste0(round(lower_param_pod_change*100, 2), "%")
    pod_change_utilization_rate[arg_names[j],paste0(years[i],' x100')] <- paste0(round(higher_param_pod_change*100, 2), "%")
  }
  print(pod_change_utilization_rate)
  write.csv(pod_change_utilization_rate, "LD11_pod_sensitivity_10000percent.csv")
}