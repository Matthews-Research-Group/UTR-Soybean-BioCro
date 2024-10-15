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
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')

# initialize lists for figures
co2_opt = '_co2_'

# names of fitted parameters
setwd('../Data/Soybean-BioCro_Parameters')
# load parameter files
param_files <- list.files(pattern = "[.]R$", recursive = TRUE)

param_files <- param_files[-c(which(param_files == "soybean_initial_values.R"))]
lapply(param_files, source)
# load parameter files
source('soybean_initial_values.R')

years <- c('2002', '2004', '2005', '2006')
for (yr in years){
  weather <- read.csv(file = paste0('../Weather_data/', yr,'_Bondville_IL_daylength_wDVI.csv'))
  emergence.ind <- which(weather$DVI>0)[1]
  hd.ind <- which(weather$doy == 289)[24]
  defoliation.ind <- which.min(abs((weather$DVI-1.5)))
  initial_state$DVI <- NULL
  
  weather.growingseason <- weather[emergence.ind:hd.ind,]
  weather.growingseason_1 <- weather[emergence.ind:defoliation.ind,]
  weather.growingseason_2 <- weather[defoliation.ind:hd.ind,]
  
  ExpBiomass <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass.csv'))
  colnames(ExpBiomass)<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  ExpBiomass.std <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
  colnames(ExpBiomass.std)<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  RootVals <- data.frame("DOY"=ExpBiomass$DOY[3], "Root"=0.17*sum(ExpBiomass[5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  numrows <- nrow(weather.growingseason)
  invwts <- ExpBiomass.std
  
  soybean_optsolver_no_dt <- partial_run_biocro(initial_state, # dt stands for defoliation and thinning
                                                parameters,
                                                weather.growingseason,
                                                steady_state_module_names,
                                                derivative_module_names,
                                                solver,
                                                arg_names,
                                                verbose = FALSE)
  updated_params <- optim_params_conversion(optim_params_short_SoyFACE)
  names(updated_params) <- arg_names
  result_no_dt <- soybean_optsolver_no_dt(updated_params)
  yield_no_dt <- max(result_no_dt['Pod'])
  
  # xyplot(data=result_no_dt, Leaf+Stem+Root+Pod~time, auto=TRUE)
  
  
  soybean_optsolver_1 <- partial_run_biocro(initial_state,
                                            parameters,
                                            weather.growingseason_1,
                                            steady_state_module_names,
                                            derivative_module_names,
                                            solver,
                                            arg_names,
                                            verbose = FALSE)
  # updated_params[28] <- 1.15
  result_1 <- soybean_optsolver_1(optim_params_conversion(optim_params_short_SoyFACE))
  # xyplot(data=result_1, Leaf+Stem+Root+Pod~time, auto=TRUE)
  # Get the final values of the differential quantities; these will be the
  # values just before defoliation
  differential_quantities_just_before_defoliation <-
    as.list(result_1[nrow(result_1), names(initial_state)])
  
  differential_quantities_just_after_defoliation <-
    differential_quantities_just_before_defoliation
  # # Now reduce the leaf mass
  remaining_leaf_percent <- 0.67
  stem_new_percentage <- 1.0
  root_new_percentage <- 1.0
  pod_new_percentage <- 1.0
  
  differential_quantities_just_after_defoliation$Leaf_substrate_carbon <-
    differential_quantities_just_before_defoliation$Leaf_substrate_carbon * remaining_leaf_percent
  differential_quantities_just_after_defoliation$Leaf_structural_carbon <-
    differential_quantities_just_before_defoliation$Leaf_structural_carbon * remaining_leaf_percent # Could be changed to a different percentage
  
  leaf_C_before_defoliation <- differential_quantities_just_before_defoliation$Leaf_substrate_carbon+
    differential_quantities_just_before_defoliation$Leaf_structural_carbon
  stem_C_before_defoliation <- differential_quantities_just_before_defoliation$Stem_substrate_carbon+
    differential_quantities_just_before_defoliation$Stem_structural_carbon
  abg_C_before_defoliation <- leaf_C_before_defoliation + stem_C_before_defoliation
  
  leaf_C_after_defoliation <- differential_quantities_just_after_defoliation$Leaf_substrate_carbon+
    differential_quantities_just_after_defoliation$Leaf_structural_carbon
  
  # If there is stem loss, then also reduce stem biomass.
  if (stem_new_percentage < 1){
    differential_quantities_just_after_defoliation$Stem_substrate_carbon <-
      differential_quantities_just_before_defoliation$Stem_substrate_carbon * stem_new_percentage
    differential_quantities_just_after_defoliation$Stem_structural_carbon <-
      differential_quantities_just_before_defoliation$Stem_structural_carbon * stem_new_percentage
  }
  
  if (root_new_percentage < 1){
    differential_quantities_just_after_defoliation$Root_substrate_carbon <-
      differential_quantities_just_before_defoliation$Root_substrate_carbon * root_new_percentage
    differential_quantities_just_after_defoliation$Root_structural_carbon <-
      differential_quantities_just_before_defoliation$Root_structural_carbon * root_new_percentage
  }
  
  if (pod_new_percentage < 1){
    differential_quantities_just_after_defoliation$Pod_substrate_carbon <-
      differential_quantities_just_before_defoliation$Pod_substrate_carbon * pod_new_percentage
    differential_quantities_just_after_defoliation$Pod_structural_carbon <-
      differential_quantities_just_before_defoliation$Pod_structural_carbon * pod_new_percentage
  }
  
  orignal_utr_params <- data.frame(optim_params_conversion(optim_params_short_SoyFACE))
  rownames(orignal_utr_params) <- arg_names
  colnames(orignal_utr_params) <- "Value"
  parameters_after_hail <- orignal_utr_params
  # parameters_after_hail['Stem_respiration_factor', 'Value'] <- orignal_utr_params['Stem_respiration_factor', 'Value'] * 2
  # parameters_after_hail['Pod_start_dvi', 'Value'] <- orignal_utr_params['Pod_start_dvi', 'Value'] + 0.1
  
  soybean_optsolver_2 <- partial_run_biocro(differential_quantities_just_after_defoliation,
                                            parameters,
                                            weather.growingseason_2,
                                            steady_state_module_names,
                                            derivative_module_names,
                                            solver,
                                            arg_names,
                                            verbose = FALSE)
  
  
  result_2 <- soybean_optsolver_2(parameters_after_hail$Value)
  result <- rbind(result_1[seq_len(nrow(result_1) - 1), ], result_2)
  # xyplot(data=result, Leaf+Stem+Root+Pod~time, auto=TRUE)
  yield_with_dt <- max(result$Pod)
  yield_percent_change <- 100 * (yield_with_dt - yield_no_dt)/yield_no_dt
  print(paste0(yr, " change: ", yield_percent_change, '%'))
}