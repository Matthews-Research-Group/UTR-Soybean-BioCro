yr <- '2003'

weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))

co2_opt <- '_ambient_'
parameters$Catm <- 372
ExpBiomass <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt , 'biomass.csv'))
colnames(ExpBiomass)<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

ExpBiomass.std <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
colnames(ExpBiomass.std)<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

sd.idx <- which(weather$doy == 147)[12] # Morgan et al. 2005
hd.ind <- which(weather$doy == max(ExpBiomass$DOY))[24]
defoliation.ind <- which(weather$doy == 198)[14]

weather.growingseason <- weather[sd.idx: hd.ind,]
weather.growingseason_1 <- weather[sd.idx:defoliation.ind,]
weather.growingseason_2 <- weather[defoliation.ind:hd.ind,]

RootVals <- data.frame("DOY"=ExpBiomass$DOY[3], "Root"=0.17*sum(ExpBiomass[5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130

# define a function to update differential values
update_differential_quantities <- function(r, updated_values, UPDATE_PARAMETERS, model){
  ## Define remaining percentages
  leaf_remaining_percentage <- 0.4
  stem_remaining_percentage <- 0.5
  
  last_row <- nrow(r)
  differential_quantities_just_before_defoliation <-
    as.list(r[last_row, updated_values])
  
  differential_quantities_just_after_defoliation <-
    differential_quantities_just_before_defoliation
  
  if (model == 'utr'){
    differential_quantities_just_after_defoliation$Leaf_substrate_carbon <-
      differential_quantities_just_before_defoliation$Leaf_substrate_carbon * leaf_remaining_percentage
    differential_quantities_just_after_defoliation$Leaf_structural_carbon <-
      differential_quantities_just_before_defoliation$Leaf_structural_carbon * leaf_remaining_percentage # Could be changed to a different percentage
    
    differential_quantities_just_after_defoliation$Stem_substrate_carbon <-
      differential_quantities_just_before_defoliation$Stem_substrate_carbon * stem_remaining_percentage
    differential_quantities_just_after_defoliation$Stem_structural_carbon <-
      differential_quantities_just_before_defoliation$Stem_structural_carbon * stem_remaining_percentage
  } 
  else if (model == 'partitioning'){
    differential_quantities_just_after_defoliation$Leaf <-
      differential_quantities_just_before_defoliation$Leaf * leaf_remaining_percentage
    
    differential_quantities_just_after_defoliation$Stem <-
      differential_quantities_just_before_defoliation$Stem * stem_remaining_percentage
  }
  
  if(UPDATE_PARAMETERS){
    differential_quantities_just_after_defoliation$DVI <-
      differential_quantities_just_before_defoliation$DVI - 0.2
  }
  
  return(differential_quantities_just_after_defoliation)
  
}
soybean$parameters$time_zone_offset <- -6
# ##### No hail scenario #####
# soybean_utr_optsolver_no_hail <- partial_run_biocro(initial_values,
#                                                     parameters,
#                                                     weather.growingseason,
#                                                     direct_modules,
#                                                     differential_modules,
#                                                     solver,
#                                                     arg_names,
#                                                     verbose = FALSE)
# 
# result_utr_no_hail <- soybean_utr_optsolver_no_hail(optim_params_conversion(optim_params_short_SoyFACE))
# result_partitioning_no_hail <- with(soybean, run_biocro(initial_values,
#                                                         parameters,
#                                                         weather.growingseason,
#                                                         direct_modules,
#                                                         differential_modules,
#                                                         solver,
#                                                         verbose = FALSE))
# 
# # plot simulation vs observation biomass and calculate mse
# fig_2003_utr_no_hail <- plot_SoyFACE_biomass(result_utr_no_hail, ExpBiomass, ExpBiomass.std, co2_opt, yr)
# result_partitioning_no_hail$Pod <- result_partitioning_no_hail$Grain + result_partitioning_no_hail$Shell
# fig_2003_partitioning_no_hail <- plot_SoyFACE_biomass(result_partitioning_no_hail, ExpBiomass, ExpBiomass.std, co2_opt, yr)
# combined_no_hail <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
#                                  arrangeGrob(fig_2003_partitioning_no_hail + theme(legend.position="none"), top = 'Partitioning Model'),
#                                  arrangeGrob(fig_2003_utr_no_hail + theme(legend.position="none"), top = 'UTR model'),
#                                  common_legend,
#                                  ncol = 4, widths = c(0.1, 5, 5, 1))

#### Before Hail #### 
soybean_utr_optsolver_hail1 <- partial_run_biocro(initial_values,
                                                  parameters,
                                                  weather.growingseason_1,
                                                  direct_modules,
                                                  differential_modules,
                                                  solver,
                                                  arg_names,
                                                  verbose = FALSE)
result_utr_hail1 <- soybean_utr_optsolver_hail1(optim_params_conversion(optim_params_short_SoyFACE))
result_partitioning_hail1 <- with(soybean, run_biocro(initial_values,
                                                       parameters,
                                                       weather.growingseason_1,
                                                       direct_modules,
                                                       differential_modules,
                                                       solver,
                                                       verbose = FALSE))
#### After Hail #### 
### Without Parameter Update ####
## UTR ##
UPDATE_PARAMETERS <- FALSE
differential_quantities_just_after_defoliation_utr <- update_differential_quantities(result_utr_hail1, 
                                                                                 names(initial_values), 
                                                                                 UPDATE_PARAMETERS, 'utr')

soybean_utr_optsolver_hail2 <- partial_run_biocro(differential_quantities_just_after_defoliation_utr,
                                                  parameters,
                                                  weather.growingseason_2,
                                                  direct_modules,
                                                  differential_modules,
                                                  solver,
                                                  arg_names,
                                                  verbose = FALSE)

orignal_utr_params <- data.frame(optim_params_conversion(optim_params_short_SoyFACE))
rownames(orignal_utr_params) <- arg_names
colnames(orignal_utr_params) <- "Value"
parameters_after_hail <- orignal_utr_params

result_utr_hail2 <- soybean_utr_optsolver_hail2(parameters_after_hail$Value)

## Partitioning model##
differential_quantities_just_after_defoliation_partitioning <- update_differential_quantities(result_partitioning_hail1, 
                                                                                 names(soybean$initial_values),
                                                                                 UPDATE_PARAMETERS, 'partitioning')

result_partitioning_hail2 <- run_biocro(differential_quantities_just_after_defoliation_partitioning,
                                         soybean$parameters,
                                         weather.growingseason_2,
                                         soybean$direct_modules,
                                         soybean$differential_modules,
                                         soybean$ode_solver,
                                         verbose = FALSE)

result_utr_hail <- rbind(result_utr_hail1[seq_len(nrow(result_utr_hail1) - 1), ], result_utr_hail2)
result_partitioning_hail <- rbind(result_partitioning_hail1[seq_len(nrow(result_partitioning_hail1) - 1), ], result_partitioning_hail2)
result_partitioning_hail$Pod <- result_partitioning_hail$Grain + result_partitioning_hail$Shell
fig_2003_utr_with_hail <- plot_SoyFACE_biomass(result_utr_hail, ExpBiomass, ExpBiomass.std, co2_opt, yr)
fig_2003_partitioning_with_hail <- plot_SoyFACE_biomass(result_partitioning_hail, ExpBiomass, ExpBiomass.std, co2_opt, yr)

### With Updated Parameters ###
## UTR ##
UPDATE_PARAMETERS <- TRUE
differential_quantities_just_after_defoliation_utr_updatedParams <- update_differential_quantities(result_utr_hail1, 
                                                                                     names(initial_values), 
                                                                                     UPDATE_PARAMETERS, 'utr')
updated_parameters <- parameters
updated_parameters$Leaf_respiration_factor <- 0.2
updated_parameters$Stem_respiration_factor <- updated_parameters$Stem_respiration_factor + 0.2

soybean_utr_optsolver_hail2_updatedParams <- partial_run_biocro(differential_quantities_just_after_defoliation_utr_updatedParams,
                                                  updated_parameters,
                                                  weather.growingseason_2,
                                                  direct_modules,
                                                  differential_modules,
                                                  solver,
                                                  arg_names,
                                                  verbose = FALSE)

result_utr_hail2_updatedParams <- soybean_utr_optsolver_hail2_updatedParams(parameters_after_hail$Value)

## Partitioning model ##
differential_quantities_just_after_defoliation_partitioning_updatedParams <- update_differential_quantities(result_partitioning_hail1, 
                                                                                              names(soybean$initial_values),
                                                                                              UPDATE_PARAMETERS, 'partitioning')

updated_soybean_parameters <- soybean$parameters
updated_soybean_parameters$mrc_stem <- updated_soybean_parameters$mrc_stem *2
result_partitioning_hail2_updatedParams <- run_biocro(differential_quantities_just_after_defoliation_partitioning_updatedParams,
                                         updated_soybean_parameters,
                                         weather.growingseason_2,
                                         soybean$direct_modules,
                                         soybean$differential_modules,
                                         soybean$ode_solver,
                                         verbose = FALSE)

result_utr_hail_updatedParams <- rbind(result_utr_hail1[seq_len(nrow(result_utr_hail1) - 1), ], 
                                       result_utr_hail2_updatedParams)
result_partitioning_hail_updatedParams <- rbind(result_partitioning_hail1[seq_len(nrow(result_partitioning_hail1) - 1), ], 
                                                result_partitioning_hail2_updatedParams)
result_partitioning_hail_updatedParams$Pod <- result_partitioning_hail_updatedParams$Grain + 
                                              result_partitioning_hail_updatedParams$Shell
fig_2003_utr_with_hail_updatedParams <- plot_SoyFACE_biomass(result_utr_hail_updatedParams, ExpBiomass, ExpBiomass.std, co2_opt, yr)
fig_2003_partitioning_with_hail_updatedParams <- plot_SoyFACE_biomass(result_partitioning_hail_updatedParams, ExpBiomass, ExpBiomass.std, co2_opt, yr)

########## Elevated CO2##########
co2_opt <- '_CO2_'
parameters$Catm <- 550
soybean$parameters$Catm <- 550

ExpBiomass.elevCO2 <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt , 'biomass.csv'))
colnames(ExpBiomass.elevCO2)<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

ExpBiomass.elevCO2.std <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
colnames(ExpBiomass.elevCO2.std)<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

RootVals<- data.frame("DOY"=ExpBiomass.elevCO2$DOY[3], "Root"=0.17*sum(ExpBiomass.elevCO2[5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130


# ## No hail
# soybean_utr_optsolver_no_hail_eCO2 <- partial_run_biocro(initial_values,
#                                                          parameters,
#                                                          weather.growingseason,
#                                                          direct_modules,
#                                                          differential_modules,
#                                                          solver,
#                                                          arg_names,
#                                                          verbose = FALSE)
# 
# result_utr_no_hail_eCO2 <- soybean_utr_optsolver_no_hail_eCO2(optim_params_conversion(optim_params_short_SoyFACE))
# 
# result_partitioning_no_hail_eCO2 <- with(soybean, run_biocro(initial_values,
#                                                              parameters,
#                                                              weather.growingseason,
#                                                              direct_modules,
#                                                              differential_modules,
#                                                              solver,
#                                                              verbose = FALSE))
# 
# # plot simulation vs observation biomass and calculate mse
# fig_2003_utr_no_hail_eCO2 <- plot_SoyFACE_biomass(result_utr_no_hail_eCO2, ExpBiomass.elevCO2, ExpBiomass.elevCO2.std, co2_opt, yr)
# result_partitioning_no_hail_eCO2$Pod <- result_partitioning_no_hail_eCO2$Grain + result_partitioning_no_hail_eCO2$Shell
# result_partitioning_no_hail_eCO2 <- plot_SoyFACE_biomass(result_partitioning_no_hail_eCO2, ExpBiomass.elevCO2, ExpBiomass.elevCO2.std, co2_opt, yr)
# combined_no_hail_eCO2 <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
#                                       arrangeGrob(result_partitioning_no_hail_eCO2 + theme(legend.position="none"), top = 'Partitioning Model'),
#                                       arrangeGrob(fig_2003_utr_no_hail_eCO2 + theme(legend.position="none"), top = 'UTR model'),
#                                       common_legend,
#                                       ncol = 4, widths = c(0.1, 5, 5, 1))

#### With Hail ####
## Before Hail ##
soybean_utr_optsolver_eCO2_hail1 <- partial_run_biocro(initial_values,
                                                       parameters,
                                                       weather.growingseason_1,
                                                       direct_modules,
                                                       differential_modules,
                                                       solver,
                                                       arg_names,
                                                       verbose = FALSE)
result_utr_eCO2_hail1 <- soybean_utr_optsolver_eCO2_hail1(optim_params_conversion(optim_params_short_SoyFACE))
result_partitioning_eCO2_hail1 <- with(soybean, run_biocro(initial_values,
                                                           parameters,
                                                           weather.growingseason_1,
                                                           direct_modules,
                                                           differential_modules,
                                                           solver,
                                                           verbose = FALSE))
## Without Parameter Update ##
# UTR
UPDATE_PARAMETERS <- FALSE
differential_quantities_just_after_defoliation_utr_eCO2 <- update_differential_quantities(result_utr_eCO2_hail1, 
                                                                                 names(initial_values), 
                                                                                 UPDATE_PARAMETERS, 'utr')

soybean_utr_optsolver_eCO2_hail2 <- partial_run_biocro(differential_quantities_just_after_defoliation_utr_eCO2,
                                                       parameters,
                                                       weather.growingseason_2,
                                                       direct_modules,
                                                       differential_modules,
                                                       solver,
                                                       arg_names,
                                                       verbose = FALSE)
result_utr_eCO2_hail2 <- soybean_utr_optsolver_eCO2_hail2(parameters_after_hail$Value)

differential_quantities_just_before_defoliation_partitioning_eCO2 <- update_differential_quantities(result_partitioning_eCO2_hail1, 
                                                                                  names(soybean$initial_values), 
                                                                                  UPDATE_PARAMETERS, 'partitioning')
# Partitioning
result_partitioning_eCO2_hail2 <- run_biocro(differential_quantities_just_before_defoliation_partitioning_eCO2,
                                             soybean$parameters,
                                             weather.growingseason_2,
                                             soybean$direct_modules,
                                             soybean$differential_modules,
                                             soybean$ode_solver,
                                             verbose = FALSE)

result_utr_eCO2_hail <- rbind(result_utr_eCO2_hail1[seq_len(nrow(result_utr_eCO2_hail1) - 1), ], result_utr_eCO2_hail2)
result_partitioning_eCO2_hail <- rbind(result_partitioning_eCO2_hail1[seq_len(nrow(result_partitioning_eCO2_hail1) - 1), ], result_partitioning_eCO2_hail2)
result_partitioning_eCO2_hail$Pod <- result_partitioning_eCO2_hail$Grain + result_partitioning_eCO2_hail$Shell
fig_2003_utr_eCO2_with_hail <- plot_SoyFACE_biomass(result_utr_eCO2_hail, ExpBiomass.elevCO2, ExpBiomass.elevCO2.std, co2_opt, yr)
fig_2003_partitioning_eCO2_with_hail <- plot_SoyFACE_biomass(result_partitioning_eCO2_hail, ExpBiomass.elevCO2, ExpBiomass.elevCO2.std, co2_opt, yr)

### With Updated Parameters ###
## UTR ##
UPDATE_PARAMETERS <- TRUE
differential_quantities_just_after_defoliation_utr_eCO2_updatedParams <- update_differential_quantities(result_utr_eCO2_hail1, 
                                                                                                        names(initial_values), 
                                                                                                        UPDATE_PARAMETERS, 'utr')
updated_parameters <- parameters
updated_parameters$Leaf_respiration_factor <- 0.2
updated_parameters$Stem_respiration_factor <- updated_parameters$Stem_respiration_factor + 0.2

soybean_utr_eCO2_optsolver_hail2_updatedParams <- partial_run_biocro(differential_quantities_just_after_defoliation_utr_eCO2_updatedParams,
                                                                     updated_parameters,
                                                                     weather.growingseason_2,
                                                                     direct_modules,
                                                                     differential_modules,
                                                                     solver,
                                                                     arg_names,
                                                                     verbose = FALSE)

result_utr_eCO2_hail2_updatedParams <- soybean_utr_eCO2_optsolver_hail2_updatedParams(parameters_after_hail$Value)

## Partitioning model ##
differential_quantities_just_after_defoliation_partitioning_eCO2_updatedParams <- update_differential_quantities(result_partitioning_eCO2_hail1, 
                                                                                                                 names(soybean$initial_values),
                                                                                                                 UPDATE_PARAMETERS, 'partitioning')

updated_soybean_parameters <- soybean$parameters
updated_soybean_parameters$mrc_stem <- updated_soybean_parameters$mrc_stem * 2
result_partitioning_eCO2_hail2_updatedParams <- run_biocro(differential_quantities_just_after_defoliation_partitioning_eCO2_updatedParams,
                                                           updated_soybean_parameters,
                                                           weather.growingseason_2,
                                                           soybean$direct_modules,
                                                           soybean$differential_modules,
                                                           soybean$ode_solver,
                                                           verbose = FALSE)

result_utr_eCO2_hail_updatedParams <- rbind(result_utr_eCO2_hail1[seq_len(nrow(result_utr_eCO2_hail1) - 1), ], 
                                            result_utr_eCO2_hail2_updatedParams)
result_partitioning_eCO2_hail_updatedParams <- rbind(result_partitioning_eCO2_hail1[seq_len(nrow(result_partitioning_eCO2_hail1) - 1), ], 
                                                     result_partitioning_eCO2_hail2_updatedParams)
result_partitioning_eCO2_hail_updatedParams$Pod <- result_partitioning_eCO2_hail_updatedParams$Grain + 
  result_partitioning_eCO2_hail_updatedParams$Shell
fig_2003_utr_eCO2_with_hail_updatedParams <- plot_SoyFACE_biomass(result_utr_eCO2_hail_updatedParams, ExpBiomass, ExpBiomass.std, co2_opt, yr)
fig_2003_partitioning_eCO2_with_hail_updatedParams <- plot_SoyFACE_biomass(result_partitioning_eCO2_hail_updatedParams, ExpBiomass, ExpBiomass.std, co2_opt, yr)

plot_partitioning(result_utr_hail, '2003') 
combined_with_hail <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
                                   arrangeGrob(
                                     arrangeGrob(
                                       arrangeGrob(fig_2003_partitioning_with_hail +
                                                     theme(axis.title.x = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           legend.position="none"),
                                                   top = 'Partitioning Model'),
                                       arrangeGrob(fig_2003_utr_with_hail +
                                                     theme(axis.title.x = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           legend.position="none"),
                                                   top = 'UTR model'),
                                       ncol = 2, top = '(A) Ambient CO2'),
                                     arrangeGrob(
                                       arrangeGrob(fig_2003_partitioning_eCO2_with_hail  +
                                                     theme(axis.title.x = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           legend.position="none")),
                                       arrangeGrob(fig_2003_utr_eCO2_with_hail  +
                                                     theme(axis.title.x = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           legend.position="none")),
                                       ncol = 2, top = '(B) Elevated CO2'),
                                     arrangeGrob(textGrob('Day of Year (2003)')),
                                     nrow = 3, heights = c(4, 4, 0.3)),
                                   common_legend,
                                   ncol = 3, widths = c(0.3, 5, 1))
# 
# combined_with_hail_v2 <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
#                                    arrangeGrob(
#                                      arrangeGrob(
#                                        arrangeGrob(fig_2003_utr_with_hail +
#                                                      theme(axis.title.x = element_blank(),
#                                                            axis.title.y = element_blank(),
#                                                            legend.position="none"),
#                                                    top = 'Ambient CO2'),
#                                        arrangeGrob(fig_2003_utr_eCO2_with_hail +
#                                                      theme(axis.title.x = element_blank(),
#                                                            axis.title.y = element_blank(),
#                                                            legend.position="none"),
#                                                    top = 'Elevated CO2'),
#                                        ncol= 2, top = 'Without Parameter Adjustment'),
#                                      arrangeGrob(
#                                        arrangeGrob(fig_2003_utr_with_hail_updatedParams  +
#                                                      theme(axis.title.x = element_blank(),
#                                                            axis.title.y = element_blank(),
#                                                            legend.position="none")),
#                                        arrangeGrob(fig_2003_utr_eCO2_with_hail_updatedParams  +
#                                                      theme(axis.title.x = element_blank(),
#                                                            axis.title.y = element_blank(),
#                                                            legend.position="none")),
#                                        ncol = 2, top = 'With Parameter Adjustment'),
#                                      arrangeGrob(textGrob('Day of Year (2003)')),
#                                      nrow = 3, heights = c(4, 4, 0.3)),
#                                    common_legend,
#                                    ncol = 3, widths = c(0.3, 5, 1))

combined_with_hail_utr <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
                                   arrangeGrob(
                                     arrangeGrob(
                                       arrangeGrob(fig_2003_utr_with_hail +
                                                     theme(axis.title.x = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           legend.position="none"),
                                                   top = '(A) Ambient CO2'),
                                       arrangeGrob(fig_2003_utr_eCO2_with_hail +
                                                     theme(axis.title.x = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           legend.position="none"),
                                                   top = '(B) Elevated CO2'),
                                       ncol= 2),
                                     arrangeGrob(textGrob('Day of Year (2003)')),
                                     nrow = 2, heights = c(4, 0.3)),
                                   common_legend,
                                   ncol = 3, widths = c(0.3, 5, 1))

combined_with_hail_partitioning <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
                                      arrangeGrob(
                                        arrangeGrob(
                                          arrangeGrob(fig_2003_partitioning_with_hail +
                                                        theme(axis.title.x = element_blank(),
                                                              axis.title.y = element_blank(),
                                                              legend.position="none"),
                                                      top = '(A) Ambient CO2'),
                                          arrangeGrob(fig_2003_partitioning_eCO2_with_hail +
                                                        theme(axis.title.x = element_blank(),
                                                              axis.title.y = element_blank(),
                                                              legend.position="none"),
                                                      top = '(B) Elevated CO2'),
                                          ncol= 2),
                                        arrangeGrob(textGrob('Day of Year (2003)')),
                                        nrow = 2, heights = c(4, 0.3)),
                                      common_legend,
                                      ncol = 3, widths = c(0.3, 5, 1))

ggsave('Fig7-hail.png',
       plot = combined_with_hail_utr,
       width = 5,
       height = 2.2,
       units = "in",
       dpi = 600
)

ggsave('FigS-hail-partitioning.png',
       plot = combined_with_hail_partitioning,
       width = 5,
       height = 2.2,
       units = "in",
       dpi = 600
)


df1 <- result_utr_hail %>%
  group_by(doy) %>%
  summarise(
    Leaf = sum(`Leaf_utilization_rate`, na.rm = TRUE),
    Stem = sum(`Stem_utilization_rate`, na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(cols = c(Leaf, Stem), names_to = "type", values_to = "value") %>%
  mutate(scenario = "Ambient")

df2 <- result_utr_eCO2_hail %>%
  group_by(doy) %>%
  summarise(
    Leaf = sum(`Leaf_utilization_rate`, na.rm = TRUE),
    Stem = sum(`Stem_utilization_rate`, na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(cols = c(Leaf, Stem), names_to = "type", values_to = "value") %>%
  mutate(scenario = "Elevated")

hail_utr_utilization_plot <-
  bind_rows(df1, df2) %>%
  ggplot(aes(x = doy, y = value, color = type, linetype = scenario)) +
  annotate("rect", xmin = 197.5, xmax = 198.5, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.4) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("Leaf" = "#117733", "Stem" = "#999933"),
    name = "Utilization Rate"
  ) +
  scale_linetype_manual(
    values = c("Ambient" = "solid", "Elevated" = "dotted"),
    name = "CO2 Level"
  ) +
  labs(
    title = "Leaf and Stem Utilization Rate by Day of Year",
    x = "2003 Day of Year (DOY)",
    y = expression(paste("Utilization Rate", "(mol m"^{-2}, "day"^{-1},")"))
  ) +
  theme_minimal()
ggsave('FigS-hail-utilization.png',
       plot = hail_utr_utilization_plot,
       width = 6,
       height = 2.5,
       units = "in",
       dpi = 600
)

