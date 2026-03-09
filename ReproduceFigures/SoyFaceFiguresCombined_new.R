library(BioCro)
library(UTRSoybeanBML)
library(BioCroWater)
library(ggplot2)
library(grid)
library(gridExtra)
library(lattice)
library(dplyr)
library(grid)

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
source('plot_lai_comparison.R')
source('plot_biomass.R')

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
initial_values$DVI <- -1

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
ExpBiomass.elevCO2.std <- list()
RootVals <- list()
numrows <- vector()
figs <- list()
figs.elevCO2 <- list()
lai.figs <- list()
pod_sensitivity <- data.frame()
allocation.figs <- list()
allocation.elevCO2.figs <- list()

for (i in 1:length(years)) {   
  yr <- years[i]
  weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))
  sd.idx <- which(weather$doy == sow.date[i])[12]
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather.growingseason[[i]] <- weather[sd.idx: hd.ind,]
  
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
  
  # plot simulation vs observation biomass and calculate mse
  figs[[i]] <- plot_SoyFACE_biomass(result, ExpBiomass[[i]], ExpBiomass.std[[i]], co2_opt, years[i])
  
  # plot partitioning
  allocation.figs[[i]] <- plot_partitioning(result, years[i])
}


# initialize lists for figures
co2_opt = '_co2_'
parameters$Catm      <- 550

# Initialize lists
soybean_optsolver <- list()
RootVals <- list()
numrows <- vector()

for (i in 1:length(years)) {  
  yr <- years[i]
  
  ExpBiomass.elevCO2[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass.csv'))
  colnames(ExpBiomass.elevCO2[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  ExpBiomass.elevCO2.std[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
  colnames(ExpBiomass.elevCO2.std[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
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
  
  # plot simulation vs observation biomass and calculate mse
  figs.elevCO2[[i]] <- plot_SoyFACE_biomass(result, ExpBiomass.elevCO2[[i]], ExpBiomass.elevCO2.std[[i]], co2_opt, years[i])
  
  # plot partitioning
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
combined_graph_SoyFACE <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
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

# Plot lai
for (i in 2:4){
  # results[[i]]$lai <- results[[i]]$lai * 1.12
  # results.elevCO2[[i]]$lai <- results.elevCO2[[i]]$lai * 1.12
  lai.figs[[i-1]] <- plot_amb_elev_lai(results[[i]], results.elevCO2[[i]], years[[i]], LAI[[i]], LAI.elevCO2[[i]])
  # print(lai.figs[[i-1]])
}

combined_graph.lai <- grid.arrange(arrangeGrob(textGrob(bquote("LAI"~(m^2~"/"~m^2)), rot = 90, gp=gpar(fontsize=12))),
                                   arrangeGrob(arrangeGrob(lai.figs[[1]] + theme(legend.position="none"),
                                                           lai.figs[[2]] + theme(legend.position="none"),
                                                           lai.figs[[3]] + theme(legend.position="none"),
                                                           ncol = 3),
                                               ncol = 1),
                                   ncol=2, widths=c(0.3, 5))

# for (i in 1:4){
#   shoot_root_ratio = data.frame(
#     fractional_doy = results[[i]]$fractional_doy,
#     ambient_ratio = with(results[[i]], (Leaf+Stem+Pod)/Root),
#     elevated_ratio = with(results.elevCO2[[i]], (Leaf+Stem+Pod)/Root))
#   print(years[i])
#   print(xyplot(data=shoot_root_ratio, ambient_ratio + elevated_ratio ~ fractional_doy, auto.key = TRUE))
# }



#################################### Plot hail event ####################################
yr <- '2003'
co2_opt <- '_ambient_'
parameters$Catm <- 372
weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))

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

RootVals<- data.frame("DOY"=ExpBiomass$DOY[3], "Root"=0.17*sum(ExpBiomass[5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130

numrows<- nrow(weather.growingseason)
invwts <- ExpBiomass.std

soybean_utr_optsolver_no_hail <- partial_run_biocro(initial_values,
                                        parameters,
                                        weather.growingseason,
                                        direct_modules,
                                        differential_modules,
                                        solver,
                                        arg_names,
                                        verbose = FALSE)

result_utr_no_hail <- soybean_utr_optsolver_no_hail(optim_params_conversion(optim_params_short_SoyFACE)) 
soybean$parameters$time_zone_offset <- -6
result_partitioning_no_hail <- with(soybean, run_biocro(initial_values,
                                                        parameters,
                                                        weather.growingseason,
                                                        direct_modules,
                                                        differential_modules,
                                                        solver,
                                                        verbose = FALSE))

# plot simulation vs observation biomass and calculate mse
fig_2003_utr_no_hail <- plot_SoyFACE_biomass(result_utr_no_hail, ExpBiomass, ExpBiomass.std, co2_opt, yr)
result_partitioning_no_hail$Pod <- result_partitioning_no_hail$Grain + result_partitioning_no_hail$Shell
fig_2003_partitioning_no_hail <- plot_SoyFACE_biomass(result_partitioning_no_hail, ExpBiomass, ExpBiomass.std, co2_opt, yr)
combined_no_hail <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
                                 arrangeGrob(fig_2003_partitioning_no_hail + theme(legend.position="none"), top = 'Partitioning Model'),
                                 arrangeGrob(fig_2003_utr_no_hail + theme(legend.position="none"), top = 'UTR model'),
                                 common_legend,
                                 ncol = 4, widths = c(0.1, 5, 5, 1))

soybean_utr_optsolver_hail1 <- partial_run_biocro(initial_values,
                                          parameters,
                                          weather.growingseason_1,
                                          direct_modules,
                                          differential_modules,
                                          solver,
                                          arg_names,
                                          verbose = FALSE)
result_utr_hail1 <- soybean_utr_optsolver_hail1(optim_params_conversion(optim_params_short_SoyFACE))
result_partitioning__hail1 <- with(soybean, run_biocro(initial_values,
                                                        parameters,
                                                       weather.growingseason_1,
                                                        direct_modules,
                                                        differential_modules,
                                                        solver,
                                                        verbose = FALSE))
# update parameters
differential_quantities_just_before_defoliation <-
  as.list(result_utr_hail1[nrow(result_utr_hail1), names(initial_values)])

differential_quantities_just_after_defoliation <-
  differential_quantities_just_before_defoliation
# Reduce the leaf mass
remaining_leaf_percent <- 0.4
# differential_quantities_just_after_defoliation$Leaf <-
#   differential_quantities_just_before_defoliation$Leaf * remaining_leaf_percent

differential_quantities_just_after_defoliation$Leaf_substrate_carbon <-
  differential_quantities_just_before_defoliation$Leaf_substrate_carbon * remaining_leaf_percent
differential_quantities_just_after_defoliation$Leaf_structural_carbon <-
  differential_quantities_just_before_defoliation$Leaf_structural_carbon * remaining_leaf_percent # Could be changed to a different percentage

# leaf_C_before_defoliation <- differential_quantities_just_before_defoliation$Leaf_substrate_carbon+
#   differential_quantities_just_before_defoliation$Leaf_structural_carbon
# stem_C_before_defoliation <- differential_quantities_just_before_defoliation$Stem_substrate_carbon+
#   differential_quantities_just_before_defoliation$Stem_structural_carbon
# abg_C_before_defoliation <- leaf_C_before_defoliation + stem_C_before_defoliation
# 
# leaf_C_after_defoliation <- differential_quantities_just_after_defoliation$Leaf_substrate_carbon+
#   differential_quantities_just_after_defoliation$Leaf_structural_carbon
# 
# remaining_stem_percent <- (abg_C_before_defoliation*(1-0.21) -
#                           leaf_C_after_defoliation)/
#   stem_C_before_defoliation

remaining_stem_percent <- 0.5


differential_quantities_just_after_defoliation$Stem_substrate_carbon <-
    differential_quantities_just_before_defoliation$Stem_substrate_carbon * remaining_stem_percent
differential_quantities_just_after_defoliation$Stem_structural_carbon <-
    differential_quantities_just_before_defoliation$Stem_structural_carbon * remaining_stem_percent

differential_quantities_just_after_defoliation$DVI <-
  differential_quantities_just_before_defoliation$DVI - 0.2

orignal_utr_params <- data.frame(optim_params_conversion(optim_params_short_SoyFACE))
rownames(orignal_utr_params) <- arg_names
colnames(orignal_utr_params) <- "Value"
parameters_after_hail <- orignal_utr_params

CHANGE_PARAMETERS <- FALSE
if(CHANGE_PARAMETERS){
  parameters_after_hail['Stem_respiration_factor', 'Value'] <- orignal_utr_params['Stem_respiration_factor', 'Value'] * 2
  # parameters_after_hail['Pod_start_dvi', 'Value'] <- orignal_utr_params['Pod_start_dvi', 'Value'] + 0.1
}

soybean_utr_optsolver_hail2 <- partial_run_biocro(differential_quantities_just_after_defoliation,
                                          parameters,
                                          weather.growingseason_2,
                                          direct_modules,
                                          differential_modules,
                                          solver,
                                          arg_names,
                                          verbose = TRUE)
result_utr_hail2 <- soybean_utr_optsolver_hail2(parameters_after_hail$Value)
# Get the final values of the differential quantities; these will be the
# values just before defoliation
differential_quantities_just_before_defoliation <-
  as.list(result_partitioning__hail1[nrow(result_partitioning__hail1), names(soybean$initial_values)])

# Now reduce the leaf mass
differential_quantities_just_after_defoliation <-
  differential_quantities_just_before_defoliation

differential_quantities_just_after_defoliation[['Leaf']] <-
  differential_quantities_just_before_defoliation[['Leaf']] * remaining_leaf_percent

differential_quantities_just_after_defoliation[['Stem']] <-
  differential_quantities_just_before_defoliation[['Stem']] * remaining_stem_percent

differential_quantities_just_after_defoliation[['DVI']] <-
  differential_quantities_just_before_defoliation[['DVI']] - 0.2

soybean$parameters$mrc_stem <- soybean$parameters$mrc_stem *2 
result_partitioning__hail2 <- run_biocro(differential_quantities_just_after_defoliation,
                                         soybean$parameters,
                                         weather.growingseason_2,
                                         soybean$direct_modules,
                                         soybean$differential_modules,
                                         soybean$ode_solver,
                                         verbose = FALSE)

result_utr_hail <- rbind(result_utr_hail1[seq_len(nrow(result_utr_hail1) - 1), ], result_utr_hail2)
result_partitioning_hail <- rbind(result_partitioning__hail1[seq_len(nrow(result_partitioning__hail2) - 1), ], result_partitioning__hail2)
result_partitioning_hail$Pod <- result_partitioning_hail$Grain + result_partitioning_hail$Shell
fig_2003_utr_with_hail <- plot_SoyFACE_biomass(result_utr_hail, ExpBiomass, ExpBiomass.std, co2_opt, yr)
fig_2003_partitioning_with_hail <- plot_SoyFACE_biomass(result_partitioning_hail, ExpBiomass, ExpBiomass.std, co2_opt, yr)


### elevated CO2
co2_opt <- '_CO2_'
parameters$Catm <- 550
soybean$parameters$Catm <- 550

ExpBiomass.elevCO2 <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt , 'biomass.csv'))
colnames(ExpBiomass.elevCO2)<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

ExpBiomass.elevCO2.std <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
colnames(ExpBiomass.elevCO2.std)<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

RootVals<- data.frame("DOY"=ExpBiomass.elevCO2$DOY[3], "Root"=0.17*sum(ExpBiomass.elevCO2[5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130

numrows<- nrow(weather.growingseason)
invwts <- ExpBiomass.elevCO2.std

soybean_utr_optsolver_no_hail_eCO2 <- partial_run_biocro(initial_values,
                                                    parameters,
                                                    weather.growingseason,
                                                    direct_modules,
                                                    differential_modules,
                                                    solver,
                                                    arg_names,
                                                    verbose = FALSE)

result_utr_no_hail_eCO2 <- soybean_utr_optsolver_no_hail_eCO2(optim_params_conversion(optim_params_short_SoyFACE)) 

result_partitioning_no_hail_eCO2 <- with(soybean, run_biocro(initial_values,
                                                        parameters,
                                                        weather.growingseason,
                                                        direct_modules,
                                                        differential_modules,
                                                        solver,
                                                        verbose = FALSE))

# plot simulation vs observation biomass and calculate mse
fig_2003_utr_no_hail_eCO2 <- plot_SoyFACE_biomass(result_utr_no_hail_eCO2, ExpBiomass.elevCO2, ExpBiomass.elevCO2.std, co2_opt, yr)
result_partitioning_no_hail_eCO2$Pod <- result_partitioning_no_hail_eCO2$Grain + result_partitioning_no_hail_eCO2$Shell
result_partitioning_no_hail_eCO2 <- plot_SoyFACE_biomass(result_partitioning_no_hail_eCO2, ExpBiomass.elevCO2, ExpBiomass.elevCO2.std, co2_opt, yr)
combined_no_hail_eCO2 <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
                                 arrangeGrob(result_partitioning_no_hail_eCO2 + theme(legend.position="none"), top = 'Partitioning Model'),
                                 arrangeGrob(fig_2003_utr_no_hail_eCO2 + theme(legend.position="none"), top = 'UTR model'),
                                 common_legend,
                                 ncol = 4, widths = c(0.1, 5, 5, 1))

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
# update parameters
differential_quantities_just_before_defoliation <-
  as.list(result_utr_eCO2_hail1[nrow(result_utr_eCO2_hail1), names(initial_values)])

differential_quantities_just_after_defoliation <-
  differential_quantities_just_before_defoliation
# Reduce the leaf mass
remaining_leaf_percent <- 0.4
# differential_quantities_just_after_defoliation$Leaf <-
#   differential_quantities_just_before_defoliation$Leaf * remaining_leaf_percent

differential_quantities_just_after_defoliation$Leaf_substrate_carbon <-
  differential_quantities_just_before_defoliation$Leaf_substrate_carbon * remaining_leaf_percent
differential_quantities_just_after_defoliation$Leaf_structural_carbon <-
  differential_quantities_just_before_defoliation$Leaf_structural_carbon * remaining_leaf_percent # Could be changed to a different percentage

# leaf_C_before_defoliation <- differential_quantities_just_before_defoliation$Leaf_substrate_carbon+
#   differential_quantities_just_before_defoliation$Leaf_structural_carbon
# stem_C_before_defoliation <- differential_quantities_just_before_defoliation$Stem_substrate_carbon+
#   differential_quantities_just_before_defoliation$Stem_structural_carbon
# abg_C_before_defoliation <- leaf_C_before_defoliation + stem_C_before_defoliation
# leaf_C_after_defoliation <- differential_quantities_just_after_defoliation$Leaf_substrate_carbon+
#   differential_quantities_just_after_defoliation$Leaf_structural_carbon
# 
# remaining_stem_percent <- (abg_C_before_defoliation*(1-0.21) -
#                              leaf_C_after_defoliation)/
#   stem_C_before_defoliation

remaining_stem_percent <- 0.5


differential_quantities_just_after_defoliation$Stem_substrate_carbon <-
  differential_quantities_just_before_defoliation$Stem_substrate_carbon * remaining_stem_percent
differential_quantities_just_after_defoliation$Stem_structural_carbon <-
  differential_quantities_just_before_defoliation$Stem_structural_carbon * remaining_stem_percent

differential_quantities_just_after_defoliation$DVI <-
  differential_quantities_just_before_defoliation$DVI - 0.2

orignal_utr_params <- data.frame(optim_params_conversion(optim_params_short_SoyFACE))
rownames(orignal_utr_params) <- arg_names
colnames(orignal_utr_params) <- "Value"
parameters_after_hail <- orignal_utr_params

CHANGE_PARAMETERS <- FALSE
if(CHANGE_PARAMETERS){
  parameters_after_hail['Stem_respiration_factor', 'Value'] <- orignal_utr_params['Stem_respiration_factor', 'Value'] * 2
  # parameters_after_hail['Pod_start_dvi', 'Value'] <- orignal_utr_params['Pod_start_dvi', 'Value'] + 0.1
}

soybean_utr_optsolver_eCO2_hail2 <- partial_run_biocro(differential_quantities_just_after_defoliation,
                                                  parameters,
                                                  weather.growingseason_2,
                                                  direct_modules,
                                                  differential_modules,
                                                  solver,
                                                  arg_names,
                                                  verbose = TRUE)
result_utr_eCO2_hail2 <- soybean_utr_optsolver_eCO2_hail2(parameters_after_hail$Value)
# Get the final values of the differential quantities; these will be the
# values just before defoliation
differential_quantities_just_before_defoliation <-
  as.list(result_partitioning_eCO2_hail1[nrow(result_partitioning_eCO2_hail1), names(soybean$initial_values)])

# Now reduce the leaf mass
differential_quantities_just_after_defoliation <-
  differential_quantities_just_before_defoliation

differential_quantities_just_after_defoliation[['Leaf']] <-
  differential_quantities_just_before_defoliation[['Leaf']] * remaining_leaf_percent

differential_quantities_just_after_defoliation[['Stem']] <-
  differential_quantities_just_before_defoliation[['Stem']] * remaining_stem_percent

differential_quantities_just_after_defoliation[['DVI']] <-
  differential_quantities_just_before_defoliation[['DVI']] - 0.2

soybean$parameters$mrc_stem <- soybean$parameters$mrc_stem *2 
result_partitioning_eCO2_hail2 <- run_biocro(differential_quantities_just_after_defoliation,
                                         soybean$parameters,
                                         weather.growingseason_2,
                                         soybean$direct_modules,
                                         soybean$differential_modules,
                                         soybean$ode_solver,
                                         verbose = FALSE)

result_utr_eCO2_hail <- rbind(result_utr_eCO2_hail1[seq_len(nrow(result_utr_eCO2_hail1) - 1), ], result_utr_eCO2_hail2)
result_partitioning_eCO2_hail <- rbind(result_partitioning_eCO2_hail1[seq_len(nrow(result_partitioning_eCO2_hail2) - 1), ], result_partitioning_eCO2_hail2)
result_partitioning_eCO2_hail$Pod <- result_partitioning_eCO2_hail$Grain + result_partitioning_eCO2_hail$Shell
fig_2003_utr_eCO2_with_hail <- plot_SoyFACE_biomass(result_utr_eCO2_hail, ExpBiomass.elevCO2, ExpBiomass.elevCO2.std, co2_opt, yr)
fig_2003_partitioning_eCO2_with_hail <- plot_SoyFACE_biomass(result_partitioning_eCO2_hail, ExpBiomass.elevCO2, ExpBiomass.elevCO2.std, co2_opt, yr)


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
                                       ncol = 2, top = 'Ambient CO2'),
                                     arrangeGrob(
                                       arrangeGrob(fig_2003_partitioning_eCO2_with_hail  + 
                                                     theme(axis.title.x = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           legend.position="none")),
                                       arrangeGrob(fig_2003_utr_eCO2_with_hail  + 
                                                     theme(axis.title.x = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           legend.position="none")), 
                                       ncol = 2, top = 'Elevated CO2'),
                                     arrangeGrob(textGrob('Day of Year (2003)')),
                                     nrow = 3, heights = c(4, 4, 0.3)),
                                 common_legend,
                                 ncol = 3, widths = c(0.3, 5, 1))
ggsave('Fig7-hail.png', 
       plot = combined_with_hail, 
       width = 6,
       height = 4.5,
       units = "in",
       dpi = 600
)

############################################# Plot LD11
years <- c('2021', '2022', '2023', '2024')
Catms <- c(414.7, 417.2, 419.3, 422.8) # from NOAA

updated_parameters <- read.csv('../Data/Soybean-BioCro_Parameters/updated-ld11-parameters.csv')
updated_idx <- match(updated_parameters$New.variable.name, names(parameters))
parameters[updated_idx] <- updated_parameters$LD11

ExpBiomass <- list()
weather.afteremergence <- list()
ld11.results <- list()
ld11.figs <- list()
lai.ld11.figs <- list()
allocation.ld11.figs <- list()

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


for (i in 1:length(years)){
  ExpBiomass[[i]] <- loadRData(paste0('../../energy-farm-biocro/soybean_ld11_biomass_', years[i],'/soybean_ld11_biomass_', years[i], '.RData'))
  sow.time <- min(ExpBiomass[[i]]$time)
  harv.time <- max(ExpBiomass[[i]]$time)
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
  cf <- 0.3 
  updated_utr_initial_values <- list(
    Leaf_substrate_carbon = sub_frac * seed_mass * leaf_frac / cf,
    Leaf_structural_carbon = str_frac * seed_mass * leaf_frac / cf,
    Stem_substrate_carbon = sub_frac * seed_mass * stem_frac / cf,
    Stem_structural_carbon = str_frac * seed_mass * stem_frac / cf,
    Root_substrate_carbon =  sub_frac * seed_mass * root_frac / cf,
    Root_structural_carbon = str_frac * seed_mass * root_frac / cf)
  
  initial_values[names(updated_utr_initial_values)] <- updated_utr_initial_values
  # weather file
  weather.aftersowing <- weather[(((weather$doy-1)*24 + weather$hour) >= sow.time) & ((weather$doy-1)*24 + weather$hour) <= harv.time, ]
  weather.aftersowing$time_zone_offset <- NULL
  parameters$RL_at_25 <- 1.28
  # parameters$par_energy_content <- 0.235
  result <- run_biocro(
    initial_values,
    parameters,
    weather.aftersowing, # weather.afteremergence[[i]],
    direct_modules,#[-1],
    differential_modules,
    solver)
  
  ld11.results[[i]] <- result
  
  # plot LD11 biomass
  ld11.figs[[i]] <- plot_ld11_biomass(result, ExpBiomass[[i]])
  
  allocation.ld11.figs[[i]] <- plot_partitioning(result, years[i])
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
       plot = combined_graph.allocation, 
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

total_precip <- sapply(1:4, function(i) {
  ld11.results[[i]] %>%
    filter(DVI >= 0 & DVI <= 1) %>%
    summarise(total = sum(precip, na.rm = TRUE)) %>%
    pull(total)
})

names(total_precip) <- years
print(total_precip)


# Plot substrate C concentration
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
r <- ld11.results[[2]]
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
# leaf_dirunal_plot <- ggplot(leaf.tnc.sampling.days, aes(hour, TNC, group = Source)) + 
#   geom_point(aes(shape=Source, color=Source, size=Source))+
#   facet_wrap("DOY_date") +
#   scale_shape_manual(values=c(17, 18, 16)) +
#   scale_size_manual(values=c(1, 2.5, 0.5)) +
#   theme_classic() +
#   theme(legend.position = c(0.85, 0.2),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent", colour = NA))+
#   # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
#   scale_x_continuous(breaks = seq(0,24,6))+
#   labs(title=element_blank(), 
#        x='Hour of the Day',
#        y='Leaf Substrate C (mol glucose eq./Mg)')
# stem_dirunal_plot <- ggplot(stem.tnc.sampling.days, aes(hour, TNC, group = Source)) + 
#   geom_point(aes(shape=Source, color=Source, size=Source))+
#   facet_wrap("DOY_date") +
#   scale_shape_manual(values=c(17, 18, 16)) +
#   scale_size_manual(values=c(1, 2.5, 0.5)) +
#   theme_classic() +
#   theme(legend.position = c(0.85, 0.2),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent", colour = NA))+
#   # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
#   scale_x_continuous(breaks = seq(0,24,6))+
#   labs(title=element_blank(), 
#        x='Hour',
#        y='Stem Substrate C (mol glucose eq./Mg)')
# # common_legend <- g_legend(stem_dirunal_plot)
# combined_graph.diurnal_C_concentration <- grid.arrange(arrangeGrob(textGrob('Substrate C Concentration (mol C / kg)', rot = 90, gp=gpar(fontsize=12))),
#                                           arrangeGrob(arrangeGrob(leaf_dirunal_plot, top = '(A) Leaf'),
#                                                       arrangeGrob(stem_dirunal_plot + theme(legend.position="none"), top = '(B) Stem')),
#                                           ncol = 2, widths=c(0.1, 10))



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


# plot assimilation by layers
times <- c(186.5, 209.5, 236.5, 258.5)
doys <- c(186, 209, 236, 258)
layer_assim <- data.frame(DOY = numeric(),
                          layer_number = numeric(),
                          layer_assimilation = numeric())# ,
                          # assim_type = character())
for (t in 1:length(times)){
  time <- times[t]
  doy <- doys[t]
  idx <- which(r$fractional_doy==time)
  for (i in 0:9){
    total_assim <- (r[idx, paste0('sunlit_Assim_layer_', i)]*
                r[idx, paste0('sunlit_fraction_layer_', i)] +
                r[idx, paste0('shaded_Assim_layer_', i)]*
                r[idx, paste0('shaded_fraction_layer_', i)]) *
      r[idx, 'lai'] /10
    # sunlit_assim <- r[idx, paste0('sunlit_Assim_layer_', i)]#  + 
    # shaded_assim <- r[idx, paste0('shaded_Assim_layer_', i)]
    # total_assim <- sunlit_assim + shaded_assim
    # new_row_sunlit <- data.frame(
    #   DOY = doy,
    #   layer_number = i,
    #   layer_assimilation = sunlit_assim,
    #   assim_type = 'sunlit')
    # new_row_shaded <- data.frame(
    #   DOY = doy,
    #   layer_number = i,
    #   layer_assimilation = shaded_assim,
    #   assim_type = 'shaded')
    # new_row <- rbind(new_row_sunlit, new_row_shaded)
    new_row <- data.frame(
      DOY = doy,
      layer_number = i,
      layer_assimilation = total_assim)
    layer_assim <- rbind(layer_assim, new_row)
  }
}

layer_assim_plot <- ggplot(layer_assim, aes(x = layer_number, y = layer_assimilation))+ #, group = assim_type)) +
  geom_point() + # aes(color=assim_type)) +
  geom_line() + #aes(color=assim_type)) +
  facet_wrap(~ DOY, nrow = 1, scales = "fixed") +
  labs(# title = "Layer Assimilation by Layer Number on Different DOYs",
    x = "Layer Number (Top:1 - Bottom: 9)",
    y = "Layer Assimilation (micromol / m^2 / s )") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(color = "grey90")
  ) +
  scale_x_continuous(breaks = unique(layer_assim$layer_number),
                     labels = as.integer(unique(layer_assim$layer_number))) 

# Display the plot
print(layer_assim_plot)
ggsave('Fig6-layer_assim.png', 
       plot = layer_assim_plot, 
       width = 10,
       height = 2.5,
       units = "in",
       dpi = 600
)

