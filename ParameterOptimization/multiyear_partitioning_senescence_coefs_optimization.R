# The work directory is hard coded since the previous one doesn't work well all the time. 
# Please change it according to your own working directory.
# BioCro version: https://github.com/pxm726/Soybean-Thornley-BioCro/tree/main/biocro-dev
# Parameter loading and plotting are put into this R code to check the variables more easily.
# To be moved to different files once everything works out.
library(BioCro)
library(UTRSoybeanBML)
library(DEoptim)
library(lattice)
library(lhs)
library(ppcor)

# Clear workspace
rm(list=ls())

save_to_pdf <- FALSE

# Set working directory to location of this file
# this.dir <- dirname(parent.frame(2)$ofile)
# setwd("C://Users//pxm72//Documents//Projects//Soybean-BioCro//ParameterOptimization")

# Cost function
source('multiyear_BioCro_optim.R')

# load parameter files
# setwd('../Data/Soybean-BioCro_Parameters')
# lapply(list.files(pattern = "[.]R$", recursive = TRUE), source)

# parameters, and weather data
steady_state_module_names <- c(
  "soil_type_selector",
  "stomata_water_stress_linear",
  "soybean_development_rate_calculator", # Soybean specific
  "thornley_utilization_calculator_lsr",
  "thornley_transport_calculator_lsr",
  "parameter_calculator",
  "soil_evaporation",
  "solar_zenith_angle", # Added
  "shortwave_atmospheric_scattering",
  "incident_shortwave_from_ground_par",
  "ten_layer_canopy_properties",
  "ten_layer_c3_canopy",
  "ten_layer_canopy_integrator"
)

derivative_module_names <- c(
  "thornley_utilization_lsr",
  "thornley_transport_lsr",
  #  "thermal_time_accumulator",
  #  "one_layer_soil_profile",
  "two_layer_soil_profile",
  "development_index",
  "thermal_time_linear"
)

initial_state <- list(
  # Variables related to the utilization growth model
  Leaf_substrate_carbon = 3e-8,
  Leaf_structural_carbon = 6e-8,
  Leaf_respiration_loss = 0.0,
  Stem_substrate_carbon = 0.5e-8,
  Stem_structural_carbon = 1e-8,
  Stem_respiration_loss = 0.0,
  Root_substrate_carbon = 0.5e-8,
  Root_structural_carbon = 1e-8,
  Root_respiration_loss = 0.0,
  
  # Other variables
  TTc = 0.0,
  soil_water_content = 0.32,
  
  cws1=                     0.32,          # dimensionless, current water status, soil layer 1
  cws2=                     0.32,          # dimensionless, current water status, soil layer 2
  DVI =                     -1             # Sowing date: DVI=-1
)



parameters <- list(
  # Parameters unrelated to any module
  timestep = 1.0,
  
  default_carbon_to_mass_factor <- 0.34,
  base_utilization_rate_constant <- 5.5e-4,
  base_utilization_km <- 50,
  base_conductance <- 3e-3,
  # Parameters related to the Thornley growth model
  
  Leaf_structural_carbon_to_mass_factor = default_carbon_to_mass_factor,
  Leaf_utilization_rate_constant = 0.4 * base_utilization_rate_constant,
  Leaf_utilization_km = 1.0 * base_utilization_km,
  Leaf_respiration_factor = 0.0,				# Leaf respiration is accounted for
  # by the canopy photosynthesis module
  
  Stem_structural_carbon_to_mass_factor = default_carbon_to_mass_factor,
  Stem_utilization_rate_constant = 1.0 * base_utilization_rate_constant,
  Stem_utilization_km = 1.2 * base_utilization_km,
  Stem_respiration_factor = 0.02,
  
  Root_structural_carbon_to_mass_factor = default_carbon_to_mass_factor,
  Root_utilization_rate_constant = 0.1 * base_utilization_rate_constant,
  Root_utilization_km = 0.5 * base_utilization_km,
  Root_respiration_factor = 0.03,
  
  substrate_conductance_Leaf_to_Stem = base_conductance,
  substrate_conductance_Stem_to_Root = 0.4 * base_conductance,
  
  
  utilization_beta_exponent = 0.333,
  
  # Parameters related to the `parameter_calculator` module
  iSp = 2.5,
  Sp_thermal_time_decay = 0.0,
  LeafN_0 = 2.0,
  LeafN = 2.0,
  vmax_n_intercept = 0,
  vmax1 = 111.2,
  alphab1 = 0,
  alpha1 = 32.5,
  
  # Parameters related to the `soil_type_selector` module
  soil_type_indicator = 6,
  
  # Parameters related to the `soil_evaporation` module
  rsec = 0.2,
  soil_clod_size = 0.04,
  soil_reflectance = 0.2,
  soil_transmission = 0.01,
  specific_heat_of_air=1010,
#  specific_heat = 1010,
  stefan_boltzman = 5.67e-8,
  
  # Parameters related to the `c3_canopy` module
  lat = 40,
  nlayers = 10,
  jmax = 213.2,
  Rd = 1.1,
  Catm = 370,
  O2 = 210,
  b0 = 0.048,
  b1 = 5,
  theta = 0.7,
  kd = 0.1,
  heightf = 3,
  kpLN = 0.2,
  lnb0 = -5,
  lnb1 = 18,
  lnfun = 0,
  chil = 0.8,
  growth_respiration_fraction = 0,
  water_stress_approach = 1,
  electrons_per_carboxylation = 4.5,
  electrons_per_oxygenation = 5.25,
  
  # Parameters related to the `thermal_time_accumulator` module
  tbase = 10,
  topt_lower = 1e4,
  topt_upper = 1e5,
  tmax = 1e6,
  
  # Parameters related to the `one_layer_soil_profile` module
  soil_depth = 1,
  phi1 = 0.01,
  phi2 = 10,
  
  # soybean_development_rate_calculator module
  maturity_group  =                        3,           # dimensionless; soybean cultivar maturity group
  Tbase_emr       =                        10,          # degrees C
  TTemr_threshold =                        60 ,         # degrees C * days
  Rmax_emrV0      =                        0.1990,      # day^-1; Setiyono et al., 2007 (https://doi.org/10.1016/j.fcr.2006.07.011), Table 2
  Tmin_emrV0      =                        5.0    ,     # degrees C; Setiyono et al., 2007, Table 2
  Topt_emrV0      =                        31.5   ,     # degrees C; Setiyono et al., 2007, Table 2
  Tmax_emrV0      =                        45.0   ,     # degrees C; Setiyono et al., 2007, Table 2
  Tmin_R0R1       =                        5.0    ,     # degrees C; Setiyono et al., 2007, Table 2 (used emrV0 values)
  Topt_R0R1       =                        31.5   ,     # degrees C; Setiyono et al., 2007, Table 2 (used emrV0 values)
  Tmax_R0R1       =                        45.0   ,     # degrees C; Setiyono et al., 2007, Table 2 (used emrV0 values)
  Tmin_R1R7       =                        0.0    ,     # degrees C; Setiyono et al., 2007, Table 2
  Topt_R1R7       =                        21.5   ,     # degrees C; Setiyono et al., 2007, Table 2
  Tmax_R1R7       =                        38.7   ,     # degrees C; Setiyono et al., 2007, Table 2
  sowing_time     =                        0      ,     # Referred to data/soybean.R: Soybean-BioCro uses the weather data to set the sowing time
  
  # incident_shortwave_from_ground_par module
#  par_energy_fraction_of_sunlight=          0.5,
  par_energy_fraction            =        0.5,
  par_energy_content             =        0.235,
  
  # two_layer_soil_profile module
  soil_depth1=                             0.0,         # meters
  soil_depth2=                             2.5,         # meters
  soil_depth3=                             10.0,        # meters
  wsFun     =                              2,           # not used, but must be defined
  hydrDist  =                              0 ,          # same as in sorghum parameter file
  rfl       =                              0.2,         # same as in sorghum parameter file
  rsdf      =                              0.44,        # same as in sorghum parameter file
  phi1      =                              0.01,
  phi2      =                              1.5,          # from Sugarcane-BioCro, Jaiswal et al. 2017 (https://doi.org/10.1038/nclimate3410)
  
  # thermal_time_linear module
  tbase     =                              10,          # degrees C
  # shortwave_atmospheric_scattering module
  atmospheric_pressure  =                  101325,
  atmospheric_transmittance =              0.85,
  atmospheric_scattering  =                0.3,

  # ten_layer_canopy_properties module
  absorptivity_par=                        0.8 ,        # Campbell and Norman, An Introduction to Environmental Biophysics, 2nd Edition
  chil            =                        0.81,        # Campbell and Norman, An Introduction to Environmental Biophysics, 2nd Edition, Table 15.1, pg 253
  kd              =                        0.7,         # Estimated from Campbell and Norman, An Introduction to Environmental Biophysics, 2nd Edition, Figure 15.4, pg 254
  heightf         =                        3  ,         # m^-1
  kpLN            =                        0  ,         # not used in Soybean-BioCro (see Note 1 at end of file)
  lnfun           =                        0,           # not used in Soybean-BioCro
  leaf_transmittance  =        0.2,                     # Referred to data/soybean.R
  leaf_reflectance    =        0.2,                     # Referred to data/soybean.R
  tpu_rate_max    =                        13,          # Fitted value based on the A-Ci data measured at UIUC in 2019-08 by Delgrado (unpublished data)
  Gs_min          =                        1e-3         # Referred to data/soybean.R
)
# years, sowing dates, and harvesting dates of growing seasons being fit to
year <- c('2002', '2005')
sow.date <- c(152, 148)
harv.date <- c(220, 216)
weather.growingseason <- list()
for (i in 1:1){
#for (i in 1:length(year)) {
  yr <- year[i]
  weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))
  sd.ind <- which(weather$doy == sow.date[i])[1]
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather.growingseason[[i]] <- weather[sd.ind:hd.ind,]
}
varying_parameters <- weather.growingseason[[1]]

# # Test the system construction for errors
# system_inputs_are_valid <- validate_system_inputs(
#   initial_state,
#   parameters,
#   varying_parameters,
#   steady_state_module_names,
#   derivative_module_names,
#   silent = TRUE
# )
# 
# # Stop if there was a problem
# if (!system_inputs_are_valid) {
#   stop("The system inputs were not valid.")
# }

# Set up basic properties for the solver
solver <- list(
#  type = 'Gro_rsnbrk',
  type = 'boost_rosenbrock',
#  type = 'Gro_rkck54',
  output_step_size = 1.0,
  adaptive_rel_error_tol = 1e-4,
  adaptive_abs_error_tole = 1e-4,
  adaptive_max_steps = 200
)

# set working directory to location of this file
#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)


# names of parameters being fit
arg_names <- c('Leaf_structural_carbon_to_mass_factor', 'Leaf_utilization_rate_constant',
               'Leaf_utilization_km', 'Leaf_respiration_factor',
               'Stem_structural_carbon_to_mass_factor', 'Stem_utilization_rate_constant',
               'Stem_utilization_km', 'Stem_respiration_factor',
               'Root_structural_carbon_to_mass_factor', 'Root_utilization_rate_constant',
               'Root_utilization_km', 'Root_respiration_factor',
               'substrate_conductance_Leaf_to_Stem', 'substrate_conductance_Stem_to_Root',
               'utilization_beta_exponent')

# years, sowing dates, and harvesting dates of growing seasons being fit to
year <- c('2002', '2005')
sow.date <- c(152, 148)
harv.date <- c(233, 233)
weather.growingseason <- list()
#for (i in 1:length(year)) {   # Remember to change back 
for (i in 1:1) {
  yr <- year[i]
  weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))
  sd.ind <- which(weather$doy == sow.date[i])[1]
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather.growingseason[[i]] <- weather[sd.ind:hd.ind,]
}
varying_parameters <- weather.growingseason[[1]]

# Initialize variables for the cost function
soybean_optsolver <- list()
ExpBiomass <- list()
ExpBiomass.std <- list()
RootVals <- list()
weights <- list()
numrows <- vector()

#for (i in 1:length(year)) {   # Remember to change back 
for (i in 1:1) {
  yr <- year[i]
  weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))
  sd.ind <- which(weather$doy == sow.date[i])[1]
  hd.ind <- which(weather$doy == harv.date[i])[24]

  weather.growingseason <- weather[sd.ind:hd.ind,]

  soybean_optsolver[[i]] <- partial_run_biocro(initial_state,
                                               parameters,
                                               weather.growingseason,
                                               steady_state_module_names,
                                               derivative_module_names,
                                               solver,
                                               arg_names,
                                               verbose = FALSE)

  ExpBiomass[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/',yr,'_ambient_biomass.csv'))[1:5 ,1:3]
  colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem")

  ExpBiomass.std[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/',yr,'_ambient_biomass_std.csv'))[1:5 ,1:3]
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem")

  RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[5], "Root"=0.17*sum(ExpBiomass[[i]][5,2:3])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130

  numrows[i] <- nrow(weather.growingseason)
  invwts <- ExpBiomass.std[[i]]
  weights[[i]] <- log(1/(invwts[,2:ncol(invwts)]+1e-5))

}

wts2 <- data.frame("Stem" = 1, "Leaf" = 1, "Root" = 0.125)

# # parameter upper limit
# upperlim<-c(70,5e-4, 60, 0.1,
#             70,7.5e-4, 72, 0.1,
#             70,2e-4, 40, 0.1,
#             5e-3, 2e-3, 0.5)
# 
# # parameter lower limit
# lowerlim<-c(25, 1e-4, 40, 0,
#             25, 2.5e-4, 48, 0,
#             25, 2.5e-5, 25, 0,
#             1e-3, 5e-4, 0.2)

# parameter upper limit
upperlim <- c(20,5e-3, 100, 0.1,
              5e-3, 5e-2, 0.2)
# parameter lower limit
lowerlim <- c(1, 1e-5, 10, 0,
              1e-3, 1e-5, 0.01)


# cost function
cost_func <- function(x){
  multiyear_BioCro_optim(x, soybean_optsolver, ExpBiomass, numrows, weights, wts2, RootVals)
}

rng.seed <- 123467 # seed for random number generator
set.seed(rng.seed)

# maximum number of iterations
max.iter <- 100

# Call DEoptim function to run optimization
parVars <- c('multiyear_BioCro_optim','soybean_optsolver','ExpBiomass','numrows','weights','wts2','RootVals')
optim_result<-DEoptim(fn=cost_func, lower=lowerlim, upper = upperlim, control=list(itermax=max.iter,parallelType=1,packages=c('BioCro'),parVar=parVars))
stopCluster() #close parallel connections


# print(optim_result)

optim_params7 = optim_result$optim$bestmem
optim_params <- list()
# carbon to mass factor
optim_params[c(1,5,9)] = optim_params7[1]

# utilization rate constant
optim_params[2] = 0.4* optim_params7[2]
optim_params[6] = optim_params7[2] # picked
optim_params[10] = 0.1* optim_params7[2]

# utilization km
optim_params[3] = optim_params7[3] # picked
optim_params[7] = 1.2 * optim_params7[3]
optim_params[11] = 0.5* optim_params7[3]

# respiration factor
optim_params[4] = 0 # Leaf respiration is accounted for by the canopy photosynthesis module
optim_params[8] = optim_params7[4] # picked
optim_params[12] = 1.5* optim_params7[4]

# substrate conductance and utilization beta component
optim_params[c(13,14,15)] = optim_params7[c(5,6,7)]

# optim_params = optim_result$optim$bestmem
result <- soybean_optsolver[[1]](optim_params)

# organize simulated data
r.lsr.doy <- reshape2::melt(result[,c("time","Root","Leaf","Stem")],id.vars="time")
r.lsr.doy$value<-r.lsr.doy$value

# organize the expirimental data (mean and std)
s.exp.leaf <- cbind(ExpBiomass[[1]][,c("DOY","Leaf")])
colnames(s.exp.leaf) <- c("time","Leaf") # DOY renamed as time
r.exp.leaf <- reshape2::melt(s.exp.leaf, id.vars = "time") # DOY renamed as time

s.exp.std.leaf <- cbind(ExpBiomass.std[[1]][,c("DOY","Leaf")])
colnames(s.exp.std.leaf) <- c("time","Leaf") # DOY renamed as time
r.exp.std.leaf <- reshape2::melt(s.exp.std.leaf, id.vars = "time") # DOY renamed as time
r.exp.std.leaf$ymin <- r.exp.leaf$value - r.exp.std.leaf$value
r.exp.std.leaf$ymax <- r.exp.leaf$value + r.exp.std.leaf$value

s.exp.stem <- cbind(ExpBiomass[[1]][,c("DOY","Stem")]) # DOY renamed as time
colnames(s.exp.stem) <- c("time","Stem") # DOY renamed as time
r.exp.stem <- reshape2::melt(s.exp.stem, id.vars = "time") # DOY renamed as time

s.exp.std.stem <- cbind(ExpBiomass.std[[1]][,c("DOY","Stem")]) # DOY renamed as time
colnames(s.exp.std.stem) <- c("time","Stem") # DOY renamed as time
r.exp.std.stem <- reshape2::melt(s.exp.std.stem, id.vars = "time") # DOY renamed as time
r.exp.std.stem$ymin <- r.exp.stem$value - r.exp.std.stem$value
r.exp.std.stem$ymax <- r.exp.stem$value + r.exp.std.stem$value

r.exp.ls <- rbind(r.exp.leaf, r.exp.stem)
r.exp.ls$source = "Observed"
r.lsr.doy$source = "Simulated"
r.all <- rbind(r.lsr.doy, r.exp.ls)
# Reverse the order as follow
r.all$variable <- factor(r.all$variable, levels = rev(levels(r.all$variable)))

# combine the simulated and expirimental data


library(ggplot2)
# Colorblind friendly color palette (https://personal.sron.nl/~pault/)
col.palette.muted <- c("#332288", "#117733", "#999933", "#882255")

size.title <- 16
size.axislabel <-14
size.axis <- 16
size.legend <- 14

f <- ggplot() + theme_classic()

f <- f + geom_point(data=r.all, aes(x=time, y=value, 
                                    color=variable,
                                    size = source, shape = source), 
                    show.legend = TRUE, stroke=0.5) + 
  scale_shape_manual(values = c(15, 16)) + 
  scale_size_manual(values = c(4, 0.5)) +
  scale_color_manual(values = col.palette.muted)

# for leaf
f <- f + geom_errorbar(data=r.exp.std.leaf, aes(x=time, ymin=ymin, ymax=ymax),  # DOY renamed as time
                       width=3.5, size=0.25, show.legend = FALSE)
# for stem
f <- f + geom_errorbar(data=r.exp.std.stem, aes(x=time, ymin=ymin, ymax=ymax),   # DOY renamed as time
                       width=3.5, size=0.25, show.legend = FALSE)

# change the plot labels and theme
f <- f + labs(title=element_blank(), x=paste0('Day of Year (',year,')'),y='Biomass (Mg / ha)')
f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
               axis.text=element_text(size=size.axis),
               axis.title=element_text(size=size.axislabel),
               legend.position = c(.1,.85), legend.title = element_blank(),
               legend.text=element_text(size=size.legend),
               legend.background = element_rect(fill = "transparent",colour = NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA))

f <- f + scale_x_continuous(breaks = seq(150,280,30))
f <- f + guides(shape = guide_legend(override.aes = list(size = 4)))
f <- f + guides(color = guide_legend(override.aes = list(size=4)))
f