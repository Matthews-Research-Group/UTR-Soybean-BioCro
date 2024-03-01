library(BioCro)

# Clear workspace
rm(list=ls())

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
  "development_index",
  "two_layer_soil_profile",
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
#  irradiance_direct_fraction = 0.2,
#  irradiance_diffuse_fraction = 0.2 
)

default_carbon_to_mass_factor <- 20
base_utilization_rate_constant <- 5.5e-4
base_utilization_km <- 50
base_conductance <- 3e-3

parameters <- list(
  # Parameters unrelated to any module
  timestep = 1.0,
  
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
  specific_heat = 1010,
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
  
  # incident_shortwave_from_ground_par module
  par_energy_fraction_of_sunlight=          0.5,
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
  # shortwave_atmospheric_scattering module
  atmospheric_pressure  =                  101325,
  atmospheric_transmittance =             0.85,
  atmospheric_scattering  =                0.3
)
# years, sowing dates, and harvesting dates of growing seasons being fit to
year <- c('2002', '2005')
sow.date <- c(152, 148)
harv.date <- c(220, 216)
weather.growingseason <- list()
for (i in 1:length(year)) {
  yr <- year[i]
  weather <- read.csv(file = paste0('Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))
  sd.ind <- which(weather$doy == sow.date[i])[1]
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather.growingseason[[i]] <- weather[sd.ind:hd.ind,]
}
varying_parameters <- weather.growingseason[[1]]

# Test the system construction for errors
system_inputs_are_valid <- validate_system_inputs(
  initial_state,
  parameters,
  varying_parameters,
  steady_state_module_names,
  derivative_module_names,
  silent = FALSE
)

# Stop if there was a problem
if (!system_inputs_are_valid) {
  stop("The system inputs were not valid.")
}

# Set up basic properties for the solver
solver <- list(
  type = 'Gro_rsnbrk',
  output_step_size = 1.0,
  adaptive_rel_error_tol = 1e-4,
  adaptive_abs_error_tole = 1e-4,
  adaptive_max_steps = 200
)

# set working directory to location of this file
#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)


# names of parameters being fit
arg_names <- c('Catm') # atmospheric CO2 parameter
params.ambient <- c(400) # ambient atmospheric CO2


# years, sowing dates, and harvesting dates of growing seasons being fit to
year <- c('2002', '2005')
sow.date <- c(152, 148)
harv.date <- c(220, 224)

## Initialize variables for the cost function
soybean_optsolver <- list()

for (i in 1:length(year)) {
  yr <- year[i]
  weather <- read.csv(file = paste0('Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))
  sd.ind <- which(weather$doy == sow.date[i])[1]
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather.growingseason <- weather[sd.ind:hd.ind,]
  
  soybean_optsolver[[i]] <- partial_gro_solver(initial_state,
                                               parameters,
                                               weather.growingseason,
                                               steady_state_module_names,
                                               derivative_module_names,
                                               arg_names,
                                               solver,
                                               verbose = FALSE)
}
result <-  soybean_optsolver[[1]](params.ambient)
