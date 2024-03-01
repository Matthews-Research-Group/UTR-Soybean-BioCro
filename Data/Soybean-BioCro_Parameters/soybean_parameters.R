# Do the calculations inside an empty list so that temporary variables are not created in .Global.
# parameters, and weather data
parameters <- list(
  # Parameters unrelated to any module
  timestep = 1.0,
  
  default_carbon_to_mass_factor <- 0.33,
  base_utilization_rate_constant <- 3e-2,
  base_utilization_km <- 1,
  base_conductance <- 0.5,
  # Parameters related to the Thornley growth model
  Leaf_carbon_to_mass_factor = default_carbon_to_mass_factor,
  Leaf_utilization_rate_constant = 0.4 * base_utilization_rate_constant,
  Leaf_utilization_km = 1.0 * base_utilization_km,
  Leaf_respiration_factor = 0.0,				# Leaf respiration is accounted for
                                        # by the canopy photosynthesis module
  Stem_carbon_to_mass_factor = default_carbon_to_mass_factor,
  Stem_utilization_rate_constant = 1.0 * base_utilization_rate_constant,
  Stem_utilization_km = 1.2 * base_utilization_km,
  Stem_respiration_factor = 0.02,
  
  Root_carbon_to_mass_factor = default_carbon_to_mass_factor,
  Root_utilization_rate_constant = 0.1 * base_utilization_rate_constant,
  Root_utilization_km = 0.5 * base_utilization_km,
  Root_respiration_factor = 0.03,
  
  Pod_carbon_to_mass_factor = default_carbon_to_mass_factor,
  Pod_utilization_rate_constant = 1.0 * base_utilization_rate_constant,
  Pod_utilization_km = 1.0 * base_utilization_km,
  Pod_respiration_factor = 0.02,
  
  substrate_conductance_Leaf_to_Stem = base_conductance,
  substrate_conductance_Stem_to_Root = 0.4 * base_conductance,
  substrate_conductance_Stem_to_Pod = 2*base_conductance,
  
  transportation_beta_exponent = 1,
  Pod_start_dvi = 1.0,
  stop_growth_dvi = 2.0,
  
  # senescence_coefficient_logistic module
  Leaf_senescence_rate_max =               0.003,
  Stem_senescence_rate_max =               0.003,
  Root_senescence_rate_max =               0.005,  # senescence of root not simulated in Soybean-BioCro (see Note 1 at end of file)
  Pod_senescence_rate_max  =               0.0,
  Leaf_senescence_alpha    =               20.0,   
  Stem_senescence_alpha    =               20.0,
  Root_senescence_alpha    =               8.0,    # 3 for 2002,2005; 4 for 2004; 3 for 2006 (not sure)
  Pod_senescence_alpha     =               110.0,
  Leaf_senescence_beta     =               -0.1,
  Stem_senescence_beta     =               -0.1,
  Root_senescence_beta     =               -2.0,
  Pod_senescence_beta      =               -0.4,
  Leaf_senescence_reuse_factor =            0.6,
  Stem_senescence_reuse_factor =            0.6,
  Root_senescence_reuse_factor =            0.6,
  Pod_senescence_reuse_factor =             0.6,
  
  # Parameters related to the `parameter_calculator` module
  iSp = 2.5,
  Sp_thermal_time_decay = 0.0,
  LeafN_0 = 2.0,
  LeafN = 2.0,
  vmax_n_intercept = 0,
  vmax1 = 111.2,
  alphab1 = 0,
  alpha1 = 32.5,
  
  # soil parameters (clay loam)
  soil_air_entry              = -2.6,
  soil_b_coefficient          = 5.2,
  soil_bulk_density           = 1.35,
  soil_clay_content           = 0.34,
  soil_field_capacity         = 0.32,
  soil_sand_content           = 0.32,
  soil_saturated_conductivity = 6.4e-05,
  soil_saturation_capacity    = 0.52,
  soil_silt_content           = 0.34,
  soil_wilting_point          = 0.2,
  
  # Parameters related to the `soil_evaporation` module
  rsec = 0.2,
  soil_clod_size = 0.04,
  soil_reflectance = 0.2,
  soil_transmission = 0.01,
  specific_heat_of_air=1010,
  #  specific_heat = 1010,
  stefan_boltzman = 5.67e-8,
  
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
  
  # solar_position_michalsky module
  lat                         =            40,
  longitude                   =            -88,
  time_zone_offset            =            -6,
  
  # shortwave_atmospheric_scattering module
  atmospheric_pressure  =                  101325,
  atmospheric_transmittance =              0.85,
  atmospheric_scattering  =                0.3,
  
  # ten_layer_canopy_properties module
  absorptivity_par            = 0.8,         # Campbell and Norman, An Introduction to Environmental Biophysics, 2nd Edition
  chil                        = 0.81,        # Campbell and Norman, An Introduction to Environmental Biophysics, 2nd Edition, Table 15.1, pg 253
  kd                          = 0.7,         # Estimated from Campbell and Norman, An Introduction to Environmental Biophysics, 2nd Edition, Figure 15.4, pg 254
  heightf                     = 3,           # m^-1
  kpLN                        = 0,           # not used in Soybean-BioCro
  leaf_reflectance            = 0.2,
  leaf_transmittance          = 0.2,
  lnfun                       = 0,           # not used in Soybean-BioCro
  
  # ten_layer_c3_canopy module
  jmax                        = 195,         # Bernacchi et al. 2005 (https://doi.org/10.1007/s00425-004-1320-8), 2002 Seasonal average
  jmax_mature                 = 195,         # Needed in the varying_Jmax25 module
  sf_jmax                     = 0.2,         # Scaling factor for jmax. Needed in the varying_Jmax25 module
  electrons_per_carboxylation = 4.5,         # Bernacchi et al. 2003 (https://doi.org/10.1046/j.0016-8025.2003.01050.x)
  electrons_per_oxygenation   = 5.25,        # Bernacchi et al. 2003 (https://doi.org/10.1046/j.0016-8025.2003.01050.x)
  tpu_rate_max                = 13,          # Fitted value based on the A-Ci data measured at UIUC in 2019-08 by Delgrado (unpublished data)
  Rd                          = 1.28,        # Davey et al. 2004 (https://doi.org/10.1104/pp.103.030569), Table 3, cv Pana, co2 368 ppm
  Catm                        = 372.59,      # micromol / mol, CO2 level in 2002
  O2                          = 210,         # millimol / mol
  b0                          = 0.008,       # Leakey et al. 2006 (https://10.1111/j.1365-3040.2006.01556.x)
  b1                          = 10.6,        # Leakey et al. 2006 (https://10.1111/j.1365-3040.2006.01556.x)
  Gs_min                      = 1e-3,
  theta                       = 0.76,        # Bernacchi et al. 2003 (https://doi.org/10.1046/j.0016-8025.2003.01050.x)
  minimum_gbw                 = 0.08,
  windspeed_height            = 5,
  beta_PSII                   = 0.5,         # Bernacchi et al. 2003 (https://doi.org/10.1046/j.0016-8025.2003.01050.x)
  
  # ten_layer_canopy_integrator module
  growth_respiration_fraction = 0
)

if (co2_opt == '_ambient_'){
  parameters$Catm <- 372 # 414.71ppm for 2021 from NOAA
}else{
  parameters$Catm <- 550
}