steady_state_module_names <- c(
  "BioCro:stomata_water_stress_linear", # stomata_water_stress = 
  "BioCro:soybean_development_rate_calculator", 
  "UTRSoybeanBML:thornley_utilization_calculator_lsrp",
  "UTRSoybeanBML:thornley_transport_calculator_lsrp",
  "UTRSoybeanBML:thornley_biomass_calculator_lsrp",
  "UTRSoybeanBML:lai_from_structural_carbon", # Originally  "BioCro:parameter_calculator", 
  "BioCro:soil_evaporation",
  "BioCro:solar_position_michalsky",# solar_coordinates = 
  "BioCro:shortwave_atmospheric_scattering",
  "BioCro:incident_shortwave_from_ground_par",
  "BioCro:ten_layer_canopy_properties",
  "BioCro:ten_layer_c3_canopy", # canopy_photosynthesis = 
  "BioCro:ten_layer_canopy_integrator"
  # leaf_water_stress = "BioCro:leaf_water_stress_exponential", (included in the most updated Soybean-BioCro)
)

derivative_module_names <- c(
  "UTRSoybeanBML:thornley_utilization_lsrp",
  "UTRSoybeanBML:thornley_transport_lsrp",
  "BioCro:two_layer_soil_profile", # soil_profile = 
  # "BioCro:development_index",
  "BioCro:thermal_time_linear" # thermal_time =
)

# Set up basic properties for the solver
solver <- list(
  type = 'boost_rkck54', # boost_rosenbrock/boost_rkck54
  output_step_size = 1.0,
  adaptive_rel_error_tol = 1e-4,
  adaptive_abs_error_tol = 1e-4,
  adaptive_max_steps = 200
)