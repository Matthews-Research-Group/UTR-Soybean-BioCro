library(BioCro)
library(BioCroWater)
weatherData <- weather$`2000` 
weather <- get_growing_season_climate(weatherData, threshold_temperature = 0)
soiltype = 3
direct_modules_new = list("BioCroWater:soil_type_selector", "BioCroWater:soil_surface_runoff",
                          "BioCroWater:soil_water_downflow", "BioCroWater:soil_water_tiledrain", 
                          "BioCroWater:soil_water_upflow","BioCroWater:multilayer_soil_profile_avg")

differential_modules_new = list("BioCroWater:soil_evaporation","BioCroWater:multi_layer_soil_profile")

init_values =   within(miscanthus_x_giganteus$initial_values,{
  soil_water_content_1 = 0.34
  soil_water_content_2 = 0.34
  soil_water_content_3 = 0.30
  soil_water_content_4 = 0.25
  soil_water_content_5 = 0.2
  soil_water_content_6 = 0.13
  sumes1               = 0.125
  sumes2               = 0
  time_factor = 0
  soil_evaporation_rate = 0
})
init_values$soil_water_content=NULL
parameters =   within(miscanthus_x_giganteus$parameters, {
  irrigation = 0.0
  swcon = 0.05/24
  curve_number = 61
  soil_type_indicator_1 = soiltype
  soil_type_indicator_2 = soiltype
  soil_type_indicator_3 = soiltype
  soil_type_indicator_4 = soiltype
  soil_type_indicator_5 = soiltype
  soil_type_indicator_6 = soiltype
  soil_depth_1 = 5
  soil_depth_2 = 15
  soil_depth_3 = 30
  soil_depth_4 = 50
  soil_depth_5 = 50
  soil_depth_6 = 50
  tile_drainage_rate          = 0.2 # dimensionless, as in DSSAT (1/d ?)
  tile_drain_depth            = 90  # cm (typical value 3-4 ft)
  skc = 0.3
  kcbmax = 0.25
  elevation                   = 219
  lai = 2
  bare_soil_albedo            = 0.15
})
parameters[c('soil_field_capacity','soil_saturated_conductivity','soil_saturation_capacity','soil_wilting_point')]=NULL

result <- run_biocro(
  init_values,
  parameters,
  weather,
  direct_modules_new,
  differential_modules_new
)

plot(result$time, result$soil_water_content, xlab='time',ylab='soil water content')
