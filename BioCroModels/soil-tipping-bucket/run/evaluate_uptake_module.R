library(BioCro)
viewlibrary(BioCroWater)
library(lattice)

result <- with(soybean, {run_biocro(
  initial_values,
  parameters,
  soybean_weather$'2002',
  direct_modules,
  differential_modules,
  ode_solver
)})

xyplot(Stem + Leaf ~ TTc, data = result, type='l', auto = TRUE)

View(BioCro::soybean$canopy_transpiration_rate)
View(canopy_transpiration_rate)
View(BioCro::soybean$results.canopy_transpiration_rate)
View(BioCro::soybean$result.canopy_transpiration_rate)

library(BioCroWater)

weatherData <- weather$`2000` 
weather <- get_growing_season_climate(weatherData, threshold_temperature = 0)
soiltype = 3
direct_modules_new = list("BioCroWater:soil_type_selector", "BioCroWater:soil_surface_runoff",
                          "BioCroWater:soil_water_downflow", "BioCroWater:soil_water_tiledrain", 
                          "BioCroWater:soil_water_upflow","BioCroWater:multilayer_soil_profile_avg")

differential_modules_new = list("BioCroWater:soil_evaporation","BioCroWater:multi_layer_soil_profile")

init_values =   within(miscanthus_x_giganteus$initial_values, {
  soil_depth_1 = 5
  soil_depth_2 = 15
  soil_depth_3 = 30
  soil_depth_4 = 50
  soil_depth_5 = 50
  soil_depth_6 = 50
  soil_water_content_1 = 0.34
  soil_water_content_2 = 0.34
  soil_water_content_3 = 0.30
  soil_water_content_4 = 0.25
  soil_water_content_5 = 0.2
  soil_water_content_6 = 0.13
  max_rooting_layer = 3
  canopy_transpiration_rate = 0.461
})

# setwd("C:/Users/natal")
evaluate_module('BioCroWater:soil_water_uptake', init_values)

                
