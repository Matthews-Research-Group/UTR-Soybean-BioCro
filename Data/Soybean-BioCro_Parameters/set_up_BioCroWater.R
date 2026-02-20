set_direct_modules<-function(current_direct_modules){
  direct_modules_soil_water = list("BioCroWater:soil_type_selector", "BioCroWater:soil_surface_runoff",
                                   "BioCroWater:soil_water_downflow", "BioCroWater:soil_water_tiledrain", 
                                   "BioCroWater:soil_water_upflow","BioCroWater:soil_water_uptake",
                                   "BioCroWater:multilayer_soil_profile_avg")
  direct_modules_new = current_direct_modules
  old_soil_evapo_index = which(direct_modules_new=="BioCro:soil_evaporation")
  direct_modules_new = direct_modules_new[-old_soil_evapo_index] #remove soil evaporation
  direct_modules_new = c(direct_modules_new[1:(old_soil_evapo_index-1)],direct_modules_soil_water,
                         direct_modules_new[old_soil_evapo_index:length(direct_modules_new)]) # insert BioCroWater
  direct_modules_new$leaf_water_stress = NULL
  return(direct_modules_new)
}

set_differential_modules<-function(current_differential_modules){
  differential_modules_soil_water = list("BioCroWater:soil_evaporation","BioCroWater:multi_layer_soil_profile")
  differential_modules_new = current_differential_modules
  old_soil_profile_index = which(differential_modules_new=="BioCro:two_layer_soil_profile")
  differential_modules_new = differential_modules_new[-old_soil_profile_index] #remove soil_profile
  differential_modules_new = c(differential_modules_new[1:(old_soil_profile_index-1)],
                               differential_modules_soil_water,
                               differential_modules_new[old_soil_profile_index:length(differential_modules_new)])
  return(differential_modules_new)
}

set_init_values<-function(current_initial_values){
  init_values =   within(current_initial_values,{
    soil_water_content_1 = current_initial_values$soil_water_content  #0-5cm
    soil_water_content_2 = current_initial_values$soil_water_content  #5-15cm
    soil_water_content_3 = current_initial_values$soil_water_content  #15-35cm
    soil_water_content_4 = current_initial_values$soil_water_content  #35-55 cm
    soil_water_content_5 = current_initial_values$soil_water_content  #55-75 cm
    soil_water_content_6 = current_initial_values$soil_water_content  #75-100cm
    sumes1               = 0.125
    sumes2               = 0
    time_factor = 0
    soil_evaporation_rate = 0
  })
  init_values$soil_water_content=NULL
  return(init_values)
}

set_parameters<-function(current_parameters){
  soiltype = 4 
  parameters =   within(current_parameters, {
    irrigation = 0.0
    swcon = 0.05/24
    curve_number = 61
    soil_type_indicator_1 = 1 
    soil_type_indicator_2 = soiltype
    soil_type_indicator_3 = soiltype
    soil_type_indicator_4 = soiltype
    soil_type_indicator_5 = soiltype
    soil_type_indicator_6 = soiltype
    soil_depth_1 = 5
    soil_depth_2 = 10
    soil_depth_3 = 20
    soil_depth_4 = 20
    soil_depth_5 = 20
    soil_depth_6 = 25
    tile_drainage_rate          = 0.2 # dimensionless, as in DSSAT (1/d ?)
    tile_drain_depth            = 90  # cm (typical value 3-4 ft)
    skc    = 0.3
    kcbmax = 0.25
    elevation                   = 219
    bare_soil_albedo            = 0.15
    max_rooting_layer = 4
  })
  parameters$LeafWS = 1
  parameters[c('soil_field_capacity','soil_saturated_conductivity','soil_saturation_capacity','soil_wilting_point')]=NULL
  parameters[c('kShell','net_assimilation_rate_shell')] = NULL
  return(parameters)
}
