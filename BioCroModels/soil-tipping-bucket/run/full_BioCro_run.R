library(BioCro)
library(BioCroWater)
rm(list=ls())
obs_weather<-read.csv('data/doi_10_5061_dryad_g0v62__v20170815/Final_Data_Deposit/ExtendedDataFig1_Weather_Data_2004thru2011/soyFACE_weather_data_2004thru2011.csv')
opt_result_name = "new_water_v128"
use_new_water    = TRUE
use_pre_opt_vars = TRUE
use_predefined   = TRUE
use_varying_sla  = FALSE
# arg_names <- c('alphaLeaf','alphaRoot','alphaStem','betaLeaf','betaRoot','betaStem',
#                'rateSeneLeaf','rateSeneStem','alphaSeneLeaf','betaSeneLeaf',
#                'alphaSeneStem','betaSeneStem','alphaShell','betaShell',
#                'grc_stem','grc_root',
#                'mrc_leaf','mrc_root')
new_par   = readRDS(paste0('optimization/opt_results/opt_result_',opt_result_name,'.rds'))
arg_names = names(new_par)
# new_par = new_par$optim$bestmem
# names(new_par) = arg_names

if(use_pre_opt_vars){
  pre_arg_names <- c('alphaRoot','betaRoot')
  
  pre_par = readRDS(paste0('optimization/opt_results/opt_result_new_water_v108.rds'))
  pre_par = pre_par$optim$bestmem
  pre_par = pre_par[c(2,5)]
}
print(new_par)

# phi2 = 0.5
soiltype = 4
years = c(2002,2004:2011)
obs_data = readRDS("data/soil_moisture_data_rearranged.rds")
obs_data$sm = obs_data$sm/100
#obs depth defition:
#layer 1: 5-15 cm; layer 2: 15-25; layer 3: 25-35.... layer 7: 65-75 cm
results = list()
ctl = list()
for (i in 1:length(years)){
  year = years[i]
  if(year==2002){
    weather_growing_season = soybean_weather$`2002`
    #get obs data
    obs_year_i = obs_data[obs_data$year==2004,]
    obs_avg = list()
    for (depth in 1:7){
      obs = obs_year_i[obs_year_i$depth==depth,] 
      obs = obs[,c('ring','doy','sm')]
      tmp = aggregate(sm ~ doy, data = obs, FUN = mean)
      obs_avg[[depth]] = tmp
    }
  }else{
    obs_weather_yeari = obs_weather[obs_weather$Year==year,]
    weatherData <- 
      read.csv(paste0("WeatherData/",year,"/",year,"_Bondville_IL_daylength.csv"))
    #further subset from the start of obs DOY
    DOY_start = obs_weather_yeari$DOY[1]
    DOY_end = tail(obs_weather_yeari$DOY,1)
    weather_growing_season = weatherData[weatherData$doy>=DOY_start & weatherData$doy<=DOY_end,]
    weather_growing_season$time_zone_offset = -6
    
    #get obs data
    obs_year_i = obs_data[obs_data$year==year,]
    obs_avg = list()
    for (depth in 1:7){
      obs = obs_year_i[obs_year_i$depth==depth,] 
      obs = obs[,c('ring','doy','sm')]
      tmp = aggregate(sm ~ doy, data = obs, FUN = mean)
      obs_avg[[depth]] = tmp
    }
  }
  
  ctl[[i]] <- run_biocro(
    soybean$initial_values,
    soybean$parameters,
    weather_growing_season,
    soybean$direct_modules,
    soybean$differential_modules,
    soybean$ode_solver
  )
  
  direct_modules_soil_water = list("BioCroWater:soil_type_selector", "BioCroWater:soil_surface_runoff",
                                    "BioCroWater:soil_water_downflow", "BioCroWater:soil_water_tiledrain", 
                                    "BioCroWater:soil_water_upflow","BioCroWater:soil_water_uptake",
                                    "BioCroWater:multilayer_soil_profile_avg")
  
  differential_modules_soil_water = list("BioCroWater:soil_evaporation","BioCroWater:multi_layer_soil_profile")
  
  direct_modules_new = soybean$direct_modules
  old_soil_evapo_index = which(direct_modules_new=="BioCro:soil_evaporation")
  direct_modules_new = direct_modules_new[-old_soil_evapo_index] #remove soil evaporation
  direct_modules_new = c(direct_modules_new[1:(old_soil_evapo_index-1)],direct_modules_soil_water,
                         direct_modules_new[old_soil_evapo_index:length(direct_modules_new)]) # insert BioCroWater
  
  # direct_modules_new$stomata_water_stress = "BioCro:stomata_water_stress_exponential"
  # direct_modules_new$stomata_water_stress = "BioCro:stomata_water_stress_sigmoid"
  
  differential_modules_new = soybean$differential_modules
  old_soil_profile_index = which(differential_modules_new=="BioCro:two_layer_soil_profile")
  differential_modules_new = differential_modules_new[-old_soil_profile_index] #remove soil_profile
  differential_modules_new = c(differential_modules_new[1:(old_soil_profile_index-1)],
                               differential_modules_soil_water,
                               differential_modules_new[old_soil_profile_index:length(differential_modules_new)])
  
  init_values =   within(soybean$initial_values,{
    soil_water_content_1 = obs_avg[[1]]$sm[1]  #0-5cm
    soil_water_content_2 = obs_avg[[1]]$sm[1]  #5-15cm
    soil_water_content_3 = obs_avg[[2]]$sm[1]  #15-35cm
    soil_water_content_4 = obs_avg[[4]]$sm[1]  #35-55 cm
    soil_water_content_5 = obs_avg[[6]]$sm[1]  #55-75 cm
    soil_water_content_6 = obs_avg[[7]]$sm[1]  #75-100cm
    sumes1               = 0.125
    sumes2               = 0
    time_factor = 0
    soil_evaporation_rate = 0
  })
  init_values$soil_water_content=NULL
  parameters =   within(soybean$parameters, {
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
  parameters[c('soil_field_capacity','soil_saturated_conductivity','soil_saturation_capacity','soil_wilting_point')]=NULL
  parameters[c('kShell','net_assimilation_rate_shell')] = NULL
  parameters$kd = parameters$k_diffuse
  parameters[arg_names] = new_par
  
  if(use_varying_sla){
    direct_modules_new$specific_leaf_area = "BioCroWater:varying_SLA"
  }
  
  if(use_predefined){
    parameters$iSp = 3.0
  }
  
  if(use_pre_opt_vars){
    parameters[pre_arg_names] = pre_par
  }
  
  # parameters$phi2 = phi2
  direct_modules_new$leaf_water_stress = NULL
  parameters$LeafWS = 1
  
  if(use_new_water){
    results[[i]] <- run_biocro(
      init_values,
      parameters,
      weather_growing_season,
      direct_modules_new,
      differential_modules_new
    )
    # doy       = results[[i]]$doy
    # kLeaf     = results[[i]]$kLeaf
    # doy_leaf  = doy[kLeaf>0.01][1]
    # print(doy_leaf-doy[1])
  }else{
    parameters =  soybean$parameters
    parameters$iSp = 3.0
    parameters[arg_names] = new_par
    results[[i]] <- run_biocro(
      soybean$initial_values,
      parameters,
      weather_growing_season,
      soybean$direct_modules,
      soybean$differential_modules
    )
  }
}

# i=1
# plot(results[[i]]$time,results[[i]]$Grain,type='l',ylim=c(0,8))
# lines(ctl[[i]]$time,ctl[[i]]$Grain,col='red')
# 
# plot(results[[i]]$time,results[[i]]$lai,type='l')
# lines(ctl[[i]]$time,ctl[[i]]$lai,col='red')
# 
# plot(results[[i]]$time,results[[i]]$soil_water_content,type='l',ylim=c(0.2,0.6))
# lines(ctl[[i]]$time,ctl[[i]]$soil_water_content,col='red')
