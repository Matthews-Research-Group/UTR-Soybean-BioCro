library(BioCro)
# library(BioCroWater)
library(dplyr)
library(tidyr)
# Define the RMSE function  
rmse <- function(observed, predicted) {
  # Check if the lengths of the observed and predicted vectors are the same
  if (length(observed) != length(predicted)) {
    stop("The lengths of observed and predicted values must be the same.")
  }                         
  
  # Calculate the squared differences
  squared_diff <- (observed - predicted)^2
  
  # Calculate the mean of the squared differences
  mean_squared_diff <- mean(squared_diff)
  
  # Calculate the square root of the mean squared differences (RMSE)
  rmse_value <- sqrt(mean_squared_diff)
  
  return(rmse_value)
} 

calc_vpd <- function(temp, rh) {
  es <- 0.6108 * exp((17.27 * temp) / (temp + 237.3))  # kPa
  ea <- es * rh  # kPa
  vpd <- es - ea  # kPa
  return(vpd)
}
# sla<-function(dvi){
#   if(dvi<0.55){
#     sla = 2.58
#   }else if(dvi>1.55){
#     sla = 2
#   }else{
#    sla = -2.392922+12.452842*dvi-6.202293*dvi^2
#   }
#   return(sla)
# }

linear_regression<-function(observed, predicted) {
  # 2. Fit linear model:
  fit <- lm(predicted ~ observed)
  
  slope <- coef(fit)[2]
  intercept <- coef(fit)[1]
  r2 <- summary(fit)$r.squared
  
  return(c(slope,intercept,r2))
}

get_CO2<-function(year){
  #co2 ppm data of SSP5-8.5 from 1970-2060. Yufeng
  years = 1970:2060
  co2_level =
    c(325.58,326.23,328.21,330.81,331.67,331.74,332.50,334.34,336.02,337.66,
      339.75,341.08,341.65,343.24,344.83,346.31,347.84,349.58,352.26,354.00,
      355.14,356.49,357.22,357.89,359.44,361.50,363.41,364.58,367.28,369.33,
      370.70,372.26,374.41,377.16,378.64,380.61,382.83,384.31,386.50,388.02,
      390.78,393.04,395.04,397.71,399.59,401.68,405.11,407.85,410.85,413.92,
      417.06,420.29,423.62,427.07,430.63,434.30,438.08,441.97,445.97,450.08,
      454.29,458.63,463.10,467.70,472.43,477.29,482.28,487.41,492.67,498.06,
      503.60,509.27,515.09,521.07,527.19,533.47,539.91,546.50,553.25,560.16,
      567.24,574.48,581.92,589.55,597.38,605.40,613.63,622.05,630.68,639.51,648.54)
  return(co2_level[years==year])
}

#Compute total daily precipitation.
#Set precipitation to 0 for all hours except 3 randomly chosen ones, 
#to which you allocate equal parts of the daily total.
redistribute_precip<-function(weather_data){
  modified_precip <- weather_data %>%
    group_by(doy) %>%
    mutate(
      daily_total = sum(precip),
      has_precip = any(precip > 0)
    ) %>%
    group_by(doy) %>%
    mutate(
      precip = if (has_precip[1]) {
        # Choose 3 random hours from the day
        hours_to_keep <- sample(hour, 3)
        ifelse(hour %in% hours_to_keep, daily_total[1] / 3, 0)
      } else {
        0
      }
    ) %>%
    select(-daily_total, -has_precip) %>%
    ungroup()
}

my_biocro_run<-function(x,use_varying_sla,use_pre_opt_vars,withDrought=FALSE){
  set.seed(123)  # for reproducibility
  use_old_soil_evapo = FALSE
  using_Gray_precip  = TRUE
  obs_weather_Gray<-
    read.csv('data/doi_10_5061_dryad_g0v62__v20170815/Final_Data_Deposit/ExtendedDataFig1_Weather_Data_2004thru2011/soyFACE_weather_data_2004thru2011.csv')
  opt_result_name = "new_water_v127"
  new_par = readRDS(paste0('optimization/opt_results/opt_result_',opt_result_name,'.rds'))
  arg_names = names(new_par)
  print(new_par)
  if(use_pre_opt_vars){
    pre_arg_names <- c('alphaRoot','betaRoot')
    
    pre_par = readRDS(paste0('optimization/opt_results/opt_result_new_water_v108.rds'))
    pre_par = pre_par$optim$bestmem
    pre_par = pre_par[c(2,5)] #'alphaRoot','betaRoot'
  }
  
  years = 2004:2011
  obs_data = readRDS("data/soil_moisture_data_rearranged.rds")
  obs_data$sm = obs_data$sm/100
  #obs depth defition:
  #layer 1: 5-15 cm; layer 2: 15-25; layer 3: 25-35.... layer 7: 65-75 cm
  results = list()
  ctl = list()
  for (i in 1:length(years)){
    year = years[i]
    obs_weather_Gray_yeari = obs_weather_Gray[obs_weather_Gray$Year==year,]
    weatherData <- 
      read.csv(paste0("WeatherData/",year,"/",year,"_Bondville_IL_daylength.csv"))
    #further subset from the start of obs DOY
    DOY_start = obs_weather_Gray_yeari$DOY[1]
    DOY_end = tail(obs_weather_Gray_yeari$DOY,1)
    print(c(year,DOY_start))
    weather_growing_season = weatherData[weatherData$doy>=DOY_start & weatherData$doy<=DOY_end,]
    weather_growing_season$time_zone_offset = -6
    
    if(using_Gray_precip){
      obs_weather_Gray_yeari_hourly <- obs_weather_Gray_yeari %>%
        # Add an hourly sequence per day
        uncount(weights = 24, .id = "hour") %>%
        # Adjust the hour (0 to 23)
        mutate(hour = hour - 1,
               precip = precip.mm. / 24)
      
      # replace precip with Gray's
      # if(year==2006){
      weather_growing_season$precip = obs_weather_Gray_yeari_hourly$precip
      # }
    }
    # 
    # weather_growing_season <- redistribute_precip(weather_growing_season)
    # mimic Gray et al's awning rainfall interception
    if(withDrought){
      print("Warning: this is a case to intercept rainfall")
      weather_growing_season$precip[weather_growing_season$solar<50 & weather_growing_season$windspeed<10] = 0
    }
    
    #get obs data
    obs_year_i = obs_data[obs_data$year==year,]
    obs_avg = list()
    for (depth in 1:7){
      obs = obs_year_i[obs_year_i$depth==depth,] 
      obs = obs[,c('ring','doy','sm')]
      tmp = aggregate(sm ~ doy, data = obs, FUN = mean)
      obs_avg[[depth]] = tmp
    }
    
    ctl[[i]] <- run_biocro(
      soybean$initial_values,
      soybean$parameters,
      weather_growing_season,
      soybean$direct_modules,
      soybean$differential_modules,
      soybean$ode_solver
    )
    

    direct_modules_soil_water = list("BioCroWater:soil_type_selector",
                                     "BioCroWater:soil_surface_runoff",  "BioCroWater:soil_water_downflow",
                                     "BioCroWater:soil_water_tiledrain", "BioCroWater:soil_water_upflow",
                                     "BioCroWater:soil_water_uptake",    "BioCroWater:multilayer_soil_profile_avg")
    if(use_old_soil_evapo){
      differential_modules_soil_water = list(#"BioCroWater:surface_excess_water",
                                           "BioCroWater:multi_layer_soil_profile")
    }else{
      differential_modules_soil_water = list(#"BioCroWater:surface_excess_water",
                                             "BioCroWater:soil_evaporation",
                                             "BioCroWater:multi_layer_soil_profile")
    }
    
    direct_modules_new = soybean$direct_modules
    old_soil_evapo_index = which(direct_modules_new=="BioCro:soil_evaporation")
    if(use_old_soil_evapo){
      direct_modules_new = c(direct_modules_new[1:(old_soil_evapo_index)],direct_modules_soil_water,
                             direct_modules_new[(old_soil_evapo_index+1):length(direct_modules_new)]) # insert BioCroWater
    }else{
      direct_modules_new = direct_modules_new[-old_soil_evapo_index] #remove soil evaporation
      direct_modules_new = c(direct_modules_new[1:(old_soil_evapo_index-1)],direct_modules_soil_water,
                           direct_modules_new[old_soil_evapo_index:length(direct_modules_new)]) # insert BioCroWater
    }
    
    differential_modules_new = soybean$differential_modules
    old_soil_profile_index = which(differential_modules_new=="BioCro:two_layer_soil_profile")
    differential_modules_new = differential_modules_new[-old_soil_profile_index] #remove soil_profile
    differential_modules_new = c(differential_modules_new[1:(old_soil_profile_index-1)],
                                 differential_modules_soil_water,
                                 differential_modules_new[old_soil_profile_index:length(differential_modules_new)])
    
    if(use_old_soil_evapo){
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
        # soil_evaporation_rate = 0
      })
    }else{
      init_swc = 0.25
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
    }
    init_values$soil_water_content=NULL
    soiltype = 4
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
      tile_drainage_rate          = 0.2 # dimensionless, as in DSSAT (1/d ?). typical 5-6 mm/day
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
    
    parameters$b0 = x[1]
    parameters$b1 = x[2]
    parameters$Catm = get_CO2(year)
    
    parameters$Vcmax_at_25 = x[3]
    parameters$Jmax_at_25  = x[4]
    
    if(length(x) == 6){
      parameters$RL_at_25        = x[5]
      parameters$Catm            = x[6]
      
      # if(year==2007){
      #   parameters$alphaLeaf = parameters$alphaLeaf * 1.1
      # }
    }
    
    if(use_varying_sla){
      parameters$iSp = parameters$iSp * 0.85
      print(c("SLA is",parameters$iSp))
      #direct_modules_new$specific_leaf_area = "BioCroWater:varying_SLA"
    }
    
    if(use_pre_opt_vars){
      parameters[pre_arg_names] = pre_par
    }
    

    
    # parameters$phi2 = 0.5
    direct_modules_new$leaf_water_stress = NULL
    parameters$LeafWS = 1
    
    # parameters$iSp = 3.0
    
    results[[i]] <- run_biocro(
      init_values,
      parameters,
      weather_growing_season,
      direct_modules_new,
      differential_modules_new
    )
  }
  return(list(results,ctl))
}  

my_biocro_run_develop<-function(x,use_varying_sla,use_pre_opt_vars,years,withDrought=FALSE){
  set.seed(123)  # for reproducibility
  using_Gray_precip  = TRUE
  opt_result_name = "new_water_v132"
  new_par = readRDS(paste0('optimization/opt_results/opt_result_',opt_result_name,'.rds'))
  arg_names = names(new_par)
  
  obs_weather_Gray<-
    read.csv('data/doi_10_5061_dryad_g0v62__v20170815/Final_Data_Deposit/ExtendedDataFig1_Weather_Data_2004thru2011/soyFACE_weather_data_2004thru2011.csv')
  
  # years = 2004:2011
  obs_data = readRDS("data/soil_moisture_data_rearranged.rds")
  obs_data$sm = obs_data$sm/100
  #obs depth defition:
  #layer 1: 5-15 cm; layer 2: 15-25; layer 3: 25-35.... layer 7: 65-75 cm
  results = list()
  ctl = list()
  for (i in 1:length(years)){
    year = years[i]
    obs_weather_Gray_yeari = obs_weather_Gray[obs_weather_Gray$Year==year,]
    weatherData <- 
      read.csv(paste0("WeatherData/",year,"/",year,"_Bondville_IL_daylength.csv"))
    #further subset from the start of obs DOY
    DOY_start = obs_weather_Gray_yeari$DOY[1]
    DOY_end = tail(obs_weather_Gray_yeari$DOY,1)
    # print(c(year,DOY_start))
    weather_growing_season = weatherData[weatherData$doy>=DOY_start & weatherData$doy<=DOY_end,]
    weather_growing_season$time_zone_offset = -6
    
    if(using_Gray_precip){
      obs_weather_Gray_yeari_hourly <- obs_weather_Gray_yeari %>%
        # Add an hourly sequence per day
        uncount(weights = 24, .id = "hour") %>%
        # Adjust the hour (0 to 23)
        mutate(hour = hour - 1,
               precip = precip.mm. / 24)
      
      # replace precip with Gray's
      # if(year==2006){
      weather_growing_season$precip = obs_weather_Gray_yeari_hourly$precip
      # }
    }
    # 
    # weather_growing_season <- redistribute_precip(weather_growing_season)
    # mimic Gray et al's awning rainfall interception
    if(withDrought){
      print("Warning: this is a case to intercept rainfall")
      weather_growing_season$precip[weather_growing_season$solar<50 & weather_growing_season$windspeed<10] = 0
    }
    
    # doy_trigger <- 1
    # weather_growing_season$precip_synth <- ifelse(weather_growing_season$doy <= doy_trigger,
    #                                         weather_growing_season$precip, 
    #                                         0)
    # weather_growing_season$precip    =  weather_growing_season$precip_synth
    
    #get obs data
    obs_year_i = obs_data[obs_data$year==year,]
    obs_avg = list()
    for (depth in 1:7){
      obs = obs_year_i[obs_year_i$depth==depth,] 
      obs = obs[,c('ring','doy','sm')]
      tmp = aggregate(sm ~ doy, data = obs, FUN = mean)
      obs_avg[[depth]] = tmp
    }
    
    ctl[[i]] <- run_biocro(
      soybean$initial_values,
      soybean$parameters,
      weather_growing_season,
      soybean$direct_modules,
      soybean$differential_modules,
      soybean$ode_solver
    )
    init_values =   within(soybean2$initial_values,{
      soil_water_content_1 = obs_avg[[1]]$sm[1]  #0-5cm
      soil_water_content_2 = obs_avg[[1]]$sm[1]  #5-15cm
      soil_water_content_3 = obs_avg[[2]]$sm[1]  #15-35cm
      soil_water_content_4 = obs_avg[[4]]$sm[1]  #35-55 cm
      soil_water_content_5 = obs_avg[[6]]$sm[1]  #55-75 cm
      soil_water_content_6 = obs_avg[[7]]$sm[1]  #75-100cm
    })
    
    parameters = soybean2$parameters
    
    parameters$b0 = x[1]
    parameters$b1 = x[2]
    parameters$Catm = get_CO2(year)
    
    parameters$Vcmax_at_25 = x[3]
    parameters$Jmax_at_25  = x[4]
    parameters$resistance_base = x[5]
    
    if(length(x) == 8){
      parameters$RL_at_25        = x[6]
      parameters$Catm            = x[7]
      parameters$resistance_amplifier = x[8]
      # if(year==2007){
      #   parameters$alphaLeaf = parameters$alphaLeaf * 1.01
      # }
    }
    
    if(use_varying_sla){
      parameters$iSp = parameters$iSp * 0.85
      # print(c("eCO2 SLA is",parameters$iSp))
      #direct_modules_new$specific_leaf_area = "BioCroWater:varying_SLA"
    }
    
    direct_modules = soybean2$direct_modules
    # direct_modules$stomata_water_stress = NULL
    # parameters$StomataWS = 1
    # 
    direct_modules$stomata_water_stress = "BioCro:stomata_water_stress_resistance"
    # init_values$uptake_laststep = 0
    
    parameters[arg_names] = new_par
    
    results[[i]] <- run_biocro(
      init_values,
      parameters,
      weather_growing_season,
      direct_modules,
      soybean2$differential_modules
    )
  }
  return(list(results,ctl))
}  

my_biocro_run_simple<-function(co2_input){
  years <- c('2002','2004','2005','2006')
  
  # sowing and harvest DOYs for each growing season
  dates <- data.frame("year" = c(2002, 2004:2006),"sow" = c(152,149,148,148), "harvest" = c(288, 289, 270, 270))
  years_obs = 2004:2011
  obs_weather_Gray<-
    read.csv('data/doi_10_5061_dryad_g0v62__v20170815/Final_Data_Deposit/ExtendedDataFig1_Weather_Data_2004thru2011/soyFACE_weather_data_2004thru2011.csv')
  obs_data = readRDS("data/soil_moisture_data_rearranged.rds")
  obs_data$sm = obs_data$sm/100
  use_pre_opt_vars  = TRUE
  using_Gray_precip = TRUE
  opt_result_name = "new_water_v127"
  new_par = readRDS(paste0('optimization/opt_results/opt_result_',opt_result_name,'.rds'))
  # new_par = new_par$optim$bestmem
  # names(new_par) = arg_names
  arg_names = names(new_par)
  print(new_par)
  if(use_pre_opt_vars){
    pre_arg_names <- c('alphaRoot','betaRoot')
    
    pre_par = readRDS(paste0('optimization/opt_results/opt_result_new_water_v108.rds'))
    pre_par = pre_par$optim$bestmem
    pre_par = pre_par[c(2,5)]
  }
  
  years = c(2002,2004:2006)
  results = list()
  for (i in 1:length(years)){
    weather_growing_season = soybean_weather[[i]]
    sowdate <- dates$sow[which(dates$year == years[i])]
    harvestdate <- dates$harvest[which(dates$year == years[i])]
    if(using_Gray_precip  & years[i]>=2004){ #Gray's does not have year 2002
      obs_weather_Gray_yeari = obs_weather_Gray[obs_weather_Gray$Year==years[i],]
      obs_weather_Gray_yeari = obs_weather_Gray_yeari[obs_weather_Gray_yeari$DOY>=sowdate &obs_weather_Gray_yeari$DOY<=harvestdate,]
      obs_weather_Gray_yeari_hourly <- obs_weather_Gray_yeari %>%
        # Add an hourly sequence per day
        uncount(weights = 24, .id = "hour") %>%
        # Adjust the hour (0 to 23)
        mutate(hour = hour - 1,
               precip = precip.mm. / 24)
      
      # replace precip with Gray's
      # If Gray's data is shorter than the growing season, we only replace the shorter part
      weather_growing_season$precip[weather_growing_season$doy>=obs_weather_Gray_yeari$DOY[1] &weather_growing_season$doy<=tail(obs_weather_Gray_yeari$DOY,1)] = obs_weather_Gray_yeari_hourly$precip
    }
    
    soiltype = 4
    direct_modules_soil_water = list("BioCroWater:soil_type_selector",
                                     "BioCroWater:soil_surface_runoff",  "BioCroWater:soil_water_downflow",
                                     "BioCroWater:soil_water_tiledrain", "BioCroWater:soil_water_upflow",
                                     "BioCroWater:soil_water_uptake",    "BioCroWater:multilayer_soil_profile_avg")
    differential_modules_soil_water = list(#"BioCroWater:surface_excess_water",
        "BioCroWater:soil_evaporation",
        "BioCroWater:multi_layer_soil_profile")
    
    direct_modules_new = soybean$direct_modules
    old_soil_evapo_index = which(direct_modules_new=="BioCro:soil_evaporation")
    direct_modules_new = direct_modules_new[-old_soil_evapo_index] #remove soil evaporation
    direct_modules_new = c(direct_modules_new[1:(old_soil_evapo_index-1)],direct_modules_soil_water,
                             direct_modules_new[old_soil_evapo_index:length(direct_modules_new)]) # insert BioCroWater
    
    differential_modules_new = soybean$differential_modules
    old_soil_profile_index = which(differential_modules_new=="BioCro:two_layer_soil_profile")
    differential_modules_new = differential_modules_new[-old_soil_profile_index] #remove soil_profile
    differential_modules_new = c(differential_modules_new[1:(old_soil_profile_index-1)],
                                 differential_modules_soil_water,
                                 differential_modules_new[old_soil_profile_index:length(differential_modules_new)])
    #get obs data
    if(years[i]!=2002){
      obs_year_i = obs_data[obs_data$year==years[i],]
    }else{
      obs_year_i = obs_data[obs_data$year==2004,]
    }
    obs_avg = list()
    for (depth in 1:7){
      obs = obs_year_i[obs_year_i$depth==depth,] 
      obs = obs[,c('ring','doy','sm')]
      tmp = aggregate(sm ~ doy, data = obs, FUN = mean)
      obs_avg[[depth]] = tmp
    }
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
      tile_drainage_rate          = 0.2 # dimensionless, as in DSSAT (1/d ?). typical 5-6 mm/day
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
    
    parameters$Catm = co2_input
    
    direct_modules_new$leaf_water_stress = NULL
    parameters$LeafWS = 1
    
    if(use_pre_opt_vars){
      parameters[pre_arg_names] = pre_par
    }
    results[[i]] <- run_biocro(
      init_values,
      parameters,
      weather_growing_season,
      direct_modules_new,
      differential_modules_new
    )
  }
  return(results)
}  