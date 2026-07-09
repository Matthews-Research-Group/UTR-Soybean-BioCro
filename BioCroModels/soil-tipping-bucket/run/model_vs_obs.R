library(BioCro)
library(BioCroWater)
water_content_avg<-function(water_content_list,layer_thickness){
  
}
years = 2004:2011
obs_data = readRDS("data/soil_moisture_data_rearranged.rds")
obs_data$sm = obs_data$sm/100
#obs depth defition:
#layer 1: 5-15 cm; layer 2: 15-25; layer 3: 25-35.... layer 7: 65-75 cm
results = list()
for (i in 1:length(years)){
  year = years[i]
  weatherData <- 
    read.csv(paste0("WeatherData/",year,"/",year,"_Bondville_IL_daylength.csv"))
  weather_growing_season <- get_growing_season_climate(weatherData, threshold_temperature = 0)
  weather_growing_season$time_zone_offset = -6
  
  ctl <- run_biocro(
    soybean$initial_values,
    soybean$parameters,
    weather_growing_season,
    soybean$direct_modules,
    soybean$differential_modules,
    soybean$ode_solver
  )
  #prescribed et with BioCro's predictions
  weather_growing_season$canopy_transpiration_rate = ctl$canopy_transpiration_rate
  
  soiltype = 3
  direct_modules_new = list("BioCroWater:soil_type_selector", "BioCroWater:soil_surface_runoff",
                            "BioCroWater:soil_water_downflow", "BioCroWater:soil_water_tiledrain", 
                            "BioCroWater:soil_water_upflow","BioCroWater:soil_water_uptake",
                            "BioCroWater:multilayer_soil_profile_avg")
  
  differential_modules_new = list("BioCroWater:soil_evaporation","BioCroWater:multi_layer_soil_profile")
  
  obs_year_i = obs_data[obs_data$year==year,]
  
  obs_avg = list()
  for (depth in 1:7){
    obs = obs_year_i[obs_year_i$depth==depth,] 
    obs = obs[,c('ring','doy','sm')]
    tmp = aggregate(sm ~ doy, data = obs, FUN = mean)
    obs_avg[[depth]] = tmp
  }
  
  #further subset from the start of obs DOY
  weather_growing_season = weather_growing_season[weather_growing_season$doy>=obs_year_i$doy[1],]
  
  init_values =   within(miscanthus_x_giganteus$initial_values,{
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
    lai = 5
    bare_soil_albedo            = 0.15
    max_rooting_layer = 4
  })
  parameters[c('soil_field_capacity','soil_saturated_conductivity','soil_saturation_capacity','soil_wilting_point')]=NULL
  
  results[[i]] <- run_biocro(
    init_values,
    parameters,
    weather_growing_season,
    direct_modules_new,
    differential_modules_new
  )
}
# 
# results_new = readRDS("results_new.rds")
# y1 = results[[3]]
# y2 = results_new[[3]]
# plot( y1$time,y1$soil_water_content_1,type='l',ylim=c(0.2,0.5))
# lines(y2$time,y2$soil_water_content_1,col='red')
# legend(180,0.3,legend=c("old", "new"),
#        col=c("black","red"), lty=1, cex=1.2,bty="n",lwd = 2,
#        seg.len=0.5, x.intersp = 0.2, y.intersp = 0.3)
# stop()

#plotting
# Open a PDF device with custom size
pdf("myplot_with_uptake.pdf", width = 8, height = 6)  # Width: 8 inches, Height: 6 inches
par(mfrow=c(2,2),mar=c(4,4,3,3),mai=c(1,1,0.5,0.5),cex=1.5) #c(bottom, left, top, right)
for (i in 1:4){
  year   = years[i]
  result = results[[i]]
  obs_year_i = obs_data[obs_data$year==year,]
  weatherData <- weather[[as.character(year)]]
  linewidth=2
  axislabelspace = 2
  plot(result$time, result$soil_water_content,type='l',
       lwd = linewidth,
       xlab='',ylab='',ylim=c(0.2,0.5),
       main = paste0("year ",year))
  
  obs = obs_year_i[obs_year_i$depth<=5,] #top 50cm average for rings and depths
  
  obs_avg = aggregate(sm ~ doy, data = obs, FUN = mean)
  lines(obs_avg$doy,obs_avg$sm,col="red",lwd = linewidth,)
  #add precipitation
  # lines(weather_growing_season$doy,weather_growing_season$precip/max(weather_growing_season$precip)/2,col="blue")
  if(i==1){
  # legend(x=240,y=0.6, legend=c("BioCro", "OBS","Precip"),
  #        col=c("black","red","blue"), lty=1, cex=0.8,bty="n",lwd = linewidth,
  #        seg.len=0.5, x.intersp = 0.2, y.intersp = 0.2)
  legend(x=230,y=0.5, legend=c("BioCro", "OBS"),
           col=c("black","red"), lty=1, cex=0.8,bty="n",lwd = linewidth,
           seg.len=0.5, x.intersp = 0.2, y.intersp = 0.2)
  }
  title(xlab = "time", line = axislabelspace)            # Add x-axis text
  title(ylab = "soil water content", line = axislabelspace)            # Add y-axis text
  # lines(result$time,result$soil_evaporation_rate/2,col="blue")
}

dev.off()

stop("stopped here")

#single layer
par(mfrow=c(2,2),mar=c(4,4,3,3),mai=c(1,1,0.5,0.5),cex=1.5) #c(bottom, left, top, right)
layer_name = "soil_water_content_2"
obs_layer  = 1
for (i in 1:4){
  year   = years[i]
  result = results[[i]]
  obs_year_i = obs_data[obs_data$year==year,]
  linewidth=2
  axislabelspace = 2
  plot(result$time, result[,layer_name],type='l',
       lwd = linewidth,
       xlab='',ylab='',ylim=c(0.,0.5),
       main = paste0("year ",year))
  
  obs = obs_year_i[obs_year_i$depth==obs_layer,] 
  obs = obs[,c('ring','doy','sm')]
  
  #top x cm average for rings
  obs_avg = aggregate(sm ~ doy, data = obs, FUN = mean)
  lines(obs_avg$doy,obs_avg$sm,col="red",lwd = linewidth)
  #add precipitation
  # lines(weather_growing_season$doy,weather_growing_season$precip/max(weather_growing_season$precip)/2,col="blue")
  if(i==1){
    # legend(x=240,y=0.6, legend=c("BioCro", "OBS","Precip"),
    #        col=c("black","red","blue"), lty=1, cex=0.8,bty="n",lwd = linewidth,
    #        seg.len=0.5, x.intersp = 0.2, y.intersp = 0.2)
    legend(x=230,y=0.5, legend=c("BioCro", "OBS"),
           col=c("black","red"), lty=1, cex=0.8,bty="n",lwd = linewidth,
           seg.len=0.5, x.intersp = 0.2, y.intersp = 0.2)
  }
  title(xlab = "time", line = axislabelspace)            # Add x-axis text
  title(ylab = layer_name, line = axislabelspace)            # Add y-axis text
  # lines(result$time,result$soil_evaporation_rate/2,col="blue")
}
