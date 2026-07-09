library(BioCro)
library(BioCroWater)
rm(list=ls())
source('functions_BioCro_run.R')
para = c(soybean$parameters$b0,soybean$parameters$b1,110,195)
use_varying_sla  = FALSE
use_pre_opt_vars = TRUE  #use root partition coefs
results_all = my_biocro_run(para,use_varying_sla,use_pre_opt_vars)
ctls     = results_all[[2]]
results  = results_all[[1]]
obs_data = readRDS("data/soil_moisture_data_rearranged.rds")
obs_data$sm = obs_data$sm/100

depth2compare = c(1,4,7)  #10 cm, 40 cm, 70 cm
model_num_layers = c(2,3,5)

y=results[[1]]
# plot(y$doy+y$hour/24,y$soil_evaporation_rate)
# plot(y$doy+y$hour/24,y$surface_runoff)

for (i in 1:length(depth2compare)){
  dep = depth2compare[i]
  #top x-cm average for rings and depths
  obs_data_sub = obs_data[obs_data$depth<=dep,]
  obs_avg = aggregate(sm ~ doy+year, data = obs_data_sub, FUN = mean)
  obs_sd  = aggregate(sm ~ doy+year, data = obs_data_sub, FUN = sd)
  obs     = cbind(obs_avg,obs_sd$sm)
  #plotting
  years = 2004:2011
  
  source('myfunctions.R')
  out_fn = paste0("figs/swc_TimeSeries_v127_",dep*10,"cm.png")
  var2plot = 'soil_water_content'
  plot_general_time_series(out_fn,var2plot,years,obs,results,model_num_layers[i])
  
  out_fn = paste0("figs/swc_scatterplot_v127_",dep*10,"cm.png")
  var2plot = 'soil_water_content'
  plot_general_scatter(out_fn,years,obs,results,var2plot,model_num_layers[i])
}

var2plot = "StomataWS"
results_i = results[[1]]
mytime = results_i$doy + results_i$hour/24
plot(mytime,results_i[,var2plot],type='l',col='green',
     main=paste0(var2plot),ylab=paste0(var2plot),ylim=c(0,1))
results_i = results[[2]]
mytime = results_i$doy + results_i$hour/24
lines(mytime,results_i[,var2plot],col='red')
results_i = results[[3]]
mytime = results_i$doy + results_i$hour/24
lines(mytime,results_i[,var2plot],col='blue')
# lines(mytime,rep(0.33,length(mytime)))
legend(x=160,y=0.5, legend=c("2004","2005","2006"),
       col=c("green","red","blue"), lty=1, cex=0.8,bty="n",lwd = 2,
       seg.len=1.0, x.intersp = 0.5, y.intersp = 0.5)


if(FALSE){
  your_data = results[[1]]
  daily_evap <- your_data %>%
    group_by(doy) %>%
    summarise(evaporation_daily = sum(soil_evaporation_rate, na.rm = TRUE))
# Open a PDF device with custom size
for (j in 1:2){
  pdf(paste0("figs/soybean_water_v119_p",j,".pdf"), width = 8, height = 6)  # Width: 8 inches, Height: 6 inches
  par(mfrow=c(2,2),mar=c(4,4,3,3),mai=c(1,1,0.5,0.5),cex=1.5) #c(bottom, left, top, right)
  i0 = (j-1)*4+1
  i1 = j*4
  for (i in i0:i1){
    year   = years[i]
    result = results[[i]]
    ctl    = ctls[[i]]
    biocro_time = result$doy + result$hour/24
    obs_year_i = obs_data[obs_data$year==year,]
    weatherData <- weather[[as.character(year)]]
    linewidth=2
    axislabelspace = 2
    plot(biocro_time, result$soil_water_content,type='l',col='blue',
         lwd = linewidth,
         xlab='',ylab='',ylim=c(0.2,0.5),
         main = paste0("year ",year))

    # lines(biocro_time, ctl$soil_water_content,col = 'red',lwd = linewidth)
    
    obs = obs_year_i[obs_year_i$depth<=5,] #top 50cm average for rings and depths

    obs_avg = aggregate(sm ~ doy, data = obs, FUN = mean)
    lines(obs_avg$doy,obs_avg$sm,col="black",lwd = linewidth)
    #add precipitation
    # lines(weather_growing_season$doy,weather_growing_season$precip/max(weather_growing_season$precip)/2,col="blue")
    if(i==1){
      # legend(x=240,y=0.6, legend=c("BioCro", "OBS","Precip"),
      #        col=c("black","red","blue"), lty=1, cex=0.8,bty="n",lwd = linewidth,
      #        seg.len=0.5, x.intersp = 0.2, y.intersp = 0.2)
      legend(x=230,y=0.5, legend=c("BioCro", "OBS"),
             col=c("blue","black"), lty=1, cex=0.8,bty="n",lwd = linewidth,
             seg.len=0.5, x.intersp = 0.2, y.intersp = 0.7)
    }
    title(xlab = "time", line = axislabelspace)            # Add x-axis text
    title(ylab = "soil water content", line = axislabelspace)            # Add y-axis text
    # lines(result$time,result$soil_evaporation_rate/2,col="blue")
  }

  dev.off()
}
}
