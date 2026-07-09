years = 2004:2011
weather_obs  = read.csv("data/doi_10_5061_dryad_g0v62__v20170815/Final_Data_Deposit/ExtendedDataFig1_Weather_Data_2004thru2011/soyFACE_weather_data_2004thru2011.csv")
soil_moisture_obs = c()
for (year in years){
    soil_moisture_file = read.csv(paste0("data/doi_10_5061_dryad_g0v62__v20170815/Final_Data_Deposit/Fig_2_Table_S1_Soil_Moisture/",
                                         year,"_soil_moisture_ringmeans.csv"))
#use ambient CO2 only. the data are not consistent in this column. Some years it uses text "amb" instead of actual CO2 number
    soil_moisture_sub  = soil_moisture_file[soil_moisture_file$CO2<550 | soil_moisture_file$CO2=="amb",]
#total 7 layers
#1: 5-15cm; 2: 15-25;...;7: 65-75cm
    for (i in 1:7){ 
        if ("doy" %in% colnames(soil_moisture_sub)){
          tmp =  data.frame(ring = soil_moisture_sub$ring,year=year, doy=soil_moisture_sub$doy,
                          depth=i,sm = soil_moisture_sub[,i+5])
        }else{
          tmp =  data.frame(ring = soil_moisture_sub$ring,year=year, doy=soil_moisture_sub$DOY,
                            depth=i,sm = soil_moisture_sub[,i+5])  
        }
        soil_moisture_obs  = rbind(soil_moisture_obs,tmp)
    }
}
soil_moisture_obs$sm  = as.numeric(soil_moisture_obs$sm)
soil_moisture_obs$ring=as.character(soil_moisture_obs$ring)

#top 50cm average
soil_moisture_avg = aggregate(sm ~ doy,
                              data=soil_moisture_obs[soil_moisture_obs$year==2004 & soil_moisture_obs$depth<=5,],
                              FUN=mean)
# test plot
# ggplot(soil_moisture_obs[soil_moisture_obs$depth==1 & soil_moisture_obs$year==2004,], aes(x=doy, y=sm,colour=ring)) + geom_line()
ggplot(soil_moisture_avg, aes(x=doy, y=sm)) + geom_line()
# weather_obs_2004 = weather_obs[weather_obs$Year==2004,]
# plot(weather_obs_2004$DOY,weather_obs_2004$precip.mm.,type='l')

saveRDS(soil_moisture_obs,file = "data/soil_moisture_data_rearranged.rds")
