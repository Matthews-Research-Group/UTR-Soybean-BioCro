

### First part of this code was taken from SoybeanBioCro-Figures.R ####
# set working directory to Soybean-BioCro/ReproduceFigures

filepath <- # path to Soybean-BioCro/ReproduceFigures
setwd(filepath)

years <- c('2002','2004','2005','2006')

# sowing and harvest DOYs for each growing season
dates <- data.frame("year" = c(2002, 2004:2006),"sow" = c(152,149,148,148), "harvest" = c(288, 289, 270, 270))

# initialize variables
results <- list()
results.elevCO2 <- list()
weather.growingseason <- list()
ExpBiomass <- list()
ExpBiomass.std <- list()
ExpBiomass.elevCO2 <- list()
ExpBiomass.std.elevCO2 <- list()
LAI <- list()
LAI.elevCO2 <- list()

# for partial_run_biocro
arg_names <- c('Catm') # atmospheric CO2 parameter
params.ambient <- c(372) # ambient atmospheric CO2
params.elevCO2 <- c(550) # elevated atmospheric CO2

for (i in 1:length(years)) {
  
  yr <- years[i]
  
  weather <- read.csv(file = paste0('../Data/Weather_data/', yr, '_Bondville_IL_daylength.csv'))
  
  sowdate <- dates$sow[which(dates$year == yr)]
  harvestdate <- dates$harvest[which(dates$year == yr)]
  sd.ind <- which(weather$doy == sowdate)[1] # start of sowing day
  hd.ind <- which(weather$doy == harvestdate)[24] # end of harvest day
  
  weather.growingseason[[i]] <- weather[sd.ind:hd.ind,]
  
  soybean_solver <- partial_run_biocro(soybean_initial_values, soybean_parameters, weather.growingseason[[i]],
                                       soybean_direct_modules, soybean_differential_modules,
                                       soybean_ode_solver, arg_names)
  
  results[[i]] <- soybean_solver(params.ambient)
  results.elevCO2[[i]] <- soybean_solver(params.elevCO2)
  
  # Load SoyFACE ambient and elevated CO2 biomass means and standard deviations
  ExpBiomass[[i]]<-read.csv(file=paste0('../Data/SoyFACE_data/', yr,'_ambient_biomass.csv'))
  ExpBiomass.std[[i]]<-read.csv(file=paste0('../Data/SoyFACE_data/', yr,'_ambient_biomass_std.csv'))
  colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem","Pod")
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Pod")
  
  ExpBiomass.elevCO2[[i]]<-read.csv(file = paste0('../Data/SoyFACE_data/',yr,'_co2_biomass.csv'))
  ExpBiomass.std.elevCO2[[i]]<-read.csv(file = paste0('../Data/SoyFACE_data/',yr,'_co2_biomass_std.csv'))
  colnames(ExpBiomass.elevCO2[[i]])<-c("DOY","Leaf","Stem","Pod")
  colnames(ExpBiomass.std.elevCO2[[i]])<-c("DOY","Leaf","Stem","Pod")
  
  # Load LAI data
  if(i>1){
    LAI[[i]]<-read.csv(file=paste0('../Data/SoyFACE_data/',yr,'_ambient_lai.csv'))
    LAI.elevCO2[[i]]<-read.csv(file=paste0('../Data/SoyFACE_data/',yr,'_elevated_lai.csv'))
  }
  
}

##### 



## calculate the error

i <- 1

r <- results[[i]]

doys <- ExpBiomass[[i]]$DOY
leaf <- ExpBiomass[[i]]$Leaf
stem <- ExpBiomass[[i]]$Stem
pod <- ExpBiomass[[i]]$Pod

inds <- NULL
for (j in 1:length(doys)){
  
  inds[j] <- which(r$time == doys[j]+0.5)
  
}


err.leaf <- (r$Leaf[inds]-leaf)
err.stem <- (r$Stem[inds]-stem)
err.pod <- (r$Grain[inds]-pod)

mse.leaf <- mean(err.leaf^2)
mse.stem <- mean(err.stem^2)
mse.pod <- mean(err.pod^2)

rmse.leaf <- sqrt(mse.leaf)
rmse.stem <- sqrt(mse.stem)
rmse.pod <- sqrt(mse.pod)

mse.biomass <- mean(c(err.leaf^2, err.stem^2, err.pod^2))
rmse.biomass <- sqrt(mse.biomass)
