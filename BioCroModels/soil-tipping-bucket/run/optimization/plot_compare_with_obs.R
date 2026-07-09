library(BioCro)
library(ggplot2) # To ggplot functions
library(reshape2) # For melt function
library(data.table) # For first and last functions
library(gridExtra)
rm(list=ls())
source('set_modules_for_soil_water.R')
source('plot_funcs.R')
use_new_water   = FALSE 
#only need to change this name. 
plot_name = 'old_water_v102'

years <- c('2002','2004','2005','2006')

# sowing and harvest DOYs for each growing season
dates <- data.frame("year" = c(2002, 2004:2006),"sow" = c(152,149,148,148), "harvest" = c(288, 289, 270, 270))

if(use_new_water){
  soybean_steadystate_modules0 = set_direct_modules() 
  soybean_derivative_modules0  = set_differential_modules() 
  soybean_initial_state0       = set_init_values() 
  soybean_parameters0          = set_parameters() 
  soybean_parameters0$kd = soybean_parameters0$k_diffuse
  soybean_solver_params  = soybean$ode_solver
}else{
  soybean_steadystate_modules0 = soybean$direct_modules 
  soybean_derivative_modules0  = soybean$differential_modules 
  soybean_initial_state0       = soybean$initial_values 
  soybean_parameters0          = soybean$parameters 
  soybean_solver_params        = soybean$ode_solver
}

# initialize variables
results <- list()
results_CTL <- list()
weather_growing_season <- list()
ExpBiomass <- list()
ExpBiomass.std <- list()
soybean_solver <- list()
RootVals <- list()
weights <- list()
numrows <- vector()

for (i in 1:length(years)) {
  yr <- years[i]
  weather <- read.csv(file = paste0('Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))

  sowdate <- dates$sow[which(dates$year == yr)]
  harvestdate <- dates$harvest[which(dates$year == yr)]
  sd.ind <- which(weather$doy == sowdate)[1]
  hd.ind <- which(weather$doy == harvestdate)[24]
  
  weather_growing_season <- weather[sd.ind:hd.ind,]
  weather_growing_season$time_zone_offset = -6

#calculate CTL for diagnostic purposes
  results[[i]] <- run_biocro(
    soybean$initial_values,
    soybean$parameters,
    weather_growing_season,
    soybean$direct_modules,
    soybean$differential_modules,
    soybean$ode_solver
  )
  results[[i]] <- run_biocro(soybean_initial_state0, 
 			     soybean_parameters0, 
			     weather_growing_season,
                             soybean_steadystate_modules0,
			     soybean_derivative_modules0, 
                             soybean_solver_params)
  results[[i]]$time = results[[i]]$doy + results[[i]]$hour/24
  
  ExpBiomass[[i]] <- read.csv(file=paste0('Data/biomasses_with_seed/',yr,'_ambient_biomass.csv'))
  colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem","Shell0","Seed","Litter","CumLitter")
  Shell = ExpBiomass[[i]]$Shell0 - ExpBiomass[[i]]$Seed
#  Shell[which.max(Shell):length(Shell)] = max(Shell) #make Shell not decline
  ExpBiomass[[i]]$Shell = Shell 

  
  ExpBiomass.std[[i]] <- read.csv(file=paste0('Data/biomasses_with_seed/',yr,'_ambient_biomass_std.csv'))
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Shell0","Seed","Litter","CumLitter")
  ExpBiomass.std[[i]]$Shell = sqrt((ExpBiomass.std[[i]]$Shell0)^2 + (ExpBiomass.std[[i]]$Seed)^2)
  
  RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[5], "Root"=0.17*sum(ExpBiomass[[i]][5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  numrows[i] <- nrow(weather_growing_season)
  invwts <- ExpBiomass.std[[i]]
  weights[[i]] <- log(1/(invwts[,2:ncol(invwts)]+1e-5))

}

plot_list = list()
plot_list_litter = list()
for (i in 1:length(years)){
	yr = years[i]
	FigA <- plot_all_tissues(results[[i]], yr, ExpBiomass[[i]], ExpBiomass.std[[i]])
	plot_list[[i]] = FigA

        
        results[[i]]$TotalLitter = results[[i]]$LeafLitter + results[[i]]$StemLitter

	FigB <- plot_litters(results[[i]],yr, ExpBiomass[[i]], ExpBiomass.std[[i]])
	plot_list_litter[[i]] = FigB
}

pdf_name = paste0("figs/fig_",plot_name,".pdf")
pdf(pdf_name,height = 8, width=8,bg='transparent')
grid.arrange(grobs = plot_list,nrow=2,ncol=2)
dev.off()

#litter biomass plot
pdf_name = paste0("figs/fig_",plot_name,"_litter.pdf")
pdf(pdf_name,height = 8, width=8,bg='transparent')
grid.arrange(grobs = plot_list_litter,nrow=2,ncol=2)
dev.off()
