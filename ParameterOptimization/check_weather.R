
save_to_pdf <- FALSE

library(BioCro)
library(DEoptim)
library(lattice)
library(latticeExtra)
library(ggplot2)
library(rstudioapi)
library(plantecophys)

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Cost function
source('soybean_parameter_expansion.R')
source('multiyear_BioCro_optim_lsrg.R')

# set year and CO2 level
co2_opt = '_ambient_' # '_ambient_' or '_co2_'
# years, sowing dates, and harvesting dates of growing seasons being fit to
year <- c('2002', '2005', '2004', '2006')
t.sow.date <- c(152, 148, 149, 148)
sow.date <- c(179, 170, 170, 172) #27th day, 22th day, 21st day,  24th day
harv.date <- c(288, 270, 289, 270)
weather.growingseason <- list()
# names of fitted parameters
arg_names <- c('Leaf_structural_carbon_to_mass_factor', 'Leaf_utilization_rate_constant',
               'Leaf_utilization_km', 'Leaf_respiration_factor',
               'Stem_structural_carbon_to_mass_factor', 'Stem_utilization_rate_constant',
               'Stem_utilization_km', 'Stem_respiration_factor',
               'Root_structural_carbon_to_mass_factor', 'Root_utilization_rate_constant',
               'Root_utilization_km', 'Root_respiration_factor',
               'Pod_structural_carbon_to_mass_factor', 'Pod_utilization_rate_constant',
               'Pod_utilization_km', 'Pod_respiration_factor',
               'substrate_conductance_Leaf_to_Stem', 'substrate_conductance_Stem_to_Root',
               'substrate_conductance_Stem_to_Pod', 'transportation_beta_exponent',
               'Leaf_senescence_rate_max','Stem_senescence_rate_max',
               'Leaf_senescence_alpha','Leaf_senescence_beta',
               'Stem_senescence_alpha','Stem_senescence_beta', 
               'Pod_start_dvi', 'stop_growth_dvi')

# Initialize variables for the cost function
soybean_optsolver <- list()
ExpBiomass <- list()
ExpBiomass.std <- list()
RootVals <- list()
weights <- list()
numrows <- vector()

setwd('../Data/Soybean-BioCro_Parameters')
# load parameter files
param_files <- list.files(pattern = "[.]R$", recursive = TRUE)
param_files <- param_files[-c(which(param_files == "soybean_initial_values.R"))]
lapply(param_files, source)

for (i in 1:length(year)) {   # Remember to change back 
  yr <- year[i]
  weather <- read.csv(file = paste0('../Weather_data/', yr,'_Bondville_IL_daylength_wDVI.csv'))
  sd.ind <- which(weather$doy == sow.date[i])[1]
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather.growingseason[[i]] <- weather[sd.ind:hd.ind,]
  # if(i == 1){
  #   weather.growingseason[[1]]$DVI = weather.growingseason[[1]]$DVI - 0.07
  # }else if(i == 2){
  #   weather.growingseason[[2]]$DVI = weather.growingseason[[2]]$DVI + 0.07
  # }
  
  ExpBiomass[[i]] <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass.csv'))
  colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  ExpBiomass.std[[i]] <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[3], "Root"=0.17*sum(ExpBiomass[[i]][5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  numrows[i] <- nrow(weather.growingseason[[i]])
  invwts <- ExpBiomass.std[[i]]
  weights[[i]] <- log(1/(invwts[,2:ncol(invwts)]+1e-5))
  
  # load parameter files
  source('soybean_initial_values.R')
  soybean_optsolver[[i]] <- partial_run_biocro(initial_state,
                                               parameters,
                                               weather.growingseason[[i]],
                                               steady_state_module_names,
                                               derivative_module_names,
                                               solver,
                                               arg_names,
                                               verbose = FALSE)
  VPD <- RHtoVPD(weather.growingseason[[i]]$h, weather.growingseason[[i]]$temp, Pa = 101)
  xyplot(VPD~weather.growingseason[[i]]$doy)
}

VPD <- RHtoVPD(weather.growingseason[[i]]$rh, weather.growingseason[[i]]$temp, Pa = 101)
plot(VPD~weather.growingseason[[i]]$doy)

i = 1
weather.growingseason[[i]]$too_cold[weather.growingseason[[i]]$temp<13]='Blue'
weather.growingseason[[i]]$too_cold[which(weather.growingseason[[i]]$temp>=13)]='Black'
weather.growingseason[[i]]$too_cold[weather.growingseason[[i]]$temp>=18]='darkgreen'
weather.growingseason[[i]]$too_cold[which(weather.growingseason[[i]]$temp>=30)]='Red'
plot(data = weather.growingseason[[i]], temp~doy, col = too_cold, xlim = c(200,280), ylim=c(5,35))

not_really_cold <- unique(weather.growingseason[[i]]$doy[weather.growingseason[[i]]$doy == 
                                     unique(weather.growingseason[[i]]$doy
                                            [weather.growingseason[[i]]$temp<13]) 
                                   & weather.growingseason[[i]]$temp > 18])
nights_cold <- unique(weather.growingseason[[i]]$doy[weather.growingseason[[i]]$too_cold=='Blue'])
print(night_cold[!(night_cold %in% not_really_cold)])
