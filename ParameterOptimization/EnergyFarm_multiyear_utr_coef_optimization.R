# The work directory is hard coded since the previous one doesn't work well all the time. 
# Please change it according to your own working directory.
# Parameter loading and plotting are put into this R code to check the variables more easily.
# To be moved to different files once everything works out.
# Clear workspace
rm(list=ls())

library(BioCro)
library(DEoptim)

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source files
source('soybean_parameter_expansion.R')
source('EnergyFarm_multiyear_utr_optim.R')
co2_opt = '_tbd_'
years <- c('2021', '2022')
Catms <- c(414.7, 417.2) # from NOAA
source('../Data/Soybean-BioCro_Parameters/soybean_parameters.R')
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')
source('../Data/Soybean-BioCro_Parameters/soybean_modules.R')
# Load LD11 data
load('../Data/Soybean-BioCro_Parameters/full_soybean_ld11.RData')
# Load 2022 TNC data
TNC.data <- read.csv('../Data/2022_Carb_data/2022_LD11_TNC_new.csv')

# Modify the modules
full_soybean_ld11$direct_modules <- steady_state_module_names
full_soybean_ld11$differential_modules <- derivative_module_names

# Set up basic properties for the solver
full_soybean_ld11$ode_solver <- solver
# Update UTR parameters
fitted.utr.params <- optim_params_conversion(optim_params_short)
names(fitted.utr.params) <- arg_names
parameters <-c(parameters, fitted.utr.params)[!duplicated(c(names(parameters), 
                                                            names(fitted.utr.params)), 
                                                          fromLast = TRUE)]

update_parameters = TRUE
if (update_parameters){
  updated.params <- full_soybean_ld11$parameters[names(full_soybean_ld11$parameters) %in% names(parameters)]
  full_soybean_ld11$parameters <-c(parameters, 
                                   updated.params)[!duplicated(c(names(parameters), 
                                                                 names(updated.params)), 
                                                               fromLast = TRUE)]
}else{
  full_soybean_ld11$parameters <- parameters
}

ExpBiomass <- list()
weather.afteremergence <- list()
soybean_optsolver <- list()
numrows <- vector()

SLA_AS_DRIVER <- TRUE

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

for (i in 1:length(years)){
  ExpBiomass[[i]] <- loadRData(paste0('../../energy-farm-biocro/soybean_ld11_biomass_', years[i],'/soybean_ld11_biomass_', years[i], '.RData')) 
  weather <- loadRData(paste0('../../energy-farm-biocro/weather_', years[i], '/weather', years[i], '_hourly.RData'))
  weather.supplement <- loadRData(paste0('../Data/Weather_data/weather', years[i], 'supplement.RData'))
  full_soybean_ld11$parameters$Catm <- Catms[i]
  
  # update initial values
  sub_frac <- 0.1           # substrate_fraction
  str_frac <- 1 - sub_frac  # structural_fraction
  seed_mass <- ExpBiomass[[i]]$initial_seed[1]
  j <- 2
  mass_t <- sum(ExpBiomass[[i]][j, c('leaf', 'stem', 'root')])
  leaf_frac <- ExpBiomass[[i]]$leaf[j]/mass_t
  stem_frac <- ExpBiomass[[i]]$stem[j]/mass_t
  root_frac <- ExpBiomass[[i]]$root[j]/mass_t
  cf <- optim_params_short[1]
  
  initial_state <- list(
    Leaf_respiration_loss = 0.0,
    Stem_respiration_loss = 0.0,
    Root_respiration_loss = 0.0,
    Pod_respiration_loss = 0.0,
    # senescence related initial state
    Leaf_senescence_loss = 0.0,
    Stem_senescence_loss = 0.0,
    Root_senescence_loss = 0.0,
    Pod_senescence_loss = 0.0,
    # Other variables
    TTc = 0.0,
    soil_water_content = 0.32,
    cws1=                     0.32,          # dimensionless, current water status, soil layer 1
    cws2=                     0.32,          # dimensionless, current water status, soil layer 2
    Rhizome =                 0.0000001,     # Mg / ha
    RhizomeLitter =           0,               # Mg / ha
    Leaf_substrate_carbon = sub_frac * seed_mass * leaf_frac / cf,
    Leaf_structural_carbon = str_frac * seed_mass * leaf_frac / cf,
    Stem_substrate_carbon = sub_frac * seed_mass * stem_frac / cf, 
    Stem_structural_carbon = str_frac * seed_mass * stem_frac / cf,
    Root_substrate_carbon =  sub_frac * seed_mass * root_frac / cf,
    Root_structural_carbon = str_frac * seed_mass * root_frac / cf,
    Pod_substrate_carbon = 1e-4,  
    Pod_structural_carbon = 9e-4)
  
  update_init = TRUE
  if (update_init){
    updated.init <- full_soybean_ld11$initial_values[names(full_soybean_ld11$initial_values) %in% names(initial_state)]
    full_soybean_ld11$initial_values <-c(initial_state, 
                                         updated.init)[!duplicated(c(names(initial_state), 
                                                                     names(updated.init)), 
                                                                   fromLast = TRUE)]
  }else{
    full_soybean_ld11$initial_values <- initial_state
  }
  
  first_data_time <- ExpBiomass[[i]]$time[1]
  weather.aftersowing <- weather[weather$time >= first_data_time, ]
  weather.aftersowing$DVI <- weather.supplement$DVI
  # start from emergence time
  weather.afteremergence[[i]] <- weather.aftersowing[-(1:which.min(abs(weather.aftersowing$DVI))),]
  
  if (SLA_AS_DRIVER) {
    # The experimental data indicates a non-monotonic dependence of SLA on
    # time. As a test, try including it as a driver instead of a parameter,
    # where the values are interpolated from the experimental ones.
    ExpBiomass[[i]]$SLA[1] <- ExpBiomass[[i]]$SLA[2]-
      (ExpBiomass[[i]]$time[2]-ExpBiomass[[i]]$time[1])* 
      (ExpBiomass[[i]]$SLA[3]-ExpBiomass[[i]]$SLA[2])/
      (ExpBiomass[[i]]$time[3]-ExpBiomass[[i]]$time[2])
    sla_func <- approxfun(
      ExpBiomass[[i]]$time,
      ExpBiomass[[i]]$SLA,
      rule = 2,
      na.rm = TRUE,
      method = 'linear'
    )
    weather.afteremergence[[i]]$iSp <- sla_func(weather.afteremergence[[i]]$time)
    
    full_soybean_ld11$parameters$iSp <- NULL
  }
  
  full_soybean_ld11$parameters$timestep <- 1
  full_soybean_ld11$parameters$time_zone_offset <- NULL
  full_soybean_ld11$parameters$Rd = 1.28
  
  soybean_optsolver[[i]] <- with(full_soybean_ld11, {partial_run_biocro(
    initial_values,
    parameters,
    weather.afteremergence[[i]],
    direct_modules,
    differential_modules,
    ode_solver,
    arg_names,
    verbose = FALSE
  )})
  # delete the first data point since it is not included in the simulation.
  ExpBiomass[[i]] <- ExpBiomass[[i]][-1, c('doy','leaf', 'stem', 'pod', 'root','leaf_litter','stem_litter')]
  names(ExpBiomass[[i]]) <- c('DOY', 'Leaf', 'Stem', 'Pod', 'Root','LeafLitter','StemLitter')
  numrows[i] <- nrow(weather.afteremergence[[i]])
}

# Optimization
# cost function
wts <- data.frame("Leaf" = 1, "Stem" = 1,"Pod" = 2, "Root" = 0.5, "Litter" = 0.5, "TNC" = 1, "Pod_start" = 100)
cost_func <- function(x){
  EF_utr_optim(x, soybean_optsolver, ExpBiomass, TNC.data, numrows, wts)
}

## testing
i <- 2
result <- soybean_optsolver[[i]](optim_params_conversion(optim_params_short))

library(lattice)

xyplot(data = result,
       substrate_transport_Leaf_to_Stem+
         Leaf_utilization_rate~
         time,
       auto.key = TRUE,
       main = as.character(i))

print(xyplot(data = result,
             Leaf_utilization_rate/Leaf_structural_carbon
             +Stem_utilization_rate/Stem_structural_carbon
             +Root_utilization_rate/Root_structural_carbon
             +Pod_utilization_rate/Pod
             ~time,
             auto.key = TRUE))
xyplot(data = result,
       Leaf+Stem+Root+Pod~time,
       auto.key = TRUE)
cost_func(optim_params_short)

xyplot(data = result,
       Leaf_substrate_carbon/Leaf+
         Stem_substrate_carbon/Stem+
         Root_substrate_carbon/Root+
         Pod_substrate_carbon/Pod~
         time,
       auto.key = TRUE)

# Parameter ranges
upperlim <- c(0.5,  # 1: carbon to mass factor [(Mg / ha) / (mol / m^2)]
              0.1, 0.1, 0.1, 1.0, # 2，3，4，5： utilization rate constant [/hr]
              0.5, 0.5, 0.5, 0.5, # 6，7, 8, 9： Km [/]
              0.8, # 10: respiration factor [/]
              0.5, 0.5, 5, # 11, 12, 13: substrate conductance [Mg / hr / [Mg / ha]^beta]
              0.1, 0.1, 0.1, # 14,15,16: senescence rate max, LSR
              10.0, 10.0, 2.0, # 17,18,19: senescence alpha, LSR [dimensionless]
              2.0, 2.0, 2.0,
              1.2, 2.2) # 20,21,22: senescence beta, LSR [/dvi]


lowerlim <- c(0.2,  # 1: carbon to mass factor [(Mg / ha) / (mol / m^2)]
              0.0, 0.0, 0.0, 0.1, # 2，3，4，5： utilization rate constant [/hr]
              0.005, 0.005, 0.005, 0.1, # 6，7, 8, 9： Km [mol / Mg]
              0.1, # 10: respiration factor [dimensionless]
              0.005, 0.005 ,0.01, # 11，12，13:substrate conductance [Mg / hr / [Mg / ha]^beta]
              0.0, 0.0, 0.0, # 14,15,16: senescence rate max, LSR
              4.0, 4.0, 1.0,# 17, 18, 19: senescence alpha, LSR [dimensionless]
              1.5, 1.5, 1.5,
              0.6, 1.8) # 20, 21, 22: senescence beta, LSR [/dvi]

rng.seed <- 1234 # seed for random number generator
set.seed(rng.seed)
# maximum number of iterations
max.iter <- 1000

# Call DEoptim function to run optimization
parVars <- c('optim_params_conversion', 'EF_utr_optim','soybean_optsolver','ExpBiomass','TNC.data', 'numrows','wts')
cl <- makeCluster(8)
clusterExport(cl, parVars,envir=environment())
sink(paste0('Optmization_EF_multiyear_output_', Sys.Date(), '.txt'))
optim_result <- DEoptim(fn=cost_func, lower=lowerlim, upper = upperlim, 
                        control=list(itermax=max.iter,parallelType=1,
                                     packages=c('BioCro', 'UTRSoybeanBML'),
                                     parVar=parVars))
optim_params_short = optim_result$optim$bestmem
print(optim_params_short)
sink()

