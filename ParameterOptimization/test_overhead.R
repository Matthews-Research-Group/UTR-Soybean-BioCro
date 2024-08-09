rm(list=ls())

library(BioCro)
library(lhs)
library(parallel)
library(pbapply)
library(microbenchmark)

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source files
source('soybean_parameter_expansion.R')
source('EnergyFarm_multiyear_utr_optim.R')
co2_opt = '_tbd_'
years <- c('2022')
Catms <- c(417.2) # from NOAA
source('../Data/Soybean-BioCro_Parameters/soybean_parameters.R')
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')
source('../Data/Soybean-BioCro_Parameters/soybean_modules.R')
# Load LD11 data
load('../Data/Soybean-BioCro_Parameters/full_soybean_ld11.RData')
# Load 2022 TNC data
TNC.data <- read.csv('../Data/2022_Carb_data/2022_LD11_TNC_new.csv')
TNC.data$Leaf <- TNC.data$Leaf * 6
TNC.data$Stem <- TNC.data$Stem * 6

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

# cost function
wts <- data.frame("Leaf" = 1, "Stem" = 1,"Pod" = 2, "Root" = 0.5, "Litter" = 0.5, "TNC" = 1, "Pod_start" = 100)
cost_func <- function(x){
  EF_utr_optim(x, soybean_optsolver, ExpBiomass, TNC.data, numrows, wts)
}

# Setup parameters
n_samples <- 1000
n_params <- 24
n_cores <- max(1, detectCores() - 1)

# Generate sample data
set.seed(123)
test_data <- matrix(runif(n_samples * n_params), nrow = n_samples)

# Function to run without progress bar
run_without_progress <- function() {
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("cost_func", "EF_utr_optim", "optim_params_conversion", "soybean_optsolver", 
                      "numrows", "ExpBiomass", "wts", "TNC.data"))
  on.exit(stopCluster(cl))
  parApply(cl, test_data, 1, cost_func)
}

# Function to run with progress bar
run_with_progress <- function() {
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("cost_func", "EF_utr_optim", "optim_params_conversion", "soybean_optsolver", 
                      "numrows", "ExpBiomass", "wts", "TNC.data"))
  on.exit(stopCluster(cl))
  pbapply(test_data, 1, cost_func, cl = cl)
}

# Run the benchmark
benchmark_results <- microbenchmark(
  with_progress = run_with_progress(),
  without_progress = run_without_progress(),
  times = 1  # Number of repetitions
)

# Print results
print(benchmark_results)

# Calculate percentage slowdown
mean_with <- mean(benchmark_results$time[benchmark_results$expr == "with_progress"])
mean_without <- mean(benchmark_results$time[benchmark_results$expr == "without_progress"])
slowdown_percent <- (mean_with - mean_without) / mean_without * 100

cat(sprintf("Estimated slowdown due to progress bar: %.2f%%\n", slowdown_percent))

# Optionally, plot the results
if (require(ggplot2)) {
  plot <- ggplot(benchmark_results, aes(x = expr, y = time / 1e9)) +
    geom_boxplot() +
    labs(title = "Execution Time: With vs Without Progress Bar",
         x = "Method", y = "Time (seconds)") +
    theme_minimal()
  print(plot)
}
