# The work directory is hard coded since the previous one doesn't work well all the time. 
# Please change it according to your own working directory.
# Parameter loading and plotting are put into this R code to check the variables more easily.
# To be moved to different files once everything works out.
# Clear workspace
rm(list=ls())

library(BioCro)
library(DEoptim)
library(ggplot2)

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Cost function
source('soybean_parameter_expansion.R')
source('SoyFACE_EF_optim.R')

# Source files
co2_opt = 'ambient'
source('../Data/Soybean-BioCro_Parameters/soybean_parameters.R')
source('../Data/Soybean-BioCro_Parameters/soybean_modules.R')

year <- c('2002', '2004', '2005', '2006')
sow.date <- c(152, 149, 148, 148)
harv.date <- c(288, 289, 270, 270)
weather.growingseason <- list()

# Load saved data
load('../../energy-farm-biocro/soybean_ld11_biomass_2022/soybean_ld11_biomass_2022.RData') 
load('../../energy-farm-biocro/ball_berry_curves_2021/soybean_ld11_ball_berry_parameters_2021.RData')
load('../../energy-farm-biocro/aci_curves_2021/soybean_ld11_fvcb_parameters_2021.RData')
load('../../energy-farm-biocro/weather_2022/weather2022_hourly.RData')
load('../Data/Soybean-BioCro_Parameters/full_soybean_ld11.RData') 
load('../Data/Weather_data/weather2022supplement.RData')
TNC.data <- read.csv('../Data/2022_Carb_data/2022_LD11_TNC.csv')

# Modify the modules
full_soybean_ld11$direct_modules <- steady_state_module_names
full_soybean_ld11$differential_modules <- derivative_module_names

# Set up basic properties for the solver
full_soybean_ld11$ode_solver <- solver
# Reset the Rd to correct the incorrect estimation
full_soybean_ld11$parameters$Rd = 1.28

optim_params_short <- c(0.3,  # 1: carbon to mass factor [(Mg / ha) / (mol / m^2)]
                        0.01, 0.015, 0.015, 0.5, # 2，3，4，5： utilization rate constant [/hr]
                        1, 0.75, 0.5, 0.25, # 6，7, 8, 9： Km [10^-4 mol / Mg]
                        0.2, # 10: respiration factor [dimensionless]
                        0.6, 0.4, 1, # 11，12，13: substrate conductance [Mg / m2 /hr / [Mg / ha]
                        0.01, 0.01, # 14, 15: 'Leaf_senescence_rate_max', 'Stem_senescence_rate_max'
                        11, 15, # 16, 17: senescence alpha [dimensionless]
                        -5.5, -7, # 18, 19: senescence beta [/dvi]
                        1, 2) # 20, 21: DVI switches

optim_params_short <-c(0.242485,    
                       0.004276,    0.032644,    0.053587,    0.586933,    
                       0.062188,    0.902778,    0.832786,    0.485468,    
                       0.357419,    
                       0.352924,    0.373693,    1.157740,
                       0.008441,    0.020550,   
                       16.118925,   20.463494,   -9.096414,   -9.807364,    
                       0.825496,    1.978349) 

fitted.thornley.params <- optim_params_conversion(optim_params_short)
arg_names <- c('Leaf_carbon_to_mass_factor', 'Stem_carbon_to_mass_factor', # 1, 2 
               'Root_carbon_to_mass_factor', 'Pod_carbon_to_mass_factor',  # 3, 4
               'Leaf_utilization_rate_constant', 'Stem_utilization_rate_constant', # 5, 6
               'Root_utilization_rate_constant', 'Pod_utilization_rate_constant',  # 7, 8
               'Leaf_utilization_km', 'Stem_utilization_km', 'Root_utilization_km', 'Pod_utilization_km', # 9,10,11,12
               'Leaf_respiration_factor', 'Stem_respiration_factor', # 13, 14
               'Root_respiration_factor', 'Pod_respiration_factor', # 15, 16
               'substrate_conductance_Leaf_to_Stem', 'substrate_conductance_Stem_to_Root', # 17, 18
               'substrate_conductance_Stem_to_Pod', # 19
               'Leaf_senescence_rate_max','Stem_senescence_rate_max', # 20, 21
               'Leaf_senescence_alpha', 'Stem_senescence_alpha',# 22, 23
               'Leaf_senescence_beta', 'Stem_senescence_beta', # 24, 25
               'Pod_start_dvi', 'stop_growth_dvi') # 26, 27
names(fitted.thornley.params) <- arg_names
# fitted.thornley.params$Stem_carbon_to_mass_factor = 0.15
parameters <-c(parameters, fitted.thornley.params)[!duplicated(c(names(parameters), 
                                                                 names(fitted.thornley.params)), 
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

# full_soybean_ld11[['2021']]$parameters$Catm <- 414.7 # 2021 value from NOAA
full_soybean_ld11$parameters$Catm <- 417.2 # 2022 value from NOAA


# update initial values
sub_frac <- 0.1           # substrate_fraction
str_frac <- 1 - sub_frac  # structural_fraction
seed_mass <- soybean_ld11_biomass_2022$initial_seed[1]
i <- 2
mass_t <- sum(soybean_ld11_biomass_2022[i, c('leaf', 'stem', 'root')])
leaf_frac <- soybean_ld11_biomass_2022$leaf[i]/mass_t # Mg/ha
stem_frac <- soybean_ld11_biomass_2022$stem[i]/mass_t # Mg/ha
root_frac <- soybean_ld11_biomass_2022$root[i]/mass_t # Mg/ha
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
  # Variables related to the utilization growth model starting from first datapoint
  # Biomass
  Leaf = seed_mass * leaf_frac,
  Stem = seed_mass * stem_frac,
  Root = seed_mass * root_frac,
  Pod = 1e-8 * cf, 
  Leaf_substrate_carbon = sub_frac * seed_mass * leaf_frac / cf,
  Leaf_structural_carbon = str_frac * seed_mass * leaf_frac / cf,
  Stem_substrate_carbon = sub_frac * seed_mass * stem_frac / cf, 
  Stem_structural_carbon = str_frac * seed_mass * stem_frac / cf,
  Root_substrate_carbon =  sub_frac * seed_mass * root_frac / cf,
  Root_structural_carbon = str_frac * seed_mass * root_frac / cf,
  Pod_substrate_carbon = 1e-9,  
  Pod_structural_carbon = 9e-9
)
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




# Make some decisions about what to do
SET_NEW_PARAMETER_VALUES <- TRUE
VERBOSE_MODEL_VALIDATION <- FALSE
SLA_AS_DRIVER <- TRUE

sowing_time <- soybean_ld11_biomass_2022$time[1]
first_data_time <- soybean_ld11_biomass_2022$time[1]
idx_diff <- which(weather2022_hourly$time == first_data_time) - which(weather2022_hourly$time == sowing_time)

# Add DVI to weather file
weather2022.aftersowing <- weather2022_hourly[weather2022_hourly$time >= first_data_time, ]
# weather2022.aftersowing$DVI <- weather2022.supplement$DVI[-(1:idx_diff)]
weather2022.aftersowing$DVI <- weather2022.supplement$DVI

# start from emergence time
weather2022.aftersowing <- weather2022.aftersowing[-(1:which.min(abs(weather2022.aftersowing$DVI))),]

if (SLA_AS_DRIVER) {
  # The experimental data indicates a non-monotonic dependence of SLA on
  # time. As a test, try including it as a driver instead of a parameter,
  # where the values are interpolated from the experimental ones.
  soybean_ld11_biomass_2022$SLA[1] <- soybean_ld11_biomass_2022$SLA[2]-
    (soybean_ld11_biomass_2022$time[2]-soybean_ld11_biomass_2022$time[1])* 
    (soybean_ld11_biomass_2022$SLA[3]-soybean_ld11_biomass_2022$SLA[2])/
    (soybean_ld11_biomass_2022$time[3]-soybean_ld11_biomass_2022$time[2])
  sla_func <- approxfun(
    soybean_ld11_biomass_2022$time,
    soybean_ld11_biomass_2022$SLA,
    rule = 2,
    na.rm = TRUE,
    method = 'linear'
  )
  
  weather2022.aftersowing$iSp <- sla_func(weather2022.aftersowing$time)
  
  full_soybean_ld11$parameters$iSp <- NULL
}

full_soybean_ld11$parameters$timestep <- 1

# Run the soybean simulation starting at noon on June 17 (DOY 168)
soybean_optsolver <- with(full_soybean_ld11, {partial_run_biocro(
  initial_values,
  parameters,
  weather2022.aftersowing,
  direct_modules,
  differential_modules,
  ode_solver,
  arg_names,
  verbose = FALSE
  )})

# Optimization
upperlim <- c(0.5,  # 1: carbon to mass factor [(Mg / ha) / (mol / m^2)]
              0.1, 0.1, 0.1, 1, # 2，3，4，5： utilization rate constant [/hr]
              2, 2, 2, 2, # 6，7, 8, 9： Km [mol / Mg]
              0.5, # 10: respiration factor [dimensionless]
              1, 1, 2, # 11, 12, 13: substrate conductance [Mg / hr / [Mg / ha]^beta]
              0.05, 0.05, # 14,15: 'Leaf_senescence_rate_max'
              30, 30, # 16, 17: senescence alpha [dimensionless]
              -3, -3, # 18, 19: senescence beta [/dvi]
              1.2, 2.2)# 20, 21: DVI switches


lowerlim <- c(0.2,  # 1: carbon to mass factor [(Mg / ha) / (mol / m^2)]
              0, 0, 0, 0, # 2，3，4，5： utilization rate constant [/hr]
              0, 0, 0, 0, # 6，7, 8, 9： Km [mol / Mg]
              0, # 10: respiration factor [dimensionless]
              0, 0 ,0, # 11，12，13:substrate conductance [Mg / hr / [Mg / ha]^beta]
              0,  0, # 14,15 'Leaf_senescence_rate_max'
              5, 5, # 16, 17: senescence alpha [dimensionless]
              -10, -10, # 18, 19: senescence beta [/dvi]
              0.8, 1.8)# 20, 21: DVI switches


# cost function
biomass <- soybean_ld11_biomass_2022[-1, c('doy','leaf', 'stem', 'pod', 'root','leaf_litter','stem_litter')]
names(biomass) <- c('DOY', 'Leaf', 'Stem', 'Pod', 'Root','LeafLitter','StemLitter')
numrows <- nrow(weather2022.aftersowing)
wts2 <- data.frame("Leaf" = 1, "Stem" = 1,"Pod" = 1, "Root" = 0.125, "Litter" = 0.125, "TNC" = 0.5)
cost_func <- function(x){
  thornley_2022_optim(x, soybean_optsolver, 
                      biomass,
                      TNC.data,
                      numrows, 
                      wts2)
}

rng.seed <- 1234 # seed for random number generator
set.seed(rng.seed)
# maximum number of iterations
max.iter <- 1000

## testing
result <- soybean_optsolver(optim_params_conversion(optim_params_short))

library(lattice)

xyplot(data = result,
       substrate_transport_Leaf_to_Stem+
         Leaf_utilization_rate~
         time,
       auto.key = TRUE,
       main = as.character(i))
library(latticeExtra)
biomass_plot <- xyplot(data = result, Leaf+Stem+Root+Pod~time, ylab = 'Biomass (Mg/ha)', auto.key = TRUE)
measured_biomass_plot <- xyplot(data = soybean_ld11_biomass_2022,
                                leaf + stem + root + pod~time,
                                type = 'l',
                                lwd = 3,
                                auto.key = TRUE)
biomass_plot + as.layer(measured_biomass_plot)

print(xyplot(data = result, Leaf_total_C_change_per_m2+
               Stem_total_C_change_per_m2+
               Root_total_C_change_per_m2+
               Pod_total_C_change_per_m2~time, auto.key = TRUE))
print(xyplot(data = result[1000:1100,], 
             Pod_substrate_carbon+
               Stem_substrate_carbon~
               time, auto.key = TRUE))
print(xyplot(data = result, 
             Leaf_substrate_carbon+
               Stem_substrate_carbon~
               time, auto.key = TRUE))
print(xyplot(data = result, 
             Leaf_structural_carbon+
               Stem_structural_carbon~
               time, auto.key = TRUE))
simulated.substrate <- xyplot(data = result[1:(which.min(abs(result$DVI-2))),],
             Leaf_substrate_carbon * 10^4 /Leaf+
             Stem_substrate_carbon * 10^4 /Stem+
             Root_substrate_carbon * 10^4 /Root+
             Pod_substrate_carbon * 10^4 /Pod
             ~time,
             ylab = 'Substrate Carbon (mol/Mg)',
             auto.key = TRUE)
measured.substrate <- xyplot(data = TNC.data,
                             Leaf + Stem ~ time,
                             type = 'l',
                             lwd = 3,
                             auto.key = TRUE)
print(simulated.substrate)
print(measured.substrate)
print(simulated.substrate + as.layer(measured.substrate))


print(xyplot(data = result,
             Leaf_structural_carbon * cf /Leaf+
             Stem_structural_carbon * cf /Stem+
             Root_structural_carbon * cf /Root+
             Pod_structural_carbon * cf /Pod~time,
             auto.key = TRUE))


print(xyplot(data = result[1:(which.min(abs(result$DVI-2))-5),],
             (Leaf_substrate_carbon+Leaf_structural_carbon) * cf / Leaf+
               (Stem_substrate_carbon+Stem_structural_carbon) * cf / Stem+
               (Pod_substrate_carbon+Pod_structural_carbon) * cf / Pod+
               (Root_substrate_carbon+Root_structural_carbon) * cf / Root~time,
             auto.key = TRUE))

print(xyplot(data = result[1:(which.min(abs(result$DVI-2))-5),],
             Leaf_utilization_rate/Leaf_structural_carbon
             +Stem_utilization_rate/Stem_structural_carbon
             +Root_utilization_rate/Root_structural_carbon
             # +Pod_utilization_rate/Pod
             ~time,
             auto.key = TRUE))
cost_func(optim_params_short)

# aggregate(10^4*result$Leaf_total_C_change_per_m2, by=list(Category=result$doy), FUN=sum)



# Call DEoptim function to run optimization
# parVars <- c('optim_params_conversion', 'thornley_2022_optim','soybean_optsolver','biomass','TNC.data', 'numrows','wts2')
# optim_result <- DEoptim(fn=cost_func, lower=lowerlim, upper = upperlim, control=list(itermax=max.iter,parallelType=1,packages=c('BioCro'),parVar=parVars))
# optim_params_short = optim_result$optim$bestmem
# print(optim_params_short)

