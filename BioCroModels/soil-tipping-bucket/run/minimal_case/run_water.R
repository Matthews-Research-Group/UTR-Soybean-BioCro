library(BioCro)
library(BioCroWater)
use_pre_opt_vars = TRUE
source('set_modules_for_soil_water.R')

opt_filename = "opt_result_new_water_v127.rds"
new_par = readRDS(opt_filename)
arg_names = names(new_par)

if(use_pre_opt_vars){
  pre_arg_names <- c('alphaRoot','betaRoot')
  
  pre_par = readRDS('opt_result_new_water_v108.rds')
  pre_par = pre_par$optim$bestmem
  pre_par = pre_par[c(2,5)] #'alphaRoot','betaRoot'
  names(pre_par) = pre_arg_names
  
  new_par = c(new_par,pre_par)
  arg_names = names(new_par)
}

soybean_steadystate_modules0 = set_direct_modules() 
soybean_derivative_modules0  = set_differential_modules() 
soybean_initial_state0       = set_init_values() 
soybean_parameters0          = set_parameters() 
soybean_parameters0$kd       = soybean_parameters0$k_diffuse

#define the gs weather for a year
weather_growing_season       = soybean_weather$`2002`

soybean_solver_params=soybean$ode_solver
#use euler to speed up the optimization process.
#you may use the non-euler in the evaluation process. 
#I have not seen any sig difference between solvers
soybean_solver_params$type="homemade_euler"

soybean_optsolver <- partial_run_biocro(soybean_initial_state0, 
 					       soybean_parameters0, 
					       weather_growing_season,
                                               soybean_steadystate_modules0,
					       soybean_derivative_modules0, 
                                               soybean_solver_params,arg_names)
results <- soybean_optsolver(new_par)
