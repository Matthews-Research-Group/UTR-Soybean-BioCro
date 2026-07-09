library(BioCro)
library(ggplot2) # To ggplot functions
library(reshape2) # For melt function
library(data.table) # For first and last functions
library(gridExtra)
library(dplyr)
library(tidyr)

source('plot_funcs.R')
#source("~/GitHub/biocro/inst/extdata/outputs/soybean_optim_seed_3456/soybean2.R")
my_soybean_parameter = soybean2
#only need to change this name. 
opt_result_names = c('new_water_v151_rng1234','new_water_v151_rng2345','new_water_v151_rng3456')

for (opt_result_name in opt_result_names){
   CO2_status = TRUE
   multiyear_biomass_plot(CO2_status,my_soybean_parameter,opt_result_name)
   CO2_status = FALSE 
   multiyear_biomass_plot(CO2_status,my_soybean_parameter,opt_result_name)
}
