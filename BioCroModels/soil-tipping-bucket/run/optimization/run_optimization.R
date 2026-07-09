library(BioCro)
library(BioCroWater)
library(DEoptim)

# Cost function
source('obj_function.R')
source('set_modules_for_soil_water.R')

use_pre_opt_vars = TRUE 
predefine        = FALSE 
use_new_water    = TRUE 
use_parallel     = TRUE 
use_varyingSLA   = FALSE 
# Names of parameters being fit
arg_names_all <- c('alphaLeaf','alphaRoot','alphaStem','betaLeaf','betaRoot','betaStem',
                   'rateSeneLeaf','rateSeneStem','alphaSeneLeaf','betaSeneLeaf',
                   'alphaSeneStem','betaSeneStem','alphaShell','betaShell',
                   'grc_stem','grc_root',
                   'mrc_leaf','mrc_stem','mrc_root','mrc_grain','iSp'
                  )
arg_names     <- c('alphaLeaf','alphaStem','betaLeaf','betaStem',
                     'rateSeneLeaf','rateSeneStem','alphaSeneLeaf','betaSeneLeaf',
                     'alphaSeneStem','betaSeneStem','alphaShell','betaShell',
                     'grc_stem','grc_root',
                     'mrc_leaf','mrc_root','iSp')

if(use_pre_opt_vars){
  pre_arg_names <- c('alphaRoot','betaRoot')
  
  pre_par = readRDS(paste0('opt_results/opt_result_new_water_v108.rds'))
  pre_par = pre_par$optim$bestmem
  pre_par = pre_par[c(2,5)]
  print(pre_par)
}

#these weights are picked without specific reasons
#the current values work fine for my case
wts_on_errors <- data.frame("Stem" = 1, "Leaf" = 1, "Shell" = 0.5, "Seed" = 1,"TotalLitter" = 0.1,"Root" = 0.5,"LAI" = 1)

output_filename = "opt_result_new_water_v119.rds"

print(arg_names)
print(wts_on_errors)
print(output_filename)

# years, sowing dates, and harvesting dates of growing seasons being fit to
# optimization uses these two years as Matthews et al (2022). 
# The other two years (2004&2006) can be considered as validation sets
year      <- c('2002', '2005')
sow.date  <- c(152, 148)
harv.date <- c(288, 270)

obs_lai_years <- c(2002,2004:2010)
obs_lai <- readRDS('Data/obs_lai_2002_2010.rds')
#get the LAI for the calibration years
obs_lai_calib <- obs_lai[obs_lai_years %in% year]

## Initialize variables for the cost function
soybean_optsolver <- list()
ExpBiomass <- list()
ExpBiomass.std <- list()
RootVals <- list()
weights <- list()
numrows <- vector()

ctl <- list()

if(use_new_water){
  soybean_steadystate_modules0 = set_direct_modules() 
  soybean_derivative_modules0  = set_differential_modules() 
  soybean_initial_state0       = set_init_values() 
  soybean_parameters0          = set_parameters() 
  soybean_parameters0$kd = soybean_parameters0$k_diffuse
}else{
  soybean_steadystate_modules0 = soybean$direct_modules 
  soybean_derivative_modules0  = soybean$differential_modules 
  soybean_initial_state0       = soybean$initial_values 
  soybean_parameters0          = soybean$parameters 
}

if(predefine){
  #pre-define
  #soybean_parameters0$mrc1 =  0.006 
  #soybean_parameters0$rateSeneStem =  0.001 
  #soybean_parameters0$phi2 =  0.5 
  #soybean_parameters0$electrons_per_carboxylation = 4
  #soybean_parameters0$electrons_per_oxygenation   = 4
  soybean_parameters0$iSp = 2.9 
}

if(use_pre_opt_vars){
  soybean_parameters0[pre_arg_names] = pre_par 
}

if(use_varyingSLA){
  soybean_steadystate_modules0$specific_leaf_area = "BioCroWater:varying_SLA"
}


soybean_solver_params=soybean$ode_solver
#use euler to speed up the optimization process.
#you may use the non-euler in the evaluation process. 
#I have not seen any sig difference between solvers
soybean_solver_params$type="homemade_euler"

for (i in 1:length(year)) {
  yr <- year[i]
  weather <- read.csv(file = paste0('Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))

  sd.ind <- which(weather$doy == sow.date[i])[1]
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather_growing_season <- weather[sd.ind:hd.ind,]
  weather_growing_season$time_zone_offset = -6
  
#calculate CTL for diagnostic purposes
  ctl[[i]] <- run_biocro(
    soybean$initial_values,
    soybean$parameters,
    weather_growing_season,
    soybean$direct_modules,
    soybean$differential_modules,
    soybean$ode_solver
  )

  soybean_optsolver[[i]] <- partial_run_biocro(soybean_initial_state0, 
 					       soybean_parameters0, 
					       weather_growing_season,
                                               soybean_steadystate_modules0,
					       soybean_derivative_modules0, 
                                               soybean_solver_params,arg_names)
  
  ExpBiomass[[i]] <- read.csv(file=paste0('Data/biomasses_with_seed/',yr,'_ambient_biomass.csv'))
  colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem","Shell0","Seed","Litter","CumLitter")
  Shell = ExpBiomass[[i]]$Shell0 - ExpBiomass[[i]]$Seed
#make Shell not decline. When shell reaches peak, we assume it does not decline after 
#because in BioCro, we do not have shell senescence 
#  Shell[which.max(Shell):length(Shell)] = max(Shell) 
  ExpBiomass[[i]]$Shell = Shell 

  ExpBiomass[[i]]$Shell0 = NULL
  
  ExpBiomass.std[[i]] <- read.csv(file=paste0('Data/biomasses_with_seed/',yr,'_ambient_biomass_std.csv'))
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Shell0","Seed","Litter","CumLitter")
  ExpBiomass.std[[i]]$Shell = sqrt((ExpBiomass.std[[i]]$Shell0)^2 + (ExpBiomass.std[[i]]$Seed)^2)
  ExpBiomass.std[[i]]$Shell0 = NULL
#print(ExpBiomass[[i]]$Leaf+ExpBiomass[[i]]$Stem+ExpBiomass[[i]]$Seed+ExpBiomass[[i]]$Shell)  
#why use the 5th DOY?
  sum_shoot  = rowSums(ExpBiomass[[i]][,c('Leaf','Stem','Seed','Shell')])
#  sum_shoot  = rowSums(ExpBiomass[[i]][,c('Leaf','Stem')])
  row_at_max = which.max(sum_shoot)
  Root_estimate = 0.17 * sum_shoot[row_at_max]
  RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[row_at_max], "Root"=Root_estimate) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  print(RootVals[[i]])
  
  numrows[i] <- nrow(weather_growing_season)
  invwts <- ExpBiomass.std[[i]]
  #this weight is to put importance based on std
  #the larger the error bar is, the less the weight it has
  weights[[i]] <- log(1/(invwts[,2:ncol(invwts)]+1e-5))

}

## Optimization settings
ul = 50 
ll = -50

# parameter upper limit
upperlim_all<-c(ul,ul,ul,
            0,0,0,
            0.0125,.005,
            ul,0,ul,0,
            ul,0,
            0.08,0.075,#grc_stem,grc_root
            1e-2,1e-2,1e-2,1e-2, #mrc_leaf,mrc_stem,mrc_root,mrc_grain
            3.0 #iSp
            ) 

# parameter lower limit
lowerlim_all<-c(0,0,0,
            ll,ll,ll,
            0,0,
            0,ll,0,ll,
            0,ll,
            8e-4,0.0025,
            1e-5,1e-5,1e-5,1e-5,
            2.5
            )

upperlim = upperlim_all[match(arg_names,arg_names_all)]
lowerlim = lowerlim_all[match(arg_names,arg_names_all)]

if(length(arg_names) != length(upperlim)) stop("number of vars does not match the number of bounds")

# cost function
cost_func <- function(x){
  multiyear_BioCro_optim_obj(x, soybean_optsolver, ExpBiomass, numrows, weights, wts_on_errors, RootVals,obs_lai_calib)
}

#fix a random seed for repeatability 
rng.seed <- 1234 
set.seed(rng.seed)

# initialize parameters being fitted as randome values from a uniform distribution
opt_pars <- runif(length(arg_names), min = lowerlim, max = upperlim)
# maximum number of iterations
max.iter <- 2500
if(use_parallel){
  # number of cores for parallelization
  # you probably should pick a number less than the total cores on your machine
  nc = 8
  cl <- makeCluster(nc)
  
  ## Call DEoptim function to run optimization
  parVars <- c('multiyear_BioCro_optim_obj','soybean_optsolver','ExpBiomass','numrows','weights','wts_on_errors','RootVals','calc_cost','obs_lai_calib')
  clusterExport(cl, parVars,envir=environment())
  optim_result<-DEoptim(fn=cost_func, lower=lowerlim, upper = upperlim, 
		     control=list(itermax=max.iter,parallelType=1,
                     packages=c('BioCro'),parVar=parVars,cluster=cl))
}else{
  optim_result<-DEoptim(fn=cost_func, lower=lowerlim, upper = upperlim, 
       		     control=list(itermax=max.iter,parallelType=0,
                      packages=c('BioCro')))
}

#save the optimized parameters to a file, which will be used for plotting later
new_par = optim_result$optim$bestmem
names(new_par) = arg_names
saveRDS(new_par, paste0('opt_results/',output_filename))
