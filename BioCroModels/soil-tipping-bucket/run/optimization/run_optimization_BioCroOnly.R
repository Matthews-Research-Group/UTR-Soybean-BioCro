library(BioCro)
library(DEoptim)
library(dplyr)
library(tidyr)

# Cost function
source('obj_function.R')

predefine        = FALSE 
use_parallel     = TRUE 

using_Gray_precip  = TRUE
obs_weather_Gray<-
    read.csv('../data/doi_10_5061_dryad_g0v62__v20170815/Final_Data_Deposit/ExtendedDataFig1_Weather_Data_2004thru2011/soyFACE_weather_data_2004thru2011.csv')
 
# Names of parameters being fit
arg_names_all <- c('alphaLeaf','alphaRoot','alphaStem','betaLeaf','betaRoot','betaStem',
                   'rateSeneLeaf','rateSeneStem','alphaSeneLeaf','betaSeneLeaf',
                   'alphaSeneStem','betaSeneStem','alphaShell','betaShell',
                   'grc_stem','grc_root',
                   'mrc_leaf','mrc_stem','mrc_root','mrc_grain','iSp')

arg_names     <- c('alphaLeaf','alphaRoot','alphaStem','betaLeaf','betaRoot','betaStem',
                   'rateSeneLeaf','rateSeneStem','alphaSeneLeaf','betaSeneLeaf',
                   'alphaSeneStem','betaSeneStem','alphaShell','betaShell',
                   'grc_stem','grc_root',
                   'mrc_leaf','mrc_root','iSp')

#arg_names     <- c('grc_stem','grc_root',
#                   'mrc_leaf','mrc_stem','mrc_root',
#                   'iSp')

#these weights are picked without specific reasons
#the current values work fine for my case
wts_on_errors <- data.frame("Stem" = 1, "Leaf" = 1, "Shell" = 0.5, "Seed" = 1,"TotalLitter" = 0.1,"Root" = 0.1,"LAI" = 1)

#fix a random seed for repeatability 
rng.seed <- 3456 
set.seed(rng.seed)

output_filename = paste0("opt_result_new_water_v150_rng",rng.seed,".rds")

print(arg_names)
print(wts_on_errors)
print(output_filename)

# years, sowing dates, and harvesting dates of growing seasons being fit to
# optimization uses these two years as Matthews et al (2022). 
# The other two years (2004&2006) can be considered as validation sets
year      <- c('2002', '2005')
sow.date  <- c(152, 148)
harv.date <- c(288, 270)
#year      <- c('2002', '2004')
#sow.date  <- c(152, 149)
#harv.date <- c(288, 289)

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

soybean_steadystate_modules0 = soybean2$direct_modules 
soybean_derivative_modules0  = soybean2$differential_modules 
soybean_initial_state0       = soybean2$initial_values 
soybean_parameters0          = soybean2$parameters 
soybean_solver_params        = soybean2$ode_solver

if(predefine){
  #pre-define
  #soybean_parameters0$mrc1 =  0.006 
  #soybean_parameters0$rateSeneStem =  0.001 
  #soybean_parameters0$phi2 =  0.5 
  #soybean_parameters0$electrons_per_carboxylation = 4
  #soybean_parameters0$electrons_per_oxygenation   = 4
  #soybean_parameters0$iSp = 3.0 
  soybean_parameters0$resistance_base     = 0.001196
  soybean_parameters0$resistance_amplifier= 6.452998
}

print(paste("iSp is",soybean_parameters0$iSp))
#print(paste("resistance are",soybean_parameters0$resistance_base,soybean_parameters0$resistance_amplifier))

for (i in 1:length(year)) {
  yr <- year[i]
  weather <- read.csv(file = paste0('Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))

  sd.ind <- which(weather$doy == sow.date[i])[1]
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather_growing_season <- weather[sd.ind:hd.ind,]
  weather_growing_season$time_zone_offset = -6

  if(using_Gray_precip  & yr>=2004){ #Gray's does not have year 2002
    obs_weather_Gray_yeari = obs_weather_Gray[obs_weather_Gray$Year==yr,]
    obs_weather_Gray_yeari = obs_weather_Gray_yeari[obs_weather_Gray_yeari$DOY>=sow.date[i] &obs_weather_Gray_yeari$DOY<=harv.date[i],]
    obs_weather_Gray_yeari_hourly <- obs_weather_Gray_yeari %>%
      # Add an hourly sequence per day
      uncount(weights = 24, .id = "hour") %>%
      # Adjust the hour (0 to 23)
      mutate(hour = hour - 1,
             precip = precip.mm. / 24)
    
    # replace precip with Gray's
#    weather_growing_season$precip = obs_weather_Gray_yeari_hourly$precip
    # replace precip with Gray's
    # If Gray's data is shorter than the growing season we hard-coded, we only replace the shorter part
    weather_growing_season$precip[weather_growing_season$doy>=obs_weather_Gray_yeari$DOY[1] &weather_growing_season$doy<=tail(obs_weather_Gray_yeari$DOY,1)] = obs_weather_Gray_yeari_hourly$precip
  }
  
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
            1e-2,1e-2,1e-2,1e-2,#mrc_leaf,mrc_stem,mrc_root,mrc_grain
            3.5 #iSp 
#YH: I lowered this not to create artificial bound. When I tested 3.5, the opt result came close to 3. To speed up optimization, I thus lowered this bound. Practically, we can simply fix this to reduce the number of parameters that need to be optimized
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


# initialize parameters being fitted as randome values from a uniform distribution
opt_pars <- runif(length(arg_names), min = lowerlim, max = upperlim)
# maximum number of iterations
max.iter <- 2000
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
