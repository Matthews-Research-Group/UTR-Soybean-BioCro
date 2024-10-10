# The work directory is hard coded since the previous one doesn't work well all the time. 
# Please change it according to your own working directory.
# Parameter loading and plotting are put into this R code to check the variables more easily.
# To be moved to different files once everything works out.
# Clear workspace
rm(list=ls())

library(BioCro)
library(UTRSoybeanBML)
library(DEoptim)

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Cost function
source('soybean_parameter_expansion.R')
source('SoyFACE_multiyear_utr_optim_tnc_constraints.R')
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')

# set year and CO2 level
co2_opt = '_ambient_' # '_ambient_' or '_co2_'

# years, sowing dates, and harvesting dates of growing seasons being fit to
year <- c('2002', '2004', '2005', '2006')
sow.date <- c(152, 149, 148, 148)
harv.date <- c(288, 289, 270, 270)
weather.growingseason <- list()

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
  emergence.idx <- which(weather$DVI>0)[1]
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather.growingseason[[i]] <- weather[emergence.idx:hd.ind,]
  
  ExpBiomass[[i]] <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass.csv'))
  colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CummulativeLitter")
  
  ExpBiomass.std[[i]] <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CummulativeLitter")
  
  # RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[3], 
  #                             "Root"=0.17*sum(ExpBiomass[[i]][5,2:4])) 
  RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[-c(1,2)], 
                              "Root"=0.17*sum(ExpBiomass[[i]][5,2:4])) 
  # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  numrows[i] <- nrow(weather.growingseason[[i]])
  invwts <- ExpBiomass.std[[i]]
  weights[[i]] <- log(1/(invwts[,2:ncol(invwts)]+1e-5))
  
  # load parameter files
  source('soybean_initial_values.R')
  initial_state$DVI <- NULL
  soybean_optsolver[[i]] <- partial_run_biocro(initial_state,
                                             parameters,
                                             weather.growingseason[[i]],
                                             steady_state_module_names,
                                             derivative_module_names,
                                             solver,
                                             arg_names,
                                             verbose = FALSE)
}

# Optimization
wts2 <- data.frame("Stem" = 1, "Leaf" = 1, "Pod" = 5, "Root" = 0.75, "CummulativeLitter" = 0.5)

# cost function
cost_func <- function(x){
  multiyear_BioCro_optim(x, soybean_optsolver[c(1,2)], ExpBiomass[c(1,2)], 
                         numrows[c(1,2)], weights[c(1,2)], wts2, RootVals[c(1,2)])
}

cost_func(optim_params_short_SoyFACE)

# Parameter ranges
upperlim <- c(# 0.35,  # 1: carbon to mass factor [(Mg / ha) / (mol / m^2)]
              0.1, 0.1, 0.1, 1.0, # 2，3，4，5： utilization rate constant [/hr]
              0.5, 0.5, 0.5, 0.5, # 6，7, 8, 9： Km [/]
              0.8, # 10: respiration factor [/]
              0.5, 0.5, 5, # 11, 12, 13: substrate conductance [Mg / hr / [Mg / ha]^beta]
              0.1, 0.1, 0.1, # 14,15,16: senescence rate max, LSR
              10.0, 10.0, 2.0, # 17,18,19: senescence alpha, LSR [dimensionless]
              2.0, 2.0, 2.0,
              1.0,
              1.2, 2.2) # 20,21,22: senescence beta, LSR [/dvi]
lowerlim <- c(# 0.25,  # 1: carbon to mass factor [(Mg / ha) / (mol / m^2)]
              0.0, 0.0, 0.0, 0.1, # 2，3，4，5： utilization rate constant [/hr]
              0.0, 0.0, 0.0, 0.0, # 6，7, 8, 9： Km [mol / Mg]
              0.1, # 10: respiration factor [dimensionless]
              0.0, 0.0 ,0.0, # 11，12，13:substrate conductance [Mg / hr / [Mg / ha]^beta] 0.005, 0.005 ,0.01
              0.0, 0.0, 0.0, # 14,15,16: senescence rate max, LSR
              0.0, 0.0, 0.0,# 17, 18, 19: senescence alpha, LSR [dimensionless]
              1.5, 1.5, 1.5,
              0.0,
              0.8, 1.8) # 20, 21, 22: senescence beta, LSR [/dvi]


rng.seed <- 1234 # seed for random number generator
set.seed(rng.seed)
# maximum number of iterations
max.iter <- 1000
# Call DEoptim function to run optimization
cl <- makeCluster(8)
parVars <- c('optim_params_conversion', 'multiyear_BioCro_optim','soybean_optsolver','ExpBiomass','numrows','weights','wts2','RootVals')
clusterExport(cl, parVars,envir=environment())
sink(paste0('Optmization_output_', Sys.Date(), '.txt'))
optim_result <- DEoptim(fn=cost_func, lower=lowerlim, upper = upperlim, 
                        control=list(itermax=max.iter,parallelType=1,
                                     packages=c('BioCro', 'UTRSoybeanBML'),
                                     parVar=parVars,
                                     cluster=cl))
optim_params_short = optim_result$optim$bestmem
print(optim_params_short)
sink()
sink()

