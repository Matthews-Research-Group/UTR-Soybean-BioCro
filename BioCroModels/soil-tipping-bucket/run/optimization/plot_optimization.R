library(BioCro)
library(BioCroWater)
library(ggplot2) # To ggplot functions
library(reshape2) # For melt function
library(data.table) # For first and last functions
library(gridExtra)

source('set_modules_for_soil_water.R')
source('plot_funcs.R')
#only need to change this name. 
opt_result_name = 'new_water_v119'
use_new_water    = TRUE 
predefine        = FALSE 

use_varyingSLA   = FALSE 
use_pre_opt_vars = TRUE 

if(use_pre_opt_vars){
  pre_arg_names <- c('alphaRoot','betaRoot')
  
  pre_par = readRDS(paste0('opt_results/opt_result_new_water_v108.rds'))
  pre_par = pre_par$optim$bestmem
  pre_par = pre_par[c(2,5)]
  print(pre_par)
}
#
# Cost function
#arg_names <- c('alphaLeaf','alphaRoot','alphaStem','betaLeaf','betaRoot','betaStem',
#                     'rateSeneLeaf','rateSeneStem','alphaSeneLeaf','betaSeneLeaf',
#                     'alphaSeneStem','betaSeneStem','alphaShell','betaShell',
#                     'grc_stem','grc_root',
#                     'mrc_leaf','mrc_root'
#                    )

#arg_names     <- c('alphaLeaf','alphaStem','betaLeaf','betaStem',
#                     'rateSeneLeaf','rateSeneStem','alphaSeneLeaf','betaSeneLeaf',
#                     'alphaSeneStem','betaSeneStem','alphaShell','betaShell',
#                     'grc_stem','grc_root',
#                     'mrc_leaf','mrc_root')

#new_par = read.table('parameters_old_water_opt_v10.txt')
new_par   = readRDS(paste0('opt_results/opt_result_',opt_result_name,'.rds'))
arg_names = names(new_par)
#new_par = new_par$optim$bestmem
#names(new_par) = arg_names
soybean_parameters0          = soybean$parameters
df=data.frame(arg_names,NEW=format(as.numeric(new_par),scientific=F),CTL=as.numeric(soybean_parameters0[arg_names]))
print(df)

#since BioCro seperate mrcs for stem and leaf, here we make mrc_stem = mrc_leaf
#soybean_parameters0$mrc_stem = as.numeric(new_par['mrc_leaf'])

years <- c('2002','2004','2005','2006')

# sowing and harvest DOYs for each growing season
dates <- data.frame("year" = c(2002, 2004:2006),"sow" = c(152,149,148,148), "harvest" = c(288, 289, 270, 270))

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
  soybean_parameters0$iSp  = 2.9
}

if(use_pre_opt_vars){
  soybean_parameters0[pre_arg_names] = pre_par 
}

if(use_varyingSLA){
  soybean_steadystate_modules0$specific_leaf_area = "BioCroWater:varying_SLA"
}

soybean_solver_params=soybean$ode_solver
soybean_solver_params$type="homemade_euler" #comment this out to use rkck54 solver

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
  results_CTL[[i]] <- run_biocro(
    soybean$initial_values,
    soybean$parameters,
    weather_growing_season,
    soybean$direct_modules,
    soybean$differential_modules,
    soybean$ode_solver
  )
  
  soybean_solver[[i]] <- partial_run_biocro(soybean_initial_state0, 
 					       soybean_parameters0, 
					       weather_growing_season,
                                               soybean_steadystate_modules0,
					       soybean_derivative_modules0, 
                                               soybean_solver_params,arg_names)
  
  soybean_solver_i = soybean_solver[[i]] 
  results[[i]] <- soybean_solver_i(new_par)
#  results[[i]] <- run_biocro(soybean_initial_state0, 
# 			                       soybean_parameters0, 
#			                       weather_growing_season,
#                             soybean_steadystate_modules0,
#			                       soybean_derivative_modules0, 
  results[[i]]$time = results[[i]]$doy + results[[i]]$hour/24
  print(c(max(results[[i]]$Leaf),max(results[[i]]$lai)))

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

#   doy = results_CTL[[i]]$doy
#   kLeaf   = results_CTL[[i]]$kLeaf
##get the first doys when these terms become non-zero
#   doy_leaf  = doy[kLeaf>0.01][1]
#
#   doy_new = results[[i]]$doy
#   kLeaf_new   = results[[i]]$kLeaf
##get the first doys when these terms become non-zero
#   doy_leaf_new  = doy_new[kLeaf_new>0.01][1]
##   print(c(doy_leaf-doy[1],doy_leaf_new-doy_new[1]))
#  if(yr==2005){
#    print(results_CTL[[i]]$kLeaf[doy==180]) 
#    print(results[[i]]$kLeaf[doy_new==180]) 
#  }
}

#saveRDS(results,paste0("rds_results/result_",opt_result_name,".rds"))
#stop()

plot_list = list()
plot_list_litter = list()
for (i in 1:length(years)){
	yr = years[i]
#	FigA <- plot_all_tissues_TWO(results_CTL[[i]],results[[i]], yr, ExpBiomass[[i]], ExpBiomass.std[[i]])
	FigA <- plot_all_tissues(results[[i]], yr, ExpBiomass[[i]], ExpBiomass.std[[i]])
	plot_list[[i]] = FigA

        
        results[[i]]$TotalLitter = results[[i]]$LeafLitter + results[[i]]$StemLitter

	FigB <- plot_litters(results[[i]],yr, ExpBiomass[[i]], ExpBiomass.std[[i]])
	plot_list_litter[[i]] = FigB
}

pdf_name = paste0("figs/fig_",opt_result_name,".pdf")
pdf(pdf_name,height = 8, width=8,bg='transparent')
grid.arrange(grobs = plot_list,nrow=2,ncol=2)
dev.off()

#litter biomass plot
pdf_name = paste0("figs/fig_",opt_result_name,"_litter.pdf")
pdf(pdf_name,height = 8, width=8,bg='transparent')
grid.arrange(grobs = plot_list_litter,nrow=2,ncol=2)
dev.off()
