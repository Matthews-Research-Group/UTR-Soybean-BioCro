library(BioCro)
library(ggplot2)
library(grid)
library(gridExtra)
library(lattice)

# Clear workspace
rm(list=ls())

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Cost function
source('../ParameterOptimization/soybean_parameter_expansion.R')
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')

# initialize lists for figures
co2_opt = '_ambient_'

# names of fitted parameters
setwd('../Data/Soybean-BioCro_Parameters')
# load parameter files
param_files <- list.files(pattern = "[.]R$", recursive = TRUE)

param_files <- param_files[-c(which(param_files == "soybean_initial_values.R"))]
lapply(param_files, source)
# load parameter files
source('soybean_initial_values.R')

yr <- '2003'
weather <- read.csv(file = paste0('../Weather_data/', yr,'_Bondville_IL_daylength_wDVI.csv'))
emergence.ind <- which(weather$DVI>0)[1]
hd.ind <- which(weather$doy == 289)[24]
defoliation.ind <- which(weather$doy == 198 )[14]
initial_state$DVI <- NULL

weather.growingseason <- weather[emergence.ind:hd.ind,]
weather.growingseason_1 <- weather[emergence.ind:defoliation.ind,]
weather.growingseason_2 <- weather[defoliation.ind:hd.ind,]

ExpBiomass <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass.csv'))
colnames(ExpBiomass)<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

ExpBiomass.std <- read.csv(file=paste0('../SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
colnames(ExpBiomass.std)<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

RootVals <- data.frame("DOY"=ExpBiomass$DOY[3], "Root"=0.17*sum(ExpBiomass[5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130

numrows <- nrow(weather.growingseason)
invwts <- ExpBiomass.std

soybean_optsolver_no_hail <- partial_run_biocro(initial_state,
                                                parameters,
                                                weather.growingseason,
                                                steady_state_module_names,
                                                derivative_module_names,
                                                solver,
                                                arg_names,
                                                verbose = FALSE)
updated_params <- optim_params_conversion(optim_params_short_SoyFACE)
names(updated_params) <- arg_names
result_no_hail <- soybean_optsolver_no_hail(updated_params)

xyplot(data=result_no_hail, Leaf+Stem+Root+Pod~time, auto=TRUE)


soybean_optsolver_1 <- partial_run_biocro(initial_state,
                                          parameters,
                                          weather.growingseason_1,
                                          steady_state_module_names,
                                          derivative_module_names,
                                          solver,
                                          arg_names,
                                          verbose = FALSE)
# updated_params[28] <- 1.15
result_1 <- soybean_optsolver_1(optim_params_conversion(optim_params_short_SoyFACE))
# xyplot(data=result_1, Leaf+Stem+Root+Pod~time, auto=TRUE)
# Get the final values of the differential quantities; these will be the
# values just before defoliation
differential_quantities_just_before_defoliation <-
  as.list(result_1[nrow(result_1), names(initial_state)])

differential_quantities_just_after_defoliation <-
  differential_quantities_just_before_defoliation
# # Now reduce the leaf mass
remaining_leaf_percent <- 0.4
# differential_quantities_just_after_defoliation$Leaf <-
#   differential_quantities_just_before_defoliation$Leaf * remaining_leaf_percent

differential_quantities_just_after_defoliation$Leaf_substrate_carbon <-
  differential_quantities_just_before_defoliation$Leaf_substrate_carbon * remaining_leaf_percent
differential_quantities_just_after_defoliation$Leaf_structural_carbon <-
  differential_quantities_just_before_defoliation$Leaf_structural_carbon * remaining_leaf_percent # Could be changed to a different percentage

leaf_C_before_defoliation <- differential_quantities_just_before_defoliation$Leaf_substrate_carbon+
  differential_quantities_just_before_defoliation$Leaf_structural_carbon
stem_C_before_defoliation <- differential_quantities_just_before_defoliation$Stem_substrate_carbon+
  differential_quantities_just_before_defoliation$Stem_structural_carbon
abg_C_before_defoliation <- leaf_C_before_defoliation + stem_C_before_defoliation

leaf_C_after_defoliation <- differential_quantities_just_after_defoliation$Leaf_substrate_carbon+
  differential_quantities_just_after_defoliation$Leaf_structural_carbon

stem_new_percentage <- (abg_C_before_defoliation*(1-0.21) -
                          leaf_C_after_defoliation)/
                          stem_C_before_defoliation

stem_new_percentage <- 0.5

# If there is stem loss, then also reduce stem biomass.
if (stem_new_percentage < 1){
  # differential_quantities_just_after_defoliation$Stem <-
  #   differential_quantities_just_before_defoliation$Stem * stem_new_percentage
  differential_quantities_just_after_defoliation$Stem_substrate_carbon <-
    differential_quantities_just_before_defoliation$Stem_substrate_carbon * stem_new_percentage
  differential_quantities_just_after_defoliation$Stem_structural_carbon <-
    differential_quantities_just_before_defoliation$Stem_structural_carbon * stem_new_percentage
}

orignal_utr_params <- data.frame(optim_params_conversion(optim_params_short_SoyFACE))
rownames(orignal_utr_params) <- arg_names
colnames(orignal_utr_params) <- "Value"
parameters_after_hail <- orignal_utr_params
parameters_after_hail['Stem_respiration_factor', 'Value'] <- orignal_utr_params['Stem_respiration_factor', 'Value'] * 2
parameters_after_hail['Pod_start_dvi', 'Value'] <- orignal_utr_params['Pod_start_dvi', 'Value'] + 0.1

parameters$Leaf_respiration_factor <- orignal_utr_params['Stem_respiration_factor', 'Value']

soybean_optsolver_2 <- partial_run_biocro(differential_quantities_just_after_defoliation,
                                          parameters,
                                          weather.growingseason_2,
                                          steady_state_module_names,
                                          derivative_module_names,
                                          solver,
                                          arg_names,
                                          verbose = TRUE)


result_2 <- soybean_optsolver_2(parameters_after_hail$Value)
result <- rbind(result_1[seq_len(nrow(result_1) - 1), ], result_2)
xyplot(data=result, Leaf+Stem+Root+Pod~time, auto=TRUE)

# organize simulated data
r.lsrp.doy <- reshape2::melt(result[,c("time","Root","Leaf","Stem","Pod")],id.vars="time")
r.lsrp.doy$value<-r.lsrp.doy$value

# Leaf
# organize the expirimental data (mean and std)
s.exp.leaf <- cbind(ExpBiomass[,c("DOY","Leaf")])
colnames(s.exp.leaf) <- c("time","Leaf") # DOY renamed as time
r.exp.leaf <- reshape2::melt(s.exp.leaf, id.vars = "time") # DOY renamed as time

s.exp.std.leaf <- cbind(ExpBiomass.std[,c("DOY","Leaf")])
colnames(s.exp.std.leaf) <- c("time","Leaf") # DOY renamed as time
r.exp.std.leaf <- reshape2::melt(s.exp.std.leaf, id.vars = "time") # DOY renamed as time
r.exp.std.leaf$ymin <- r.exp.leaf$value - r.exp.std.leaf$value
r.exp.std.leaf$ymax <- r.exp.leaf$value + r.exp.std.leaf$value

# Stem
s.exp.stem <- cbind(ExpBiomass[,c("DOY","Stem")]) # DOY renamed as time
colnames(s.exp.stem) <- c("time","Stem") # DOY renamed as time
r.exp.stem <- reshape2::melt(s.exp.stem, id.vars = "time") # DOY renamed as time

s.exp.std.stem <- cbind(ExpBiomass.std[,c("DOY","Stem")]) # DOY renamed as time
colnames(s.exp.std.stem) <- c("time","Stem") # DOY renamed as time
r.exp.std.stem <- reshape2::melt(s.exp.std.stem, id.vars = "time") # DOY renamed as time
r.exp.std.stem$ymin <- r.exp.stem$value - r.exp.std.stem$value
r.exp.std.stem$ymax <- r.exp.stem$value + r.exp.std.stem$value

# Pod
s.exp.pod <- cbind(ExpBiomass[,c("DOY","Pod")]) # DOY renamed as time
colnames(s.exp.pod) <- c("time","Pod") # DOY renamed as time
r.exp.pod <- reshape2::melt(s.exp.pod, id.vars = "time") # DOY renamed as time

s.exp.std.pod <- cbind(ExpBiomass.std[,c("DOY","Pod")]) # DOY renamed as time
colnames(s.exp.std.pod) <- c("time","Pod") # DOY renamed as time
r.exp.std.pod <- reshape2::melt(s.exp.std.pod, id.vars = "time") # DOY renamed as time
r.exp.std.pod$ymin <- r.exp.pod$value - r.exp.std.pod$value
r.exp.std.pod$ymax <- r.exp.pod$value + r.exp.std.pod$value

# Combine
r.exp.ls <- rbind(r.exp.leaf, r.exp.stem, r.exp.pod)
r.exp.ls$Source = "Observed"
r.lsrp.doy$Source = "Simulated"
r.all <- rbind(r.lsrp.doy, r.exp.ls)
# Reverse the order as follow
r.all$Organ <- factor(r.all$variable, levels = rev(levels(r.all$variable)))

# combine the simulated and experimental data

# Colorblind friendly color palette (https://personal.sron.nl/~pault/)
col.palette.muted <- c( "#882255", "#999933", "#117733", "#332288")

size.title <- 14
size.axislabel <- 12
size.axis <- 12
size.legend <- 14

f <- ggplot() + theme_classic()

f <- f + geom_point(data=r.all, aes(x=time, y=value,
                                    color=Organ,
                                    size = Source, shape = Source),
                    show.legend = FALSE, stroke=0.5) +
  scale_shape_manual(values = c(15, 16)) +
  scale_size_manual(values = c(3, 0.5)) +
  scale_color_manual(values = col.palette.muted)

# for leaf
f <- f + geom_errorbar(data=r.exp.std.leaf, aes(x=time, ymin=ymin, ymax=ymax),  # DOY renamed as time
                       width=3.5, size=0.25, show.legend = FALSE)
# for stem
f <- f + geom_errorbar(data=r.exp.std.stem, aes(x=time, ymin=ymin, ymax=ymax),   # DOY renamed as time
                       width=3.5, size=0.25, show.legend = FALSE)

# for pod
f <- f + geom_errorbar(data=r.exp.std.pod, aes(x=time, ymin=ymin, ymax=ymax),   # DOY renamed as time
                       width=3.5, size=0.25, show.legend = FALSE)

# change the plot labels and theme
print(yr)
f <- f + labs(title=element_blank(), x=paste0('Day of Year (',yr,')'),y='Biomass (Mg / ha)')
f <- f + theme(plot.title=element_text(size=size.title, hjust=0.5),
               axis.text=element_text(size=size.axis),
               axis.title=element_text(size=size.axislabel),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA))

f <- f + scale_x_continuous(breaks = seq(150,280,30))
f <- f + scale_y_continuous(limits = c(0, 10), breaks = seq(0,10,2))
f

# xyplot(result_no_hail$substrate_transport_Leaf_to_Stem+
#          result$substrate_transport_Leaf_to_Stem~result$time,
#        auto.key = list(space = 'top'))
# xyplot(data = result[954:1024,], substrate_transport_Leaf_to_Stem~time)
# xyplot(data=result[954:1024,], Leaf_mass_fraction+Stem_mass_fraction~time, auto=TRUE)
# 
# canopy_assim_daily <- aggregate(result$canopy_assimilation_rate,list(result$doy), FUN=sum) * 10 / 3
# leaf_export_daily <- aggregate(result$substrate_transport_Leaf_to_Stem,list(result$doy), FUN=sum)
# leaf_utilization_daily <- aggregate(result$Leaf_utilization_rate,list(result$doy), FUN=sum)
# stem_utilization_daily <- aggregate(result$Stem_utilization_rate,list(result$doy), FUN=sum)
# df <- data.frame(DOY = leaf_export_daily$Group.1,
#                  Canopy_Assimilation_mol_per_day = canopy_assim_daily$x,
#                  Leaf_Export_mol_per_day = leaf_export_daily$x,
#                  Leaf_Utilization_mol_per_day = leaf_utilization_daily$x,
#                  Stem_Utilization_mol_per_day = stem_utilization_daily$x)
# xyplot(data = df,
#        Canopy_Assimilation_mol_per_day+
#        Leaf_Export_mol_per_day+
#        Leaf_Utilization_mol_per_day+
#        Stem_Utilization_mol_per_day~
#        DOY, type = c('p','l'), auto=TRUE)
# sum(result$substrate_transport_Leaf_to_Stem)
# sum(result$Stem_utilization_rate)
# 
# bottom_layer_daily_sunlit_Assim_no_hail <- aggregate(result_no_hail$sunlit_Assim_layer_9,
#                                              list(result$doy), FUN=sum)
# bottom_layer_daily_sunlit_Assim_w_hail <- aggregate(result$sunlit_Assim_layer_9,
#                                              list(result$doy), FUN=sum)
# xyplot(bottom_layer_daily_sunlit_Assim_no_hail$x +
#        bottom_layer_daily_sunlit_Assim_w_hail$x ~
#        bottom_layer_daily_sunlit_Assim_no_hail$Group.1,
#        type = c('p','l'),
#        main = 'Bottom Layer Sunlit',
#        xlab = 'DOY',
#        ylab = 'Net Assimilation (micromole / m^2 /s)',
#        auto.key = list(text = c("No Hail", "Hail")))
# 
# bottom_layer_daily_shaded_Assim_no_hail <- aggregate(result_no_hail$shaded_Assim_layer_9,
#                                                      list(result$doy), FUN=sum)
# bottom_layer_daily_shaded_Assim_w_hail <- aggregate(result$shaded_Assim_layer_9,
#                                                     list(result$doy), FUN=sum)
# xyplot(bottom_layer_daily_shaded_Assim_no_hail$x +
#          bottom_layer_daily_shaded_Assim_w_hail$x ~
#          bottom_layer_daily_shaded_Assim_w_hail$Group.1,
#        type = c('p','l'),
#        main = 'Bottom Layer Shaded',
#        xlab = 'DOY',
#        ylab = 'Net Assimilation (micromole / m^2 /s)',
#        auto.key = list(text = c("No Hail", "Hail")))
# 
# canopy_assim_daily_no_hail <- aggregate(result_no_hail$canopy_assimilation_rate,
#                                                      list(result$doy), FUN=sum)
# canopy_assim_daily_w_hail <- aggregate(result$canopy_assimilation_rate,
#                                                     list(result$doy), FUN=sum)
# 
# canopy_assim_plot <- xyplot(canopy_assim_daily_no_hail$x+
#                               canopy_assim_daily_w_hail$x
#                               ~canopy_assim_daily_w_hail$Group.1,
#                             type = c('p','l'),
#                             main = 'Canopy Net Assimilation',
#                             type = c('p','l'),
#                             xlab = 'DOY',
#                             type = c('p','l'),
#                             type = c('p','l'),
#                             auto.key = list(text = c("No Hail", "Hail")))
