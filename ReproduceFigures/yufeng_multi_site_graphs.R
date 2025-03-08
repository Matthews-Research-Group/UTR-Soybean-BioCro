# Clear the workspace
rm(list=ls())

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source files
source('../ParameterOptimization/soybean_parameter_expansion.R')
co2_opt = '_tbd_'
source('../Data/Soybean-BioCro_Parameters/soybean_parameters.R')
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')
source('../Data/Soybean-BioCro_Parameters/soybean_modules.R')
# Load packages
library(BioCro)
library(UTRSoybeanBML)
library(reshape2)
library(ggplot2)
library(patchwork)
library(lattice)

# Load saved data
sow_har_dates <- readRDS('../Data/yufeng_multi_site/sow_har_dates.rds')
load('../Data/Soybean-BioCro_Parameters/full_soybean_ld11.RData') 
load('../../energy-farm-biocro/soybean_ld11_biomass_2022/soybean_ld11_biomass_2022.RData')
weather_site <- list()
for (i in 1:2){
  file_name <- paste0('../Data/yufeng_multi_site/Weather_data/site_',i, '_2010_2022.csv')
  weather_site[[i]] <- read.csv(file_name)
}


plots <- list()
biocro_results <- list()
biocro_results_sites <- list()
# Modify the modules
full_soybean_ld11$direct_modules <- steady_state_module_names
full_soybean_ld11$differential_modules <- derivative_module_names

# Set up basic properties for the solver
full_soybean_ld11$ode_solver <- solver

# Set UTR parameters
fitted.utr.params <- optim_params_conversion(optim_params_short_SoyFACE)
names(fitted.utr.params) <- arg_names
parameters <-c(parameters, fitted.utr.params)[!duplicated(c(names(parameters), 
                                                                 names(fitted.utr.params)), 
                                                               fromLast = TRUE)]


updated.params <- full_soybean_ld11$parameters[names(full_soybean_ld11$parameters) %in% names(parameters)]
full_soybean_ld11$parameters <-c(parameters, 
                                   updated.params)[!duplicated(c(names(parameters), 
                                                                 names(updated.params)), 
                                                               fromLast = TRUE)]

# update initial values
sub_frac <- 0.1           # substrate_fraction
str_frac <- 1 - sub_frac  # structural_fraction
seed_mass <- soybean_ld11_biomass_2022$initial_seed[1]
i <- 2
mass_t <- sum(soybean_ld11_biomass_2022[i, c('leaf', 'stem', 'root')])
leaf_frac <- soybean_ld11_biomass_2022$leaf[i]/mass_t
stem_frac <- soybean_ld11_biomass_2022$stem[i]/mass_t
root_frac <- soybean_ld11_biomass_2022$root[i]/mass_t
cf <- 0.3 # optim_params_short[1]

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

# Make some decisions about what to do
MAKE_OPTIONAL_PLOTS <- TRUE
SET_NEW_PARAMETER_VALUES <- TRUE
VERBOSE_MODEL_VALIDATION <- FALSE
SLA_AS_DRIVER <- FALSE

for (site_num in 1:2){

  years <- 2010:2019
  
  if (site_num == 2){
    full_soybean_ld11$parameters$maturity_group <- 2
    full_soybean_ld11$parameters$lat <- 44.23944444
  }
  
  
  for(i in 1:length(years)){
    year <- years[i]
    full_soybean_ld11$parameters$Catm <- BioCro::catm_data[BioCro::catm_data$year==year,'Catm'] # 2010 value from NOAA
    sowing_time <- sow_har_dates[which(sow_har_dates$year == year), 'sowdate'][site_num]
    harvest_time <- sow_har_dates[which(sow_har_dates$year == year), 'hardate'][site_num]
    weather <- weather_site[[site_num]]
    weather.aftersowing <- weather[weather$year == year & 
                                     weather$doy >= sowing_time &
                                     weather$doy <= harvest_time, ]
    
    # Obtain DVI from the original Soybean-BioCro
    # change site 2 parameters
    soybean_parameters <- soybean$parameters
    if (site_num == 2){
      soybean_parameters$maturity_group <- 2
      soybean_parameters$lat <- 44.23944444
    }
    soybean_biocro_result <- run_biocro(soybean$initial_values,
                                        soybean_parameters,
                                        weather.aftersowing,
                                        soybean$direct_modules,
                                        soybean$differential_modules,
                                        soybean$ode_solver)
    
    weather.aftersowing$DVI <- soybean_biocro_result$DVI
    weather.aftersowing$time <- weather.aftersowing$doy + 
      weather.aftersowing$hour/24
    # SLA_AS_DRIVER <- FALSE
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
      
      weather.aftersowing$iSp <- sla_func(weather.aftersowing$time)
      
      full_soybean_ld11$parameters$iSp <- NULL
    }
    # print(xyplot(data=weather.aftersowing, iSp~time))
    
    full_soybean_ld11$parameters$timestep <- 1
    full_soybean_ld11$parameters$time_zone_offset <- NULL
    # optim_params_short[19] = 1.6
    full_soybean_ld11$parameters$Rd = 1.28
    full_soybean_ld11$parameters$iSp <- 3
    
    # Run the soybean simulation 
    soybean_optsolver <- with(full_soybean_ld11, {partial_run_biocro(
      initial_values,
      parameters,
      weather.aftersowing,
      direct_modules,
      differential_modules,
      ode_solver,
      arg_names,
      verbose = FALSE
    )})
    
    biocro_result <- soybean_optsolver(optim_params_conversion(optim_params_short_SoyFACE))
    
    # Plot the biomass along with the measured values
    biocro_organ_biomass <- biocro_result[c('time', 'Leaf', 'Stem', 'Root', 'Pod')]
    biocro_organ_biomass_tall <- melt(biocro_organ_biomass, id.vars = 'time')
    names(biocro_organ_biomass_tall) <- c('time','Organ', 'biomass')
    
    # Save the data for plotting
    biocro_organ_biomass_tall_2022_ld11 <- biocro_organ_biomass_tall
    
    size.title <- 12
    size.axislabel <-10
    size.axis <- 12
    size.legend <- 10
    
    col.palette.muted <- c( "#117733", "#999933",  "#882255", "#332288")
    
    p <- ggplot() + theme_classic() +
      geom_line(data = biocro_organ_biomass_tall, aes(x = time, y = biomass, color = Organ), linewidth = 1) +
      #   geom_point(data = field_organ_biomass_tall, aes(x = time, y = biomass, color = Organ), shape = 15, size = 4)+
      theme(plot.title=element_text(size=size.title, hjust=0.5),
            axis.text=element_text(size=size.axis),
            axis.title=element_text(size=size.axislabel),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.position = "none")+
      scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
      scale_x_continuous(breaks = seq(180,280,30))+
      labs(title=element_blank(), x=paste0('Site ', site_num, ' (', year, ')'), y='Biomass (Mg/ha)')+
      scale_color_manual(values = col.palette.muted)
    # print(p)
    plots[[i]] <- p
    biocro_results[[i]] <- biocro_result
    
  }
  # Combine all plots into a single plot
  combined_plot <- wrap_plots(plots, ncol = 4)  # Adjust ncol for number of columns
  
  # Display the combined plot
  print(combined_plot)
  saveRDS(biocro_results, paste0('site_',site_num,'_multi_year.rds'))
  # biocro_results_sites[[site_num]] <- biocro_results
}

i<-2

xyplot(biocro_results[[i]]$lai~
         biocro_results[[i]]$time,
       auto=TRUE)

