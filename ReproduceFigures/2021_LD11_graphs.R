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
library(lattice)
library(reshape2)
library(ggplot2)

# Load saved data
load('../../energy-farm-biocro/soybean_ld11_development_2021/soybean_ld11_development_data_2021.RData')
load('../../energy-farm-biocro/soybean_ld11_biomass_2021/soybean_ld11_biomass_2021.RData')
load('../Data/Soybean-BioCro_Parameters/full_soybean_ld11.RData')
load('../../energy-farm-biocro/weather_2021/weather2021_hourly.RData')
load('../Data/Weather_data/weather2021supplement.RData')
# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Modify the modules
full_soybean_ld11$direct_modules <- steady_state_module_names
full_soybean_ld11$differential_modules <- derivative_module_names

# Set up basic properties for the solver
full_soybean_ld11$ode_solver <- solver

# optim_params_short[19] = 1.6
full_soybean_ld11$parameters$Rd = 1.28
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

full_soybean_ld11$parameters$Catm <- 414.7 # 2021 value from NOAA

# update initial values
sub_frac <- 0.1           # substrate_fraction
str_frac <- 1 - sub_frac  # structural_fraction
seed_mass <- soybean_ld11_biomass_2021$initial_seed[1]
i <- 2
mass_t <- sum(soybean_ld11_biomass_2021[i, c('leaf', 'stem', 'root')])
leaf_frac <- soybean_ld11_biomass_2021$leaf[i]/mass_t
stem_frac <- soybean_ld11_biomass_2021$stem[i]/mass_t
root_frac <- soybean_ld11_biomass_2021$root[i]/mass_t
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
  # Leaf = seed_mass * leaf_frac,
  # Stem = seed_mass * stem_frac,
  # Root = seed_mass * root_frac,
  # Pod = 1e-3 * cf, 
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
SLA_AS_DRIVER <- TRUE

sowing_time <- soybean_ld11_biomass_2021$time[1]
#first_data_time <- soybean_ld11_biomass_2021$time[2]
first_data_time <- sowing_time
idx_diff <- which(weather2021_hourly$time == first_data_time) - which(weather2021_hourly$time == sowing_time)

# Add DVI to weather file
weather2021.aftersowing <- weather2021_hourly[weather2021_hourly$time >= first_data_time, ]
# weather2021.supplement.dvi.rep <- rep(weather2021.supplement$DVI, each = 6)
# weather2021.supplement.dl.rep <- rep(weather2021.supplement$day_length, each = 6)
# 
# idx_diff <- idx_diff+6
weather2021.aftersowing$DVI <- weather2021.supplement$DVI
# weather2021.aftersowing$DVI <- weather2021.supplement.dvi.rep[-(1:idx_diff)]
# weather2021.aftersowing$day_length <- weather2021.supplement.dl.rep[-(1:idx_diff)]

# start from emergence time
weather2021.afteremergence <- weather2021.aftersowing[-(1:which.min(abs(weather2021.aftersowing$DVI))),]

if (SLA_AS_DRIVER) {
  # The experimental data indicates a non-monotonic dependence of SLA on
  # time. As a test, try including it as a driver instead of a parameter,
  # where the values are interpolated from the experimental ones.
  soybean_ld11_biomass_2021$SLA[1] <- soybean_ld11_biomass_2021$SLA[2]-
    (soybean_ld11_biomass_2021$time[2]-soybean_ld11_biomass_2021$time[1])* 
    (soybean_ld11_biomass_2021$SLA[3]-soybean_ld11_biomass_2021$SLA[2])/
    (soybean_ld11_biomass_2021$time[3]-soybean_ld11_biomass_2021$time[2])
  sla_func <- approxfun(
    soybean_ld11_biomass_2021$time,
    soybean_ld11_biomass_2021$SLA,
    rule = 2,
    na.rm = TRUE,
    method = 'linear'
  )
  
  weather2021.afteremergence$iSp <- sla_func(weather2021.afteremergence$time)
  
  full_soybean_ld11$parameters$iSp <- NULL
}

full_soybean_ld11$parameters$timestep <- 1
full_soybean_ld11$parameters$time_zone_offset <- NULL
# optim_params_short[19] = 1.6
full_soybean_ld11$parameters$Rd = 1.28

# Run the soybean simulation starting at noon on June 17 (DOY 168)
biocro_result <- with(full_soybean_ld11, {run_biocro(
    initial_values,
    parameters,
    weather2021.afteremergence,
    direct_modules,
    differential_modules,
    ode_solver
)})

# Plot the biomass along with the measured values
xyplot(
    Leaf + Stem + Root + Pod ~ time,
    data = biocro_result,
    type = 'l',
    auto = TRUE,
    grid = TRUE,
    xlab = 'Day of year (2021)',
    ylab = 'Soybean biomass (Mg / ha)',
    ylim = c(-1, 6),
    panel = function(...) {
        panel.xyplot(...)
        panel.points(
            soybean_ld11_biomass_2021$leaf ~ soybean_ld11_biomass_2021$time,
            type = 'b',
            col = 'darkblue',
            pch = 16,
            lty = 2
        )
        panel.points(
            soybean_ld11_biomass_2021$root ~ soybean_ld11_biomass_2021$time,
            type = 'b',
            col = 'darkgreen',
            pch = 16,
            lty = 2
        )
        panel.points(
            soybean_ld11_biomass_2021$stem ~ soybean_ld11_biomass_2021$time,
            type = 'b',
            col = 'darkmagenta',
            pch = 16,
            lty = 2
        )
        panel.points(
            soybean_ld11_biomass_2021$pod ~ soybean_ld11_biomass_2021$time,
            type = 'b',
            col = 'darkred',
            pch = 16,
            lty = 2
        )
    }
)


biocro_organ_biomass <- biocro_result[c('time', 'Leaf', 'Stem', 'Root', 'Pod')]
biocro_organ_biomass_tall <- melt(biocro_organ_biomass, id.vars = 'time')
names(biocro_organ_biomass_tall) <- c('time','Organ', 'biomass')

field_organ_biomass <- soybean_ld11_biomass_2021[c('time', 'leaf', 'stem', 'root', 'pod')]
names(field_organ_biomass)[names(field_organ_biomass) %in% c('leaf', 'stem', 'root', 'pod')] <- c('Leaf', 'Stem', 'Root', 'Pod')
field_organ_biomass_tall <- melt(field_organ_biomass, id.vars = 'time')
names(field_organ_biomass_tall) <- c('time','Organ', 'biomass')
biocro_organ_biomass_tall_2021_ld11 <- biocro_organ_biomass_tall
field_organ_biomass_tall_2021_ld11 <- field_organ_biomass_tall
save(biocro_organ_biomass_tall_2021_ld11, field_organ_biomass_tall_2021_ld11, file = 'organ_biomass_plot_2021_ld11.RData')


# size.title <- 16
# size.axislabel <-14
# size.axis <- 16
# size.legend <- 14
# 
# ggplot() + theme_classic() +
#   geom_line(data = biocro_organ_biomass_tall, aes(x = time, y = biomass, color = Organ), size = 1) +
#   geom_point(data = field_organ_biomass_tall, aes(x = time, y = biomass, color = Organ), shape = 15, size = 4)+
#   theme(plot.title=element_text(size=size.title, hjust=0.5),
#         axis.text=element_text(size=size.axis),
#         axis.title=element_text(size=size.axislabel),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent", colour = NA))+
#   scale_y_continuous(limits = c(0, 6), breaks = seq(0, 8, 2)) +
#   scale_x_continuous(breaks = seq(180,280,30))+
#   labs(title=element_blank(), x='Day of Year (2021)',y='Biomass (Mg/ha)')
# 
# # Plot the LAI along with the measured values
# xyplot(
#   lai ~ time,
#   data = biocro_result,
#   type = 'l',
#   auto = TRUE,
#   grid = TRUE,
#   xlab = 'Day of year (2021)',
#   ylab = 'Soybean LAI (dimensionless)',
#   panel = function(...) {
#     panel.xyplot(...)
#     panel.points(
#       soybean_ld11_biomass_2021$LAI ~ soybean_ld11_biomass_2021$time,
#       type = 'b',
#       col = 'darkblue',
#       pch = 16,
#       lty = 2
#     )
#   }
# )
# 
# # Plot the specific leaf area along with the measured values
# xyplot(
#   Sp ~ time,
#   data = biocro_result,
#   type = 'l',
#   auto = TRUE,
#   grid = TRUE,
#   xlab = 'Day of year (2021)',
#   ylab = 'Soybean specific leaf area (ha / Mg)',
#   ylim = c(0, 4),
#   panel = function(...) {
#     panel.xyplot(...)
#     panel.points(
#       soybean_ld11_biomass_2021$SLA ~ soybean_ld11_biomass_2021$time,
#       type = 'b',
#       col = 'darkblue',
#       pch = 16,
#       lty = 2
#     )
#   }
# )
# 
# # Plot the Leaf nonstructural C along with total nonstructural carbon
# biocro_result$total_substrate = biocro_result$Leaf_substrate_carbon + biocro_result$Stem_substrate_carbon + biocro_result$Root_substrate_carbon + biocro_result$Pod_substrate_carbon
# xyplot(
#   Sp ~ total_substrate,
#   data = biocro_result[1:which.max(biocro_result$Pod),],
#   type = 'l',
#   auto = TRUE,
#   grid = TRUE,
#   xlab = 'Soybean total substrate C (Mg / ha)',
#   ylab = 'Specific Leaf Area (ha / Mg)'
# )
# 
# 
# # Plot senescence
# s_key <- list(text = c('Leaf Litter',
#                        'Stem Litter'),
#               x = .05, y = .95)
# cf <- optim_params_short[1]
# biocro_senescence_plot <- xyplot(data = biocro_result,
#                                  Leaf_senescence_loss * cf
#                                  +Stem_senescence_loss * cf
#                                  ~time,
#                                  main = '2021 Senescence Litter', 
#                                  xlab = 'DOY', 
#                                  ylab='Litter (Mg/ha)',
#                                  ylim = c(0,2), auto.key=s_key)
# 
# experiment_plot <- xyplot(data = soybean_ld11_biomass_2021,
#                           stem_litter + leaf_litter ~time,
#                           cex = 2)
# 
# dev.new()
# biocro_senescence_plot+as.layer(experiment_plot)
# 
# weather2021.afteremergence$too_cold[which(weather2021.afteremergence$temp<13)]='Blue'
# weather2021.afteremergence$too_cold[which(weather2021.afteremergence$temp>=13)]='Black'
# weather2021.afteremergence$too_cold[which(weather2021.afteremergence$temp>=18)]='darkgreen'
# weather2021.afteremergence$too_cold[which(weather2021.afteremergence$temp>=30)]='Red'
# plot(data = weather2021.afteremergence, temp~doy, col = too_cold, xlim = c(200,280), ylim=c(5,35))
# 
# not_really_cold <- unique(weather2021.afteremergence$doy[weather2021.afteremergence$doy == unique(weather2021.afteremergence$doy[weather2021.afteremergence$temp<13]) & weather2021.afteremergence$temp > 18])
# night_cold <- unique(weather2021.afteremergence$doy[weather2021.afteremergence$too_cold=='Blue'])
# night_cold[!(night_cold %in% not_really_cold)]

source('plot_partitioning.R')
allocation_percentage_tall <- plot_partitioning(biocro_result, '2021')

