# Clear the workspace
rm(list=ls())

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source files
source('../ParameterOptimization/soybean_parameter_expansion.R')
co2_opt = 'ambient'
source('../Data/Soybean-BioCro_Parameters/soybean_parameters.R')
source('../Data/Soybean-BioCro_Parameters/soybean_modules.R')

# Load packages
library(BioCro)
library(lattice)
library(latticeExtra)
library(reshape2)
library(ggplot2)


# Load saved data
load('../../energy-farm-biocro/soybean_ld10_biomass_2022/soybean_ld10_biomass_2022.RData') 
species.ball_berry.params <- read.csv('../../energy-farm-biocro/ball_berry_curves_2022/ball_berry_parameters_by_species.csv')
load('../../energy-farm-biocro/aci_curves_2022/c3_aci_averages_2022.RData')
load('../../energy-farm-biocro/weather_2022/weather2022_hourly.RData')
load('../Data/Soybean-BioCro_Parameters/full_soybean_ld11.RData') 
load('../Data/Weather_data/weather2022supplement.RData')

# Modify the modules
full_soybean_ld10 <- full_soybean_ld11
full_soybean_ld10$direct_modules <- steady_state_module_names
full_soybean_ld10$differential_modules <- derivative_module_names

# Set up basic properties for the solver
full_soybean_ld10$ode_solver <- solver
# solver$type <- 'boost_rkck54'
# solver$adaptive_max_steps <- 10000

optim_params_short <-c(0.233784,    0.004728,    0.034942,    0.022596,
                       0.962753,    0.164531,    1.840606,    1.357177,
                       3.092379,    0.293349,    0.317670,    0.779145,
                       1.198915,    0.013672,    0.039876,
                       10.891061,   14.012313,   -5.537149,   -5.798883,
                       0.883654,    1.995483)

# optim_params_short[19] = 1.6
full_soybean_ld10$parameters$Rd = 1.28
fitted.thornley.params <- optim_params_conversion(optim_params_short)
names(fitted.thornley.params) <- c('Leaf_carbon_to_mass_factor', 'Stem_carbon_to_mass_factor', # 1, 2 
                                   'Root_carbon_to_mass_factor', 'Pod_carbon_to_mass_factor',  # 3, 4
                                   'Leaf_utilization_rate_constant', 'Stem_utilization_rate_constant', # 5, 6
                                   'Root_utilization_rate_constant', 'Pod_utilization_rate_constant',  # 7, 8
                                   'Leaf_utilization_km', 'Stem_utilization_km', 'Root_utilization_km', 'Pod_utilization_km', # 9,10,11,12
                                   'Leaf_respiration_factor', 'Stem_respiration_factor', # 13, 14
                                   'Root_respiration_factor', 'Pod_respiration_factor', # 15, 16
                                   'substrate_conductance_Leaf_to_Stem', 'substrate_conductance_Stem_to_Root', # 17, 18
                                   'substrate_conductance_Stem_to_Pod', # 'transportation_beta_exponent', # 19, 20
                                   'Leaf_senescence_rate_max','Stem_senescence_rate_max', # 21, 22
                                   'Leaf_senescence_alpha', 'Stem_senescence_alpha',# 23, 24
                                   'Leaf_senescence_beta', 'Stem_senescence_beta', # 25, 26
                                   'Pod_start_dvi', 'stop_growth_dvi') # 27, 28
# fitted.thornley.params$Stem_carbon_to_mass_factor = 0.15
parameters <-c(parameters, fitted.thornley.params)[!duplicated(c(names(parameters), 
                                                                 names(fitted.thornley.params)), 
                                                               fromLast = TRUE)]

update_parameters = TRUE
if (update_parameters){
  updated.params <- full_soybean_ld10$parameters[names(full_soybean_ld10$parameters) %in% names(parameters)]
  full_soybean_ld10$parameters <-c(parameters, 
                                   updated.params)[!duplicated(c(names(parameters), 
                                                                 names(updated.params)), 
                                                               fromLast = TRUE)]
  # update ball-berry parameters
  full_soybean_ld10$parameters$b0 <- as.numeric(species.ball_berry.params$bb_intercept_avg[species.ball_berry.params$species=='ld10'][2])
  full_soybean_ld10$parameters$b1 <- as.numeric(species.ball_berry.params$bb_slope_avg[species.ball_berry.params$species=='ld10'][2])
  # update A-Ci parameters
  full_soybean_ld10$parameters$electrons_per_carboxylation <- 4.0
  full_soybean_ld10$parameters$electrons_per_oxygenation <- 4.0
  full_soybean_ld10$parameters$vmax1 <- 
    c3_aci_averages$main_data$Vcmax_at_25_avg[c3_aci_averages$main_data$species=='ld10']
  full_soybean_ld10$parameters$jmax <- 
    c3_aci_averages$main_data$J_at_25_avg[c3_aci_averages$main_data$species=='ld10']
  full_soybean_ld10$parameters$jmax <- 
    c3_aci_averages$main_data$J_at_25_avg[c3_aci_averages$main_data$species=='ld10']
  full_soybean_ld10$parameters$tpu_rate_max <- 
    c3_aci_averages$main_data$TPU_avg[c3_aci_averages$main_data$species=='ld10']
}else{
  full_soybean_ld10$parameters <- parameters
}

# full_soybean_ld10[['2021']]$parameters$Catm <- 414.7 # 2021 value from NOAA
full_soybean_ld10$parameters$Catm <- 417.2 # 2022 value from NOAA

# update initial values
sub_frac <- 0.1           # substrate_fraction
str_frac <- 1 - sub_frac  # structural_fraction
seed_mass <- soybean_ld10_biomass_2022$initial_seed[1]
i <- 2
mass_t <- sum(soybean_ld10_biomass_2022[i, c('leaf', 'stem', 'root')])
leaf_frac <- soybean_ld10_biomass_2022$leaf[i]/mass_t
stem_frac <- soybean_ld10_biomass_2022$stem[i]/mass_t
root_frac <- soybean_ld10_biomass_2022$root[i]/mass_t
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
  Pod = 1e-3 * cf, 
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
  updated.init <- full_soybean_ld10$initial_values[names(full_soybean_ld10$initial_values) %in% names(initial_state)]
  full_soybean_ld10$initial_values <-c(initial_state, 
                                       updated.init)[!duplicated(c(names(initial_state), 
                                                                   names(updated.init)), 
                                                                 fromLast = TRUE)]
}else{
  full_soybean_ld10$initial_values <- initial_state
}


# Make some decisions about what to do
MAKE_OPTIONAL_PLOTS <- TRUE
SET_NEW_PARAMETER_VALUES <- TRUE
VERBOSE_MODEL_VALIDATION <- FALSE
SLA_AS_DRIVER <- TRUE

sowing_time <- soybean_ld10_biomass_2022$time[1]
first_data_time <- soybean_ld10_biomass_2022$time[1]
idx_diff <- which(weather2022_hourly$time == first_data_time) - which(weather2022_hourly$time == sowing_time)

# Add DVI to weather file
weather2022.aftersowing <- weather2022_hourly[weather2022_hourly$time >= first_data_time, ]
# weather2022.aftersowing$DVI <- weather2022.supplement$DVI[-(1:idx_diff)]
weather2022.aftersowing$DVI <- weather2022.supplement$DVI

# start from emergence time
weather2022.aftersowing <- weather2022.aftersowing[-(1:which.min(abs(weather2022.aftersowing$DVI))),]
# weather2022.aftersowing <- weather2022.aftersowing[-(1:640),]

# if (SLA_AS_DRIVER) {
#   # The experimental data indicates a non-monotonic dependence of SLA on
#   # time. As a test, try including it as a driver instead of a parameter,
#   # where the values are interpolated from the experimental ones.
#   soybean_ld10_biomass_2022$SLA[1] = 0
#   soybean_ld10_biomass_2022$SLA[length(soybean_ld10_biomass_2022$SLA)] = 0
#   sla_func <- approxfun(
#     soybean_ld10_biomass_2022$time,
#     soybean_ld10_biomass_2022$SLA,
#     rule = 1,
#     # na.rm = TRUE,
#     method = 'linear'
#   )
#   
#   weather2021.aftersowing$iSp <- sla_func(weather2021.aftersowing$time) * 0.65
#   
#   full_soybean_ld10$parameters$iSp <- NULL
# }

if (SLA_AS_DRIVER) {
  # The experimental data indicates a non-monotonic dependence of SLA on
  # time. As a test, try including it as a driver instead of a parameter,
  # where the values are interpolated from the experimental ones.
  soybean_ld10_biomass_2022$SLA[1] <- soybean_ld10_biomass_2022$SLA[2]-
    (soybean_ld10_biomass_2022$time[2]-soybean_ld10_biomass_2022$time[1])* 
    (soybean_ld10_biomass_2022$SLA[3]-soybean_ld10_biomass_2022$SLA[2])/
    (soybean_ld10_biomass_2022$time[3]-soybean_ld10_biomass_2022$time[2])
  sla_func <- approxfun(
    soybean_ld10_biomass_2022$time,
    soybean_ld10_biomass_2022$SLA,
    rule = 2,
    na.rm = TRUE,
    method = 'linear'
  )
  
  weather2022.aftersowing$iSp <- sla_func(weather2022.aftersowing$time)
  
  full_soybean_ld10$parameters$iSp <- NULL
}

full_soybean_ld10$parameters$timestep <- 1
# Run the soybean simulation starting at noon on June 17 (DOY 168)
biocro_result <- with(full_soybean_ld10, {run_biocro(
    initial_values,
    parameters,
    weather2022.aftersowing,
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
    xlab = 'Day of year (2022)',
    ylab = 'Soybean biomass (Mg / ha)',
    ylim = c(-1, 8),
    panel = function(...) {
        panel.xyplot(...)
        panel.points(
            soybean_ld10_biomass_2022$leaf ~ soybean_ld10_biomass_2022$time,
            type = 'b',
            col = 'darkblue',
            pch = 16,
            lty = 2
        )
        panel.points(
            soybean_ld10_biomass_2022$root ~ soybean_ld10_biomass_2022$time,
            type = 'b',
            col = 'darkgreen',
            pch = 16,
            lty = 2
        )
        panel.points(
            soybean_ld10_biomass_2022$stem ~ soybean_ld10_biomass_2022$time,
            type = 'b',
            col = 'darkmagenta',
            pch = 16,
            lty = 2
        )
        panel.points(
            soybean_ld10_biomass_2022$pod ~ soybean_ld10_biomass_2022$time,
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

field_organ_biomass <- soybean_ld10_biomass_2022[c('time', 'leaf', 'stem', 'root', 'pod')]
names(field_organ_biomass)[names(field_organ_biomass) %in% c('leaf', 'stem', 'root', 'pod')] <- c('Leaf', 'Stem', 'Root', 'Pod')
field_organ_biomass_tall <- melt(field_organ_biomass, id.vars = 'time')
names(field_organ_biomass_tall) <- c('time','Organ', 'biomass')

# Save the data for plotting
biocro_organ_biomass_tall_2022_ld10 <- biocro_organ_biomass_tall
field_organ_biomass_tall_2022_ld10 <- field_organ_biomass_tall
save(biocro_organ_biomass_tall_2022_ld10, field_organ_biomass_tall_2022_ld10, file = 'organ_biomass_plot_2022_ld10.RData')


size.title <- 16
size.axislabel <-14
size.axis <- 16
size.legend <- 14

ggplot() + theme_classic() +
  geom_line(data = biocro_organ_biomass_tall, aes(x = time, y = biomass, color = Organ), size = 1) +
  geom_point(data = field_organ_biomass_tall, aes(x = time, y = biomass, color = Organ), shape = 15, size = 4)+
  theme(plot.title=element_text(size=size.title, hjust=0.5),
        axis.text=element_text(size=size.axis),
        axis.title=element_text(size=size.axislabel),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(), x='Day of Year (2022)',y='Biomass (Mg/ha)')
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
#       soybean_ld10_biomass_2022$LAI ~ soybean_ld10_biomass_2022$time,
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
#       soybean_ld10_biomass_2022$SLA ~ soybean_ld10_biomass_2022$time,
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

# Plot senescence
# s_key <- list(text = c('Leaf Litter',
#                        'Stem Litter'),
#               x = .05, y = .95)
# cf <- optim_params_short[1]
# biocro_senescence_plot <- xyplot(data = biocro_result,
#                                  Leaf_senescence_loss * cf
#                                  +Stem_senescence_loss * cf
#                                  ~time,
#                                  main = '2022 Senescence Litter', 
#                                  xlab = 'DOY', 
#                                  ylab='Litter (Mg/ha)',
#                                  ylim = c(0,2), auto.key=s_key)
# 
# experiment_plot <- xyplot(data = soybean_ld10_biomass_2022,
#                           stem_litter + leaf_litter ~time,
#                           cex = 2)
# 
# dev.new()
# biocro_senescence_plot+as.layer(experiment_plot)
# 
# weather2022.aftersowing$too_cold[which(weather2022.aftersowing$temp<13)]='Blue'
# weather2022.aftersowing$too_cold[which(weather2022.aftersowing$temp>=13)]='Black'
# weather2022.aftersowing$too_cold[which(weather2022.aftersowing$temp>=18)]='darkgreen'
# weather2022.aftersowing$too_cold[which(weather2022.aftersowing$temp>=30)]='Red'
# plot(data = weather2021.aftersowing, temp~doy, col = too_cold, xlim = c(200,280), ylim=c(5,35))
# 
# not_really_cold <- unique(weather2022.aftersowing$doy[weather2022.aftersowing$doy == unique(weather2022.aftersowing$doy[weather2022.aftersowing$temp<13]) & weather2022.aftersowing$temp > 18])
# night_cold <- unique(weather2021.aftersowing$doy[weather2022_hourly.aftersowing$too_cold=='Blue'])
# night_cold[!(night_cold %in% not_really_cold)]
save(biocro_organ_biomass_tall_2022_ld10, field_organ_biomass_tall_2022_ld10, file = 'organ_biomass_plot_2021_ld10.RData')

