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
library(latticeExtra)
library(reshape2)
library(ggplot2)

# Load saved data
load('../../energy-farm-biocro/soybean_ld11_biomass_2022/soybean_ld11_biomass_2022.RData')
load('../../energy-farm-biocro/weather_2022/weather2022_hourly.RData')
load('../Data/Weather_data/weather2022supplement.RData')

load('../Data/Soybean-BioCro_Parameters/full_soybean_ld11.RData') 
leaf.carbs <- read.csv('../Data/2022_Carb_data/LD11_Leaf_Carbs.csv')
TNC.data <- read.csv('../Data/2022_Carb_data/2022_LD11_TNC_new.csv')
leaf.tnc.all <- read.csv('../Data/2022_Carb_data/Leaf_all.csv')
stem.tnc.all <- read.csv('../Data/2022_Carb_data/Stem_all.csv')

valid_data_long <- function(tnc.all){
  tnc.long <- data.frame()
  col.names <- colnames(tnc.all)[1:7]
  for (i in 1:nrow(tnc.all)){
    for (j in 1:4){
      if (tnc.all[i, paste0('Valid_',j)]==1){
        new_row <- c(tnc.all[i, col.names], TNC = tnc.all[i, paste0('Carb_tot_',j)])
        tnc.long <- rbind(tnc.long, new_row)
      }
    }
  }
  return (tnc.long)
}
leaf.tnc.long <- valid_data_long(leaf.tnc.all)
stem.tnc.long <- valid_data_long(stem.tnc.all)

# Modify the modules
full_soybean_ld11$direct_modules <- steady_state_module_names
full_soybean_ld11$differential_modules <- derivative_module_names

# Set up basic properties for the solver
full_soybean_ld11$ode_solver <- solver

# Set UTR parameters
fitted.utr.params <- optim_params_conversion(optim_params_short)
names(fitted.utr.params) <- arg_names
parameters <-c(parameters, fitted.utr.params)[!duplicated(c(names(parameters), 
                                                                 names(fitted.utr.params)), 
                                                               fromLast = TRUE)]


  updated.params <- full_soybean_ld11$parameters[names(full_soybean_ld11$parameters) %in% names(parameters)]
  full_soybean_ld11$parameters <-c(parameters, 
                                   updated.params)[!duplicated(c(names(parameters), 
                                                                 names(updated.params)), 
                                                               fromLast = TRUE)]
  full_soybean_ld11$parameters$Catm <- 417.2 # 2022 value from NOAA


# update initial values
sub_frac <- 0.1           # substrate_fraction
str_frac <- 1 - sub_frac  # structural_fraction
seed_mass <- soybean_ld11_biomass_2022$initial_seed[1]
i <- 2
mass_t <- sum(soybean_ld11_biomass_2022[i, c('leaf', 'stem', 'root')])
leaf_frac <- soybean_ld11_biomass_2022$leaf[i]/mass_t
stem_frac <- soybean_ld11_biomass_2022$stem[i]/mass_t
root_frac <- soybean_ld11_biomass_2022$root[i]/mass_t
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

sowing_time <- soybean_ld11_biomass_2022$time[1]
first_data_time <- soybean_ld11_biomass_2022$time[1]
idx_diff <- which(weather2022_hourly$time == first_data_time) - which(weather2022_hourly$time == sowing_time)

# Add DVI to weather file
weather2022.aftersowing <- weather2022_hourly[weather2022_hourly$time >= first_data_time, ]
# weather2022.aftersowing$DVI <- weather2022.supplement$DVI[-(1:idx_diff)]
weather2022.aftersowing$DVI <- weather2022.supplement$DVI

# start from emergence time
weather2022.afteremergence <- weather2022.aftersowing[-(1:which.min(abs(weather2022.aftersowing$DVI))),]

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
  
  weather2022.afteremergence$iSp <- sla_func(weather2022.afteremergence$time)
  
  full_soybean_ld11$parameters$iSp <- NULL
}

full_soybean_ld11$parameters$timestep <- 1
full_soybean_ld11$parameters$time_zone_offset <- NULL
# optim_params_short[19] = 1.6
full_soybean_ld11$parameters$Rd = 1.28

# Run the soybean simulation starting at noon on June 17 (DOY 168)
soybean_optsolver <- with(full_soybean_ld11, {partial_run_biocro(
  initial_values,
  parameters,
  weather2022.afteremergence,
  direct_modules,
  differential_modules,
  ode_solver,
  arg_names,
  verbose = FALSE
)})

biocro_result <- soybean_optsolver(optim_params_conversion(optim_params_short))

# Run the simulation with reduced light intensity by 20%
weather2022.afteremergence.lowlight <- weather2022.afteremergence
weather2022.afteremergence.lowlight$solar <- weather2022.afteremergence.lowlight$solar * 0.8
biocro_result_lowlight <- with(full_soybean_ld11, {run_biocro(
  initial_values,
  parameters,
  weather2022.afteremergence.lowlight,
  direct_modules,
  differential_modules,
  ode_solver
)})

# Run the simulation with increased temperature by 1 degree Celcius
weather2022.afteremergence.highT <- weather2022.afteremergence
weather2022.afteremergence.highT$temp <- weather2022.afteremergence.highT$temp + 1
biocro_result_highT <- with(full_soybean_ld11, {run_biocro(
  initial_values,
  parameters,
  weather2022.afteremergence.highT,
  direct_modules,
  differential_modules,
  ode_solver
)})

# compare the biomass reduction
print((max(biocro_result_lowlight$Leaf)-max(biocro_result$Leaf))/max(biocro_result$Leaf))
print((max(biocro_result_lowlight$Stem)-max(biocro_result$Stem))/max(biocro_result$Stem))
print((max(biocro_result_lowlight$Root)-max(biocro_result$Root))/max(biocro_result$Root))
print((max(biocro_result_lowlight$Pod)-max(biocro_result$Pod))/max(biocro_result$Pod))

# Plot the biomass along with the measured values
biocro_organ_biomass <- biocro_result[c('time', 'Leaf', 'Stem', 'Root', 'Pod')]
biocro_organ_biomass_tall <- melt(biocro_organ_biomass, id.vars = 'time')
names(biocro_organ_biomass_tall) <- c('time','Organ', 'biomass')

biocro_organ_lowlight_biomass <- biocro_result_lowlight[c('time', 'Leaf', 'Stem', 'Root', 'Pod')]
biocro_organ_biomass_lowlight_tall <- melt(biocro_organ_lowlight_biomass, id.vars = 'time')
names(biocro_organ_biomass_lowlight_tall) <- c('time','Organ', 'biomass')

field_organ_biomass <- soybean_ld11_biomass_2022[c('time', 'leaf', 'stem', 'root', 'pod')]
names(field_organ_biomass)[names(field_organ_biomass) %in% c('leaf', 'stem', 'root', 'pod')] <- c('Leaf', 'Stem', 'Root', 'Pod')
field_organ_biomass_tall <- melt(field_organ_biomass, id.vars = 'time')
names(field_organ_biomass_tall) <- c('time','Organ', 'biomass')


# Save the data for plotting
biocro_organ_biomass_tall_2022_ld11 <- biocro_organ_biomass_tall
field_organ_biomass_tall_2022_ld11 <- field_organ_biomass_tall
save(biocro_organ_biomass_tall_2022_ld11, field_organ_biomass_tall_2022_ld11, file = 'organ_biomass_plot_2022_ld11.RData')

size.title <- 16
size.axislabel <-14
size.axis <- 16
size.legend <- 14

col.palette.muted <- c( "#117733", "#999933",  "#882255", "#332288")

ggplot() + theme_classic() +
  geom_line(data = biocro_organ_biomass_tall, aes(x = time, y = biomass, color = Organ), size = 1) +
  geom_line(data = biocro_organ_biomass_lowlight_tall, aes(x = time, y = biomass, color = Organ), size = 1, linetype = 'dashed') +
  geom_point(data = field_organ_biomass_tall, aes(x = time, y = biomass, color = Organ), shape = 15, size = 4)+
  theme(plot.title=element_text(size=size.title, hjust=0.5),
        axis.text=element_text(size=size.axis),
        axis.title=element_text(size=size.axislabel),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_y_continuous(limits = c(0, 7), breaks = seq(0, 7, 2)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(), x='Day of Year (2022)',y='Biomass (Mg/ha)')+
  scale_color_manual(values = col.palette.muted)

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
#       soybean_ld11_biomass_2022$LAI ~ soybean_ld11_biomass_2022$time,
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
#       soybean_ld11_biomass_2022$SLA ~ soybean_ld11_biomass_2022$time,
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
# experiment_plot <- xyplot(data = soybean_ld11_biomass_2022,
#                           stem_litter + leaf_litter ~time,
#                           cex = 2)
# 
# dev.new()
# biocro_senescence_plot+as.layer(experiment_plot)
# 
# weather2022.afteremergence$too_cold[which(weather2022.afteremergence$temp<13)]='Blue'
# weather2022.afteremergence$too_cold[which(weather2022.afteremergence$temp>=13)]='Black'
# weather2022.afteremergence$too_cold[which(weather2022.afteremergence$temp>=18)]='darkgreen'
# weather2022.afteremergence$too_cold[which(weather2022.afteremergence$temp>=30)]='Red'
# plot(data = weather2021.aftersowing, temp~doy, col = too_cold, xlim = c(200,280), ylim=c(5,35))
# 
# not_really_cold <- unique(weather2022.afteremergence$doy[weather2022.afteremergence$doy == unique(weather2022.afteremergence$doy[weather2022.afteremergence$temp<13]) & weather2022.afteremergence$temp > 18])
# night_cold <- unique(weather2021.aftersowing$doy[weather2022_hourly.aftersowing$too_cold=='Blue'])
# night_cold[!(night_cold %in% not_really_cold)]

# # plot the Carbohydrates
leaf.carbs$carbs_total <- leaf.carbs$Starch + leaf.carbs$GFS_total
library(tidyverse)
# Reshape your data to be long
df_long <- leaf.carbs %>% 
  pivot_longer(cols = c("Starch", "GFS_total"), 
               names_to = "Measurement", 
               values_to = "Value")

df_std <- leaf.carbs %>% 
  pivot_longer(cols = c("Starch_std", "GFS_std"), 
               names_to = "Measurement", 
               values_to = "Std")

# Join the two reshaped data frames
df_joined <- data.frame(df_long$Time, df_long$DOY_date, df_long$hour, 
                        df_long$Measurement, df_long$Value, df_long$carbs_total,
                        df_std$Std)

# Create the plot
p <- ggplot(df_joined, aes(x = df_long.hour, y = df_long.Value, fill = df_long.Measurement)) + 
  geom_bar(stat = "identity", position = "stack", width = 2) +
  geom_errorbar(aes(ymin = ifelse(df_long.Measurement=='Starch', # since GFS is stacked on starch
                                  df_long.Value, df_long.carbs_total) - df_std.Std, 
                    ymax = ifelse(df_long.Measurement=='Starch', # since GFS is stacked on starch
                                  df_long.Value, df_long.carbs_total) + df_std.Std)) +  # To align the error bars correctly
  facet_wrap(~ df_long.DOY_date) +         # Creates a subplot for each date
  theme_minimal() +
  labs(fill ='Type', x = "Hour", y = "Leaf Carbohydrate Concentration (nmol/mg)",
       labels=c('GFS total', 'Starch')) +
  scale_fill_brewer(palette = "Set2") + # Choose the color palette
  scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 24))+
  scale_y_continuous(breaks = seq(0, 1500, 300), limits = c(0, 1500))+
  theme(legend.position = 'right',
        legend.text = element_text(size = 8),  # Adjust text size
        legend.key.size = unit(0.4, "cm"))

# Display the plot
print(p)


# Compare Land-based TNC
TNC.data$Source <- 'Measured'
TNC.data$time <- TNC.data$DOY + TNC.data$hour/24

data_long <- tidyr::gather(TNC.data[,c('time', 'Leaf', 'Stem')], key="Type", value="Value", -time)

# Create line plot
ggplot(data_long, aes(x=time, y=Value, color=Type)) +
  theme_classic() +
  scale_color_discrete(labels = c("Leaf", "Stem"))+
  theme(legend.position = c(0.8, 0.8), 
        panel.background = element_rect(fill = "transparent",colour = NA))+
  geom_line()+
  labs(x='Day of Year (2022)',
       y='Substrate C (mol C / Mg Dry Mass)')

# Plot measured vs simulated TNC per m2 
# Measured alone
for (i in 1:dim(TNC.data)[1]){
  TNC.data[i, 'Leaf_substrate_carbon'] <- 
    TNC.data[i, 'Leaf'] * 10^(-4) *
    soybean_ld11_biomass_2022[which(soybean_ld11_biomass_2022$doy==TNC.data[i, 'Biomass_DOY']), 'leaf']
  
  TNC.data[i, 'Stem_substrate_carbon'] <- 
    TNC.data[i, 'Stem'] * 10^(-4) *
    soybean_ld11_biomass_2022[which(soybean_ld11_biomass_2022$doy==TNC.data[i, 'Biomass_DOY']), 'stem']
}

data_long_per_m2 <- tidyr::gather(TNC.data[,c('time', 'Leaf_substrate_carbon', 'Stem_substrate_carbon')], key="Type", value="Value", -time)

# Create line plot
ggplot(data_long_per_m2, aes(x=time, y=Value, color=Type)) +
  theme_classic() +
  scale_color_discrete(labels = c("Leaf", "Stem"))+
  theme(legend.position = "bottom", 
        panel.background = element_rect(fill = "transparent",colour = NA))+
  geom_line()+
  labs(title = 'Substrate C per land area', 
       x='Day of Year (2022)',
       y='Substrate C (mol C/m^2)')

# Simulated + Measured
# Leaf
# inds <- which(biocro_result$time %in% TNC.data$time)
# TNC.data$Leaf_substrate_carbon <- biocro_result$Leaf[inds] * TNC.data$Leaf * 10^(-4)
# TNC.data$Stem_substrate_carbon <- biocro_result$Stem[inds] * TNC.data$Stem * 10^(-4)

biocro_result$Source <- 'Simulated'
Leaf.carb.data <- rbind(biocro_result[,c('time','Leaf_substrate_carbon','Source')], 
                        TNC.data[, c('time','Leaf_substrate_carbon','Source')])

ggplot(Leaf.carb.data, aes(time, Leaf_substrate_carbon, group = Source)) + 
  geom_point(aes(shape=Source, color=Source, size=Source))+
  scale_shape_manual(values=c(18, 16)) +
  scale_size_manual(values=c(2.5, 0.5)) +
  scale_color_manual(values = c(col.palette.muted[1],"grey70"))+
  theme_classic() +
  theme(legend.position = c(0.85, 0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.6, 0.1)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(), 
       x='Day of Year (2022)',
       y='Leaf Substrate C (mol C/m^2)')

# Stem
biocro_result$Source <- 'Simulated'
Stem.carb.data <- rbind(biocro_result[,c('time','Stem_substrate_carbon','Source')], 
                        TNC.data[, c('time','Stem_substrate_carbon','Source')])

ggplot(Stem.carb.data, aes(time, Stem_substrate_carbon, group = Source)) + 
  geom_point(aes(shape=Source, color=Source, size=Source))+
  scale_shape_manual(values=c(18, 16)) +
  scale_size_manual(values=c(2, 0.5)) +
  scale_color_manual(values = c(col.palette.muted[2],"grey70"))+
  theme_classic() +
  theme(legend.position = c(0.85, 0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.6, 0.1)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(), 
       x='Day of Year (2022)',
       y='Stem Substrate C (mol C/m^2)')

# Plot measured vs simulated TNC per mass 
# Simulated TNC alone
sim_substrate_C_by_mass <- data.frame(
  time = biocro_result$time,
  hour = biocro_result$hour,
  Leaf = 10^4 * biocro_result$Leaf_substrate_carbon / biocro_result$Leaf,
  Stem = 10^4 * biocro_result$Stem_substrate_carbon / biocro_result$Stem
)

data_long_per_mass <- tidyr::gather(sim_substrate_C_by_mass[1:2600,c('time', 'Leaf', 
                                                                     'Stem')], 
                                    key="Type", value="Value", -time)

# Create line plot
ggplot(data_long_per_mass, aes(x=time, y=Value/1000, color=Type)) +
  theme_classic() +
  scale_color_discrete(labels = c("Leaf", "Stem"))+
  theme(legend.position = "bottom", 
        panel.background = element_rect(fill = "transparent",colour = NA))+
  geom_line()+
  labs(title = 'Substrate C per mass', 
       x='Day of Year (2022)',
       y='Substrate C (mol C/ kg)')

# Simulated + Measured
# Leaf
sim_leaf_tnc_by_mass <- sim_substrate_C_by_mass[, c('time','Leaf','hour')]
colnames(sim_leaf_tnc_by_mass)[which(colnames(sim_leaf_tnc_by_mass)=='Leaf')] <- 'TNC'
sim_leaf_tnc_by_mass$Source <- 'Simulated'

leaf.tnc.mean <- leaf.tnc.all[,c('time','Carb_total_mean', 'hour','DOY_date')]
colnames(leaf.tnc.mean)[which(colnames(leaf.tnc.mean)=='Carb_total_mean')] <- 'TNC'
leaf.tnc.mean$Source <- 'Measured (mean)'

leaf.tnc.long$Source <- 'Measured (individual)'

Leaf.carb.data <- rbind(sim_leaf_tnc_by_mass[,c('time','TNC','Source')], 
                        leaf.tnc.mean[, c('time','TNC','Source')],
                        leaf.tnc.long[, c('time', 'TNC', 'Source')] )

ggplot(Leaf.carb.data, aes(time, TNC, group = Source)) + 
  geom_point(aes(shape=Source, color=Source, size=Source))+
  scale_shape_manual(values=c(3, 18, 16)) +
  scale_size_manual(values=c(1, 4, 0.5)) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  #scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 1)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(), 
       x='Day of Year (2022)',
       y='Leaf Substrate C (mol C/Mg)')

# Stem
sim_stem_tnc_by_mass <- sim_substrate_C_by_mass[, c('time','Stem','hour')]
colnames(sim_stem_tnc_by_mass)[which(colnames(sim_stem_tnc_by_mass)=='Stem')] <- 'TNC'
sim_stem_tnc_by_mass$Source <- 'Simulated'

stem.tnc.mean <- stem.tnc.all[,c('time','Carb_total_mean', 'hour','DOY_date')]
colnames(stem.tnc.mean)[which(colnames(stem.tnc.mean)=='Carb_total_mean')] <- 'TNC'
stem.tnc.mean$Source <- 'Measured (mean)'

stem.tnc.long$Source <- 'Measured (individual)'

Stem.carb.data <- rbind(sim_stem_tnc_by_mass[,c('time','TNC','Source')], 
                        stem.tnc.mean[, c('time','TNC','Source')],
                        stem.tnc.long[, c('time', 'TNC', 'Source')] )

ggplot(Stem.carb.data, aes(time, TNC, group = Source)) + 
  geom_point(aes(shape=Source, color=Source, size=Source))+
  scale_shape_manual(values=c(3, 18, 16)) +
  scale_size_manual(values=c(1, 4, 0.5)) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(), 
       x='Day of Year (2022)',
       y='Stem Substrate C (mol C/Mg)')

# diurnal changes of substrate C
# take out the last TNC data because in simulation the crop has stopped growing
sim_leaf_tnc_by_mass$DOY <- as.integer(sim_leaf_tnc_by_mass$time)
leaf.tnc.mean$DOY <- as.integer(leaf.tnc.mean$time)
ind.sampling.days <- which(sim_leaf_tnc_by_mass$DOY %in% leaf.tnc.all$DOY)
sim_leaf_tnc_sampling_days <- sim_leaf_tnc_by_mass[ind.sampling.days,]

date.list <- data.frame(
  DOY = c(186, 187, 209, 236, 258),
  DOY_date = c('DOY 186: 07/05/22', 'DOY 187: 07/06/22', 'DOY 209: 07/28/22', 'DOY 236: 08/24/22', 'DOY 258: 09/15/22'))

sim_leaf_tnc_sampling_days <- left_join(sim_leaf_tnc_sampling_days, date.list, by = 'DOY')

col.names <- c('DOY_date','hour','TNC', 'Source')
leaf.tnc.sampling.days <- rbind(sim_leaf_tnc_sampling_days[,col.names],
                                leaf.tnc.long[,col.names],
                                leaf.tnc.mean[,col.names])

ggplot(leaf.tnc.sampling.days, aes(hour, TNC, group = Source)) + 
  geom_point(aes(shape=Source, color=Source, size=Source))+
  facet_wrap("DOY_date") +
  scale_shape_manual(values=c(3, 18, 16)) +
  scale_size_manual(values=c(1, 2.5, 0.5)) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  scale_x_continuous(breaks = seq(0,24,6))+
  labs(title=element_blank(), 
       x='Hour of the Day',
       y='Leaf Substrate C (mol C/Mg)')

sim_stem_tnc_by_mass$DOY <- as.integer(sim_stem_tnc_by_mass$time)
stem.tnc.mean$DOY <- as.integer(stem.tnc.mean$time)
ind.sampling.days <- which(sim_stem_tnc_by_mass$DOY %in% stem.tnc.all$DOY)
sim_stem_tnc_sampling_days <- sim_stem_tnc_by_mass[ind.sampling.days,]

date.list <- data.frame(
  DOY = c(186, 187, 209, 236, 258, 278),
  DOY_date = c('DOY 186: 07/05/22', 'DOY 187: 07/06/22', 
               'DOY 209: 07/28/22', 'DOY 236: 08/24/22', 
               'DOY 258: 09/15/22', 'DOY 278: 10/05/22'))

sim_stem_tnc_sampling_days <- left_join(sim_stem_tnc_sampling_days, date.list, by = 'DOY')

col.names <- c('DOY_date','hour','TNC', 'Source')
stem.tnc.sampling.days <- rbind(sim_stem_tnc_sampling_days[,col.names],
                                stem.tnc.long[,col.names],
                                stem.tnc.mean[,col.names])

ggplot(stem.tnc.sampling.days, aes(hour, TNC, group = Source)) + 
  geom_point(aes(shape=Source, color=Source, size=Source))+
  facet_wrap("DOY_date") +
  scale_shape_manual(values=c(3, 18, 16)) +
  scale_size_manual(values=c(1, 2.5, 0.5)) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  scale_x_continuous(breaks = seq(0,24,6))+
  labs(title=element_blank(), 
       x='Hour',
       y='Stem Substrate C (mol C/Mg)')

leaf.tnc.sampling.days$Organ <- 'Leaf'
stem.tnc.sampling.days$Organ <- 'Stem'

TNC.sampling.days <- rbind(leaf.tnc.sampling.days, stem.tnc.sampling.days)
TNC.sampling.days <- TNC.sampling.days[-which(TNC.sampling.days$DOY_date=='DOY 278: 10/05/22'),]
ggplot(TNC.sampling.days, aes(hour, TNC, group = Source)) + 
  geom_point(aes(shape=Source, color=Organ, size=Source))+
  facet_wrap("DOY_date") +
  scale_shape_manual(values=c(3, 18, 16)) +
  scale_size_manual(values=c(1, 2.5, 0.5)) +
  scale_color_manual(values = col.palette.muted)+
  theme_classic() +
  theme(legend.position = c(0.84, 0.18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  scale_x_continuous(breaks = seq(0,24,6))+
  labs(title=element_blank(), 
       x='Hour',
       y='Substrate C (mol C/Mg)')

source('plot_partitioning.R')
allocation_percentage_tall <- plot_partitioning(biocro_result, '2022')
allocation_percentage_tall$Sunlight <- 'Measured'
allocation_percentage_lowlight_tall <- plot_partitioning(biocro_result_lowlight, '2022')
allocation_percentage_lowlight_tall$Sunlight <- 'Measured -20%'
allocation_percentage_tall_combined <- rbind(allocation_percentage_tall, allocation_percentage_lowlight_tall)

ggplot() + theme_classic() +
  geom_point(data = allocation_percentage_tall_combined, aes(x = DOY, y = Percentage, color = Sunlight), size = 1) +
  facet_wrap("Organ") +
  theme(plot.title=element_text(size=size.title, hjust=0.5),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_y_continuous(limits = c(-20, 140), breaks = seq(-20, 140, 20)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(), x='Day of Year (2022)',y='Allocation %')

# copied
result <- biocro_result
year <- '2022'
canopy_assim_daily <- aggregate(result$canopy_assimilation_rate,list(result$doy), 
                                FUN=sum) * 0.6 / 180.156e-3
leaf_reuse_daily <- aggregate(result$Leaf_senescence_rate * 
                                parameters$Leaf_senescence_reuse_factor,
                              list(result$doy), FUN=sum)
stem_reuse_daily <- aggregate(result$Stem_senescence_rate * 
                                parameters$Stem_senescence_reuse_factor,
                              list(result$doy), FUN=sum)
root_reuse_daily <- aggregate(result$Root_senescence_rate * 
                                parameters$Root_senescence_reuse_factor,
                              list(result$doy), FUN=sum)
net_subC_input <- data.frame(Group.1 = canopy_assim_daily$Group.1,
                             x = canopy_assim_daily$x + 
                               leaf_reuse_daily$x + 
                               stem_reuse_daily$x + 
                               root_reuse_daily$x)
leaf_export_daily <- aggregate(result$substrate_transport_Leaf_to_Stem,
                               list(result$doy), FUN=sum)
pod_allocation_daily <- aggregate(result$substrate_transport_Stem_to_Pod,
                                  list(result$doy), FUN=sum) 
root_allocation_daily <- aggregate(result$substrate_transport_Stem_to_Root,
                                   list(result$doy), FUN=sum) + root_reuse_daily
stem_allocation_daily <- leaf_export_daily + stem_reuse_daily - 
  pod_allocation_daily - root_allocation_daily
leaf_allocation_daily <- canopy_assim_daily + leaf_reuse_daily - leaf_export_daily

allocation_percentage <- data.frame(DOY = leaf_export_daily$Group.1, 
                                    Leaf = 100 * leaf_allocation_daily$x / net_subC_input$x,
                                    Stem = 100 * stem_allocation_daily$x / net_subC_input$x,
                                    Root = 100 * root_allocation_daily$x / net_subC_input$x,
                                    Pod = 100 * pod_allocation_daily$x / net_subC_input$x)
# To check if the sum is 100%
# allocation_percentage$tot_percentage <- allocation_percentage$Leaf+
#                                         allocation_percentage$Stem+
#                                         allocation_percentage$Root+
#                                         allocation_percentage$Pod
