# Generate all the plots in the Piao et al. 2026 "Integrating carbon utilization 
# and transport processes into a crop growth model enables the prediction of 
# emergent soybean carbon allocation behavior" manuscript.


library(BioCro)           # partial_run_biocro(), soybean model/parameter set
library(UTRSoybeanBML)    # UTRSoybeanBML: the utilization-transport-resistance 
                          # (UTR) C allocation BioCro Module Library (BML)
library(BioCroWater)      # BioCroWater: soil water BML
library(ggplot2)
library(grid)             # textGrob, gpar, unit, grid.newpage/draw/text
library(gridExtra)        # grid.arrange, arrangeGrob
library(dplyr)
library(tidyr)            # pivot_longer, gather
library(purrr)            # map_dfr, used in plot_allocation_by_use.R

# Clear workspace
rm(list=ls())

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Parameters
source('../ParameterOptimization/soybean_parameter_expansion.R')
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')

# Source mse calculation and partitioning/biomass plotting functions
source('CalculateScripts/calculate_mse.R')
source('PlotScripts/plot_partitioning.R')
source('PlotScripts/plot_biomass.R')

# All figures produced by this script (and by the scripts it sources) are
# written here; created fresh on every run.
FIGURE_DIR <- 'GeneratedFigures'
dir.create(FIGURE_DIR, showWarnings = FALSE)

# Toggle the RUN_SENSITIVITY / PLOT_SOURCE_SINK_MANIPULATION / PLOT_LAI flags
# below to control which optional figures to produce or  
# if sensitivity analysis is performed, which can take longer to run.
RUN_SENSITIVITY <- TRUE                # pod mass sensitivity to +/-10% inputs; slow
PLOT_SOURCE_SINK_MANIPULATION <- TRUE  # shade/pod-removal/defoliation/hail figures; slow
PLOT_LAI <- FALSE                      # SoyFACE and LD11 simulated-vs-observed LAI figures

if(PLOT_LAI){
  source('PlotScripts/plot_lai_comparison.R')
}

# years, sowing dates, and harvesting dates of growing seasons being fit to
years <- c('2002', '2004', '2005', '2006')
sow.date <- c(152, 149, 148, 148)
harv.date <- c(288, 289, 270, 270)

# Model set up
initial_values <- soybean$initial_values
direct_modules <- soybean$direct_modules
differential_modules <- soybean$differential_modules
solver <- soybean$ode_solver
parameters <- soybean$parameters

# Set up BioCroWater
source('../Data/Soybean-BioCro_Parameters/set_up_BioCroWater.R')
initial_values       <- set_init_values(initial_values)
parameters           <- set_parameters(parameters)
parameters$kd        <- parameters$k_diffuse
direct_modules       <- set_direct_modules(direct_modules)
differential_modules <- set_differential_modules(differential_modules)

# Update UTR modules
source('../Data/Soybean-BioCro_Parameters/set_up_UTRSoybeanBML.R')
initial_values       <- set_init_values(initial_values)
parameters           <- set_parameters(parameters)
direct_modules       <- set_direct_modules(direct_modules)
differential_modules <- set_differential_modules(differential_modules)

# Update UTR parameters
fitted.utr.params <- optim_params_conversion(optim_params_short_SoyFACE)
names(fitted.utr.params) <- arg_names
parameters <-c(parameters, fitted.utr.params)[!duplicated(c(names(parameters),
                                                            names(fitted.utr.params)),
                                                          fromLast = TRUE)]

# Minor adjustments
parameters$time_zone_offset <- -6
initial_values$DVI <- -1

# Initialize lists
results <- list()
results.elevCO2 <- list()
weather.growingseason <- list()
soybean_optsolver <- list()
ExpBiomass <- list()
ExpBiomass.elevCO2 <- list()
LAI <- list()
LAI.elevCO2 <- list()
ExpBiomass.std <- list()
ExpBiomass.elevCO2.std <- list()
figs <- list()
figs.elevCO2 <- list()
lai.figs <- list()
allocation.figs <- list()
allocation.elevCO2.figs <- list()

pod_sensitivity_lower <- NULL  # for -10%
pod_sensitivity_upper <- NULL  # for +10%

# =============================================================================
# Pioneer 93B15 @ ambient CO2: 2002, 2004-2006 growing seasons
# =============================================================================
co2_opt = '_ambient_'
parameters$Catm      <- 372

for (i in 1:length(years)) {
  yr <- years[i]
  weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))
  sd.idx <- which(weather$doy == sow.date[i])[12]
  hd.ind <- which(weather$doy == harv.date[i])[24]

  weather.growingseason[[i]] <- weather[sd.idx: hd.ind,]

  ExpBiomass[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, '_ambient_biomass.csv'))
  colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

  ExpBiomass.std[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, '_ambient_biomass_std.csv'))
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

  if(PLOT_LAI && i>1){
    LAI[[i]]<-read.csv(file=paste0('../Data/lai/',yr,'_ambient_lai.csv'))
    LAI.elevCO2[[i]]<-read.csv(file=paste0('../Data/SoyFACE_data/lai',yr,'_elevated_lai.csv'))
  }

  soybean_optsolver[[i]] <- partial_run_biocro(initial_values,
                                               parameters,
                                               weather.growingseason[[i]],
                                               direct_modules,
                                               differential_modules,
                                               solver,
                                               arg_names,
                                               verbose = FALSE)

  result <- soybean_optsolver[[i]](
    optim_params_conversion(optim_params_short_SoyFACE))
  results[[i]] <- result

  # figure out the DVI when shell stops growing
  final_shell <- tail(ExpBiomass[[i]]$Pod,1) - tail(ExpBiomass[[i]]$Seed, 1)
  print(paste0('Simulated pod mass has reached the measured final shell mass at DVI ',
               result$DVI[which.min(abs(result$Pod-final_shell))]))

  # plot simulation vs observation biomass and calculate mse (mse is printed)
  figs[[i]] <- plot_SoyFACE_biomass(
    result, ExpBiomass[[i]], ExpBiomass.std[[i]], co2_opt, years[i])

  # plot partitioning
  allocation.figs[[i]] <- plot_partitioning(result, years[i])

  # calculate sensitivity
  if (RUN_SENSITIVITY) {source('CalculateScripts/calculate_sensitivity.R')}

}
# =============================================================================
# Pioneer 93B15 @ elevated CO2: 2002, 2004-2006 growing seasons
# =============================================================================
co2_opt = '_co2_'
parameters$Catm <- 550

# Initialize lists
soybean_optsolver <- list()

for (i in 1:length(years)) {
  yr <- years[i]

  ExpBiomass.elevCO2[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass.csv'))
  colnames(ExpBiomass.elevCO2[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

  ExpBiomass.elevCO2.std[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
  colnames(ExpBiomass.elevCO2.std[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

  soybean_optsolver[[i]] <- partial_run_biocro(initial_values,
                                               parameters,
                                               weather.growingseason[[i]],
                                               direct_modules,
                                               differential_modules,
                                               solver,
                                               arg_names,
                                               verbose = FALSE)

  result <- soybean_optsolver[[i]](optim_params_conversion(optim_params_short_SoyFACE))
  results.elevCO2[[i]] <- result

  # plot simulation vs observation biomass and calculate mse (mse is printed)
  figs.elevCO2[[i]] <- plot_SoyFACE_biomass(result, ExpBiomass.elevCO2[[i]], ExpBiomass.elevCO2.std[[i]], co2_opt, years[i])

  # plot partitioning
  allocation.elevCO2.figs[[i]] <- plot_partitioning(result, years[i])

  # calculate sensitivity
  if (RUN_SENSITIVITY) {source('CalculateScripts/calculate_sensitivity.R')}

}



# extract legend from a ggplot object, to be shared across combined panels
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# =============================================================================
# Figure S: SoyFACE simulated vs observed LAI
# =============================================================================
if(PLOT_LAI){
  for (i in 2:4){
    lai.figs[[i-1]] <- plot_amb_elev_lai(
      results[[i]], results.elevCO2[[i]],
      years[[i]], LAI[[i]], LAI.elevCO2[[i]])
  }

  combined_graph.lai <-
    grid.arrange(
      arrangeGrob(
        textGrob(bquote("LAI"~(m^2~"/"~m^2)), rot = 90, gp=gpar(fontsize=12))),
        arrangeGrob(
          arrangeGrob(
            lai.figs[[1]] + theme(legend.position="none"),
            lai.figs[[2]] + theme(legend.position="none"),
            lai.figs[[3]] + theme(legend.position="none"),
            ncol = 3), ncol = 1),
      ncol=2, widths=c(0.3, 5))

  ggsave(file.path(FIGURE_DIR, 'FigS-soyface-lai.png'),
         plot = combined_graph.lai,
         width = 8,
         height = 2,
         units = "in",
         dpi = 600
  )
}

# =============================================================================
# Figure 5 Root-Shoot Ratios
# =============================================================================
# Build combined data frame across all years
root_shoot_all <- do.call(rbind, lapply(1:length(years), function(i) {
  data.frame(
    year           = years[i],
    fractional_doy = results[[i]]$fractional_doy,
    DVI            = results[[i]]$DVI,
    ambient_ratio  = with(results[[i]],      Root / (Leaf + Stem + Pod)),
    elevated_ratio = with(results.elevCO2[[i]], Root / (Leaf + Stem + Pod))
  )
}))

# Pivot to long format for ggplot
root_shoot_long <- root_shoot_all %>%
  pivot_longer(
    cols      = c(ambient_ratio, elevated_ratio),
    names_to  = "treatment",
    values_to = "ratio"
  ) %>%
  mutate(treatment = recode(treatment,
                            ambient_ratio  = "Ambient CO2",
                            elevated_ratio = "Elevated CO2"
  ))

# Plot
root_shoot_ratio_plot <- ggplot(root_shoot_long, 
                                aes(x = fractional_doy, 
                                    y = ratio, 
                                    colour = treatment)) +
  geom_line() +
  facet_wrap(~ year, ncol = 2) +
  scale_colour_manual(values = c("Ambient CO2" = "steelblue", 
                                 "Elevated CO2" = "tomato")) +
  labs(
    x      = "Day of Year",
    y      = "Root:Shoot Ratio",
    colour = "CO2 Level"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

print(root_shoot_ratio_plot)
ggsave(file.path(FIGURE_DIR, 'Fig-root-shoot-ratio.png'),
       plot = root_shoot_ratio_plot,
       width = 6,
       height = 4,
       units = "in",
       dpi = 600
)


root_shoot_mean <- root_shoot_long %>%
  mutate(dvi_bin = round(DVI / 0.01) * 0.01) %>%   # bin on DVI instead of fractional_doy
  group_by(dvi_bin, treatment) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    min_ratio  = min(ratio,  na.rm = TRUE),
    max_ratio  = max(ratio,  na.rm = TRUE),
    .groups = "drop"
  )

root_shoot_mean_plot <- ggplot(root_shoot_mean,
                               aes(x = dvi_bin, colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = min_ratio, ymax = max_ratio), alpha = 0.2, colour = NA) +
  geom_line(aes(y = mean_ratio), linewidth = 0.8) +
  scale_colour_manual(values = c("Ambient CO2"  = "steelblue",
                                 "Elevated CO2" = "tomato"),
                      name = "CO2 Level") +
  scale_fill_manual(values = c("Ambient CO2"  = "steelblue",
                               "Elevated CO2" = "tomato"),
                    name = "CO2 Level") +
  labs(
    x = "DVI",
    y = "Root:Shoot Ratio"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm")
  )

print(root_shoot_mean_plot)

ggsave(file.path(FIGURE_DIR, 'Fig-root-shoot-ratio-mean.png'),
       plot   = root_shoot_mean_plot,
       width  = 4,
       height = 3,
       units  = "in",
       dpi    = 600)

# =============================================================================
# Figure S10 Ratio of Stem to Root Transport to Daily Assimilation
# =============================================================================
# Build combined data frame across all years
transportation_all <- do.call(rbind, lapply(1:4, function(i) {
  
  # Ambient: daily sums then ratio
  ambient_daily <- results[[i]] %>%
    group_by(doy) %>%
    summarise(
      Stem_to_Root  = sum(substrate_transport_Stem_to_Root),
      total_assimilation  = sum(canopy_assimilation_rate),
      .groups = "drop"
    ) %>%
    mutate(ratio = Stem_to_Root / total_assimilation,
           treatment = "Ambient CO2")
  
  # Elevated: daily sums then ratio
  elevated_daily <- results.elevCO2[[i]] %>%
    group_by(doy) %>%
    summarise(
      Stem_to_Root  = sum(substrate_transport_Stem_to_Root),
      total_assimilation  = sum(canopy_assimilation_rate),
      .groups = "drop"
    ) %>%
    mutate(ratio = Stem_to_Root / total_assimilation,
           treatment = "Elevated CO2")
  
  rbind(ambient_daily, elevated_daily) %>%
    mutate(year = years[i])
}))

transport_assimilation_ratio_plot <- 
  ggplot(transportation_all, aes(x = doy, y = ratio, colour = treatment)) +
  geom_line() +
  facet_wrap(~ year, ncol = 2) +
  scale_colour_manual(values = c("Ambient CO2" = "steelblue", 
                                 "Elevated CO2" = "tomato")) +
  labs(
    x      = "Day of Year",
    y      = "Daily Stem to Root Transport / Daily Canopy Assimilation",
    colour = "Treatment"
  ) +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  )
print(transport_assimilation_ratio_plot)
ggsave(file.path(FIGURE_DIR, 'FigS-transport-assimilation-ratio.png'),
       plot = transport_assimilation_ratio_plot,
       width = 8,
       height = 5,
       units = "in",
       dpi = 600
)

# =============================================================================
# Figure S11 Utilization Saturation Ratio 
# =============================================================================
# Build combined data frame across all years
leaf_saturation_all <- do.call(rbind, lapply(1:4, function(i) {

  data.frame(
    fractional_doy = results[[i]]$fractional_doy,
    ambient        = with(results[[i]],
                          Leaf_utilization_rate / Leaf_structural_carbon / 
                            parameters$Leaf_utilization_rate_constant),
    elevated       = with(results.elevCO2[[i]],
                          Leaf_utilization_rate / Leaf_structural_carbon / 
                            parameters$Leaf_utilization_rate_constant),
    year           = years[i]
  )
})) %>%
  pivot_longer(
    cols      = c(ambient, elevated),
    names_to  = "treatment",
    values_to = "saturation"
  ) %>%
  mutate(treatment = recode(treatment,
                            ambient  = "Ambient CO2",
                            elevated = "Elevated CO2"
  ))

# Plot
utilization_saturation_plot <- 
  ggplot(leaf_saturation_all, aes(x = fractional_doy, y = saturation,
                                  colour = treatment, 
                                  linewidth = treatment, 
                                  alpha = treatment)) +
  geom_line() +
  facet_wrap(~ year, ncol = 2) +
  scale_colour_manual(values = c("Ambient CO2" = "steelblue", 
                                 "Elevated CO2" = "tomato")) +
  scale_linewidth_manual(values = c("Ambient CO2" = 0.8,         
                                    "Elevated CO2" = 0.4)) +
  scale_alpha_manual(values = c("Ambient CO2" = 1,           
                                "Elevated CO2" = 0.5)) +
  labs(x = "Day of Year", y = "Leaf Utilization Rate Saturation",
       colour = "Treatment", linewidth = "Treatment", alpha = "Treatment") +
  theme_bw() +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold"))

print(utilization_saturation_plot)
ggsave(file.path(FIGURE_DIR, 'FigS-utilization_saturation.png'),
       plot = utilization_saturation_plot,
       width = 8,
       height = 6,
       units = "in",
       dpi = 600
)

# =============================================================================
# Source Sink Manipulations
# =============================================================================

if(PLOT_SOURCE_SINK_MANIPULATION){
  # Figure 8 Shade treatment
  source('PlotScripts/plot_shade.R')
  # Figure 9 Sod removal treatment
  source('PlotScripts/plot_pod_removal.R')
  # Figure 10 Defoliation treatment
  source('PlotScripts/plot_defoliation.R')
  # Figure 11 Hail event
  source('PlotScripts/plot_2003_hail.R')
}

# =============================================================================
# Plot LD11 Biomass and substrate C concentrations
# =============================================================================
years <- c('2021', '2022', '2023', '2024')
Catms <- c(414.7, 417.2, 419.3, 422.8) # from NOAA

# Update LD11 photosynthetic parameters with 2022 LICOR measurement results
# Analyzed with PhotoGEA (Lochocki et al, )
updated_parameters <- read.csv('../Data/Soybean-BioCro_Parameters/updated-ld11-parameters-2022.csv')
updated_idx <- match(updated_parameters$New.variable.name, names(parameters))
parameters[updated_idx] <- updated_parameters$LD11

ExpBiomass <- list()
ld11.results <- list()
ld11.figs <- list()
allocation.ld11.figs <- list()

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

for (i in 1:length(years)){
  yr <- years[i]
  ExpBiomass[[i]] <- loadRData(paste0('../Data/EnergyFarm_data/soybean_ld11_biomass_', years[i], '.RData'))
  sow.time <- min(ExpBiomass[[i]]$time)
  harv.time <- max(ExpBiomass[[i]]$time)
  weather <- loadRData(paste0('../Data/Weather_data/weather', years[i], '_hourly.RData'))
  parameters$Catm <- Catms[i]
  # update initial values
  sub_frac <- 0.1           # substrate_fraction
  str_frac <- 1 - sub_frac  # structural_fraction
  seed_mass <- ExpBiomass[[i]]$initial_seed[1]
  leaf_frac <- 0.45
  stem_frac <- 0.25
  root_frac <- 0.3

  cf <- 0.3
  updated_utr_initial_values <- list(
    Leaf_substrate_carbon = sub_frac * seed_mass * leaf_frac / cf,
    Leaf_structural_carbon = str_frac * seed_mass * leaf_frac / cf,
    Stem_substrate_carbon = sub_frac * seed_mass * stem_frac / cf,
    Stem_structural_carbon = str_frac * seed_mass * stem_frac / cf,
    Root_substrate_carbon =  sub_frac * seed_mass * root_frac / cf,
    Root_structural_carbon = str_frac * seed_mass * root_frac / cf)

  initial_values[names(updated_utr_initial_values)] <- updated_utr_initial_values
  # weather file
  weather.aftersowing <- weather[(((weather$doy-1)*24 + weather$hour) >= sow.time) & ((weather$doy-1)*24 + weather$hour) <= harv.time, ]
  weather.aftersowing$time_zone_offset <- NULL
  parameters$RL_at_25 <- 1.28

  soybean_optsolver[[i]] <- partial_run_biocro(initial_values,
                                            parameters,
                                            weather.aftersowing,
                                            direct_modules,
                                            differential_modules,
                                            solver,
                                            arg_names,
                                            verbose = FALSE)

  result <- soybean_optsolver[[i]](optim_params_conversion(optim_params_short_SoyFACE))

  ld11.results[[i]] <- result

  # plot LD11 biomass
  ld11.figs[[i]] <- plot_ld11_biomass(result, ExpBiomass[[i]])

  # plot partitioning
  allocation.ld11.figs[[i]] <- plot_partitioning(result, years[i])

  # calculate sensitivity
  if (RUN_SENSITIVITY) {source('CalculateScripts/calculate_sensitivity.R')}
}

if (RUN_SENSITIVITY) {
  print(pod_sensitivity_lower)
  print(pod_sensitivity_upper)
}

# =============================================================================
# Plot LAI
# =============================================================================
if(PLOT_LAI){
  # Combine all years into one data frame
  lai_df <- do.call(rbind, lapply(1:4, function(i) {
    # Simulated LAI (line)
    sim <- data.frame(
      time  = ld11.results[[i]]$time,
      value = ld11.results[[i]]$lai,
      type  = "Simulated",
      year  = years[i]
    )

    # Observed LAI (dots) - reshape 3 columns to long format
    obs_wide <- data.frame(
      time                          = ExpBiomass[[i]]$time,
      LAI_from_LMA                  = ExpBiomass[[i]]$LAI_from_LMA,
      LAI_from_planting_density     = ExpBiomass[[i]]$LAI_from_planting_density,
      LAI_from_measured_population  = ExpBiomass[[i]]$LAI_from_measured_population
    )

    obs_long <- pivot_longer(obs_wide,
                             cols      = -time,
                             names_to  = "type",
                             values_to = "value")
    obs_long$year <- years[i]

    rbind(sim, obs_long)
  }))

  # Separate simulated and observed for different geoms
  sim_df <- subset(lai_df, type == "Simulated")
  obs_df <- subset(lai_df, type != "Simulated")

  ld11_lai_comparison <- ggplot() +
    geom_line(data = sim_df,
              aes(x = (time/24 + 1), y = value),
              color = "black", linewidth = 0.8) +
    geom_point(data = obs_df,
               aes(x = (time/24 + 1), y = value, color = type, shape = type),
               na.rm = TRUE, size = 2.5) +
    facet_wrap(~ year, ncol = 2, scales = "free_x") +
    labs(
      x      = "DOY",
      y      = "LAI",
      color  = "Observed",
      shape  = "Observed",
      title  = "Simulated vs Observed LD11 LAI"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  print(ld11_lai_comparison)
  ggsave(file.path(FIGURE_DIR, 'FigS-ld11-lai.png'),
         plot = ld11_lai_comparison,
         width = 8,
         height = 6,
         units = "in",
         dpi = 600
  )
}


# =============================================================================
# Figure S1 C allocation based on net import
# =============================================================================
common_legend <- g_legend(allocation.ld11.figs[[1]])
combined_graph.allocation <- 
  grid.arrange(
    arrangeGrob(textGrob('Allocation %', rot = 90, gp=gpar(fontsize=12))),
    arrangeGrob(arrangeGrob(allocation.figs[[1]] + theme(legend.position="none"),
                            allocation.figs[[2]] + theme(legend.position="none"),
                            allocation.figs[[3]] + theme(legend.position="none"),
                            allocation.figs[[4]] + theme(legend.position="none"),
                            ncol = 4, top = 'Pioneer 93B15 at ambient CO2'),
                arrangeGrob(allocation.elevCO2.figs[[1]] + theme(legend.position="none"),
                            allocation.elevCO2.figs[[2]] + theme(legend.position="none"),
                            allocation.elevCO2.figs[[3]] + theme(legend.position="none"),
                            allocation.elevCO2.figs[[4]] + theme(legend.position="none"),
                            ncol = 4, top = 'Pioneer 93B15 at elevated CO2'),
                arrangeGrob(allocation.ld11.figs[[1]] + theme(legend.position="none"),
                            allocation.ld11.figs[[2]] + theme(legend.position="none"),
                            allocation.ld11.figs[[3]] + theme(legend.position="none"),
                            allocation.ld11.figs[[4]] + theme(legend.position="none"),
                            ncol = 4, top = 'LD11-2170 at ambient CO2'),
              nrow = 3),
    common_legend,
    ncol=3, widths=c(0.2, 5, 0.8))

ggsave(file.path(FIGURE_DIR, 'FigS-allocation-netimport.png'),
       plot = combined_graph.allocation,
       width = 10,
       height = 7,
       units = "in",
       dpi = 600
)

# =============================================================================
# Plot usage based allocation
source('PlotScripts/plot_allocation_by_use.R') # Figures 3, 4, S2, S3, S5
# =============================================================================


# =============================================================================
# Figure 2 Biomass
# =============================================================================
get_legend_grob <- function(plot) {
  g <- ggplotGrob(plot)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

common_legend_horizontal <- get_legend_grob(
  ld11.figs[[1]] +
    theme(legend.position = "bottom",
          legend.direction = "horizontal") +
    guides(color = guide_legend(nrow = 1, title.position = "left"),
           fill = guide_legend(nrow = 1, title.position = "left"))
)

all_biomass_plot <- grid.arrange(
  arrangeGrob(
    arrangeGrob(
      textGrob('Biomass (Mg / ha)', rot = 90), # Col1: y-axis Grob
      arrangeGrob( # Col 2: Biomass figures
        arrangeGrob( # Col 2, Row 1: Pioneer 93B15 ambient CO2 Grob
              figs[[1]] + theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank(),
                                legend.position="none")
                        + labs(title = '2002'),
              figs[[2]] + theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank(),
                                legend.position="none")
                        + labs(title = '2004'),
              figs[[3]] + theme(legend.position="none",
                                axis.title.x = element_blank(),
                                axis.text.x = element_blank())
                        + labs(title = '2005'),
              figs[[4]] + theme(legend.position="none",
                                axis.title.x = element_blank(),
                                axis.text.x = element_blank())
                        + labs(title = '2006'),
              ncol = 4,
              top = textGrob("Pioneer 93B15 at Ambient CO2",
                           gp=gpar(fontface="bold", fontsize=12))
        ),
        textGrob(""), # Row 2: spacer,
        arrangeGrob( # Row 3: Pioneer 93B15 ambient CO2 Grob
          figs.elevCO2[[1]] + theme(axis.title.x = element_blank(),
                            axis.text.x = element_blank(),
                            legend.position="none")
                            + labs(title = '2002'),
          figs.elevCO2[[2]] + theme(axis.title.x = element_blank(),
                            axis.text.x = element_blank(),
                            legend.position="none")
                            + labs(title = '2004'),
          figs.elevCO2[[3]] + theme(legend.position="none",
                            axis.title.x = element_blank(),
                            axis.text.x = element_blank())
                            + labs(title = '2005'),
          figs.elevCO2[[4]] + theme(legend.position="none",
                            axis.title.x = element_blank(),
                            axis.text.x = element_blank())
                            + labs(title = '2006'),
          ncol = 4,
          top = textGrob("Pioneer 93B15 at Elevated CO2",
                       gp=gpar(fontface="bold", fontsize=12))
        ),
        textGrob(""), # Row 4: spacer,
        arrangeGrob( # Row 5: LD11 figures
          ld11.figs[[1]] + theme(legend.position="none") + labs(title = '2021'),
          ld11.figs[[2]] + theme(legend.position="none") + labs(title = '2022'),
          ld11.figs[[3]] + theme(legend.position="none") + labs(title = '2023'),
          ld11.figs[[4]] + theme(legend.position="none") + labs(title = '2024'),
          ncol = 4,
          top = textGrob("LD11-2170 at Ambient CO2",
                         gp=gpar(fontface="bold", fontsize=12))
        ),
        nrow = 5, heights = c(1, 0.05, 1, 0.05, 1.2)
      ),
      ncol = 2, widths = c(0.3, 10)
    ),
    common_legend_horizontal,
    nrow = 2, heights = c(10, 0.5))
)

png(file.path(FIGURE_DIR, "Fig-Biomass.png"), width = 10, height = 8, units = "in", res = 600)
grid.newpage()
grid.draw(all_biomass_plot)

# Layout fractions
x_left      <- 0.3 / 10.3          # left edge of the plot column (after y-axis grob)
plot_width  <- 10  / 10.3          # width of the 4-column plot area
col_width   <- plot_width / 4

total_h     <- 10 + 0.5
plot_h_frac <- 10 / total_h        # fraction of figure height used by plots

row_heights <- c(1, 0.05, 1, 0.05, 1.2)
row_total   <- sum(row_heights)    # = 3.2

# Top y of each of the 3 plot rows (rows 1, 3, 5 in the 5-row layout)
row_tops <- c(
  1,                                           # row A-D  (row 1 top)
  1 - (row_heights[1] + row_heights[2]) / row_total,  # row E-H
  1 - (sum(row_heights[1:4])) / row_total              # row I-L
) * plot_h_frac

letters_labels <- c('A','B','C','D','E','F','G','H','I','J','K','L')
label_idx <- 1

for (row in 1:3) {
  for (col in 1:4) {
    x_pos <- x_left + (col - 1) * col_width + 0.01
    y_pos <- row_tops[row] + 0.01
    grid.text(
      paste0("(", letters_labels[label_idx], ")"),
      x    = x_pos,
      y    = y_pos,
      just = c("left", "top"),
      gp   = gpar(fontsize = 12)
    )
    label_idx <- label_idx + 1
  }
}
dev.off()

# =============================================================================
# Plot substrate C concentration
# =============================================================================
leaf.tnc.all <- read.csv('../Data/2022_Carb_data/Leaf_all.csv')
stem.tnc.all <- read.csv('../Data/2022_Carb_data/Stem_all.csv')

# convert to long format
convert_data_long <- function(tnc.all){
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
leaf.tnc.long <- convert_data_long(leaf.tnc.all)
stem.tnc.long <- convert_data_long(stem.tnc.all)

# convert unit from nmol glucose / mg (mol glucose / Mg) to mol C / kg
leaf.tnc.long$TNC <- leaf.tnc.long$TNC * 6 / 1000
stem.tnc.long$TNC <- stem.tnc.long$TNC * 6 / 1000

# Simulated
r <- ld11.results[[2]]
sim_substrate_C_by_mass <- data.frame(
  time = r$fractional_doy,
  hour = r$hour,
  Leaf = 10 * r$Leaf_substrate_carbon / r$Leaf, # convert from (mol/m2) / (Mg/ha) = 10^4 * mol / Mg = 10 mol / kg
  Stem = 10 * r$Stem_substrate_carbon / r$Stem
)
data_long_per_mass <- tidyr::gather(sim_substrate_C_by_mass[1:2600,
                                                            c('time', 'Leaf', 'Stem')],
                                    key="Type", value="Value", -time)


# Simulated + Measured
# Leaf
sim_leaf_tnc_by_mass <- sim_substrate_C_by_mass[, c('time','Leaf','hour')]
colnames(sim_leaf_tnc_by_mass)[which(colnames(sim_leaf_tnc_by_mass)=='Leaf')] <- 'TNC'
sim_leaf_tnc_by_mass$Source <- 'Simulated'

leaf.tnc.mean <- leaf.tnc.all[,c('time','Carb_total_mean', 'hour','DOY_date')]
colnames(leaf.tnc.mean)[which(colnames(leaf.tnc.mean)=='Carb_total_mean')] <- 'TNC'
leaf.tnc.mean$Source <- 'Measured (mean)'
# convert unit from nmol glucose / mg (mol glucose / Mg) to mol C / kg
leaf.tnc.mean$TNC <- leaf.tnc.mean$TNC * 6 / 1000

leaf.tnc.long$Source <- 'Measured (individual)'

Leaf.carb.data <- rbind(sim_leaf_tnc_by_mass[,c('time','TNC','Source')],
                        leaf.tnc.mean[, c('time','TNC','Source')],
                        leaf.tnc.long[, c('time', 'TNC', 'Source')] )


leaf_tnc_plot <- ggplot(Leaf.carb.data, aes(time, TNC, group = Source)) +
  geom_point(aes(shape=Source, color=Source, size=Source, alpha = Source))+
  scale_shape_manual(values=c(8, 18, 16)) +
  scale_size_manual(values=c(1.2, 4, 1)) +
  scale_alpha_manual(values=c(0.6, 1, 0.4)) +
  scale_color_manual(values=c('#D81B60','#117733','grey')) +
  theme_classic() +
  theme(plot.title=element_text(size=size.title, hjust=0.5),
        axis.text=element_text(size=size.axis),
        axis.title=element_text(size=size.axislabel),
        legend.position = c(0.85, 0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(),
       x='Day of Year (2022)',
       y = NULL)

# Stem
sim_stem_tnc_by_mass <- sim_substrate_C_by_mass[, c('time','Stem','hour')]
colnames(sim_stem_tnc_by_mass)[which(colnames(sim_stem_tnc_by_mass)=='Stem')] <- 'TNC'
sim_stem_tnc_by_mass$Source <- 'Simulated'

stem.tnc.mean <- stem.tnc.all[,c('time','Carb_total_mean', 'hour','DOY_date')]
colnames(stem.tnc.mean)[which(colnames(stem.tnc.mean)=='Carb_total_mean')] <- 'TNC'
stem.tnc.mean$Source <- 'Measured (mean)'
# convert unit from nmol glucose / mg (mol glucose / Mg) to mol C / kg
stem.tnc.mean$TNC <- stem.tnc.mean$TNC * 6 / 1000
stem.tnc.long$Source <- 'Measured (individual)'

Stem.carb.data <- rbind(sim_stem_tnc_by_mass[,c('time','TNC','Source')],
                        stem.tnc.mean[, c('time','TNC','Source')],
                        stem.tnc.long[, c('time', 'TNC', 'Source')] )

stem_tnc_plot <- ggplot(Stem.carb.data, aes(time, TNC, group = Source)) +
  geom_point(aes(shape=Source, color=Source, size=Source, alpha = Source))+
  scale_shape_manual(values=c(8, 18, 16)) +
  scale_size_manual(values=c(1, 4, 1)) +
  scale_alpha_manual(values=c(0.6, 1, 0.4)) +
  scale_color_manual(values=c('#D81B60','#999933','grey')) +
  scale_y_continuous(limits=c(0, 10), n.breaks=6)+
  theme_classic() +
  theme(plot.title=element_text(size=size.title, hjust=0.5),
        axis.text=element_text(size=size.axis),
        axis.title=element_text(size=size.axislabel),
        legend.position = c(0.85, 0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5)) +
  scale_x_continuous(breaks = seq(180,280,30))+
  labs(title=element_blank(),
       x='Day of Year (2022)',
       y = NULL)
# =============================================================================
# Figure 6 seasonal substrate C concentration
# =============================================================================
png(file.path(FIGURE_DIR, "Fig-substrate-annual.png"), width = 10, height = 3.5, units = "in", res = 600)
grid.newpage()

tnc_growingseasion_plot <- grid.arrange(
  arrangeGrob(textGrob('Substrate C (mol C / kg)', rot = 90, gp=gpar(fontsize=12))),
  arrangeGrob(
    arrangeGrob(leaf_tnc_plot, top = 'Leaf'),
    arrangeGrob(stem_tnc_plot, top = 'Stem'),
    ncol = 2
  ),
  ncol=2, widths=c(0.2, 10)
)

grid.draw(tnc_growingseasion_plot)

# x positions: left edge of each panel
# 0.2 / 10.2 ≈ 0.020 is the y-axis grob width fraction
x_left      <- 0.2 / 10.2
plot_width  <- 10  / 10.2
col_width   <- plot_width / 2

grid.text("(A)", x = x_left + 0.01,              y = 0.99, just = c("left", "top"), gp = gpar(fontsize = 12))
grid.text("(B)", x = x_left + col_width + 0.01,  y = 0.99, just = c("left", "top"), gp = gpar(fontsize = 12))

dev.off()

# Diurnal changes of substrate C
# take out the last TNC data because in simulation the crop has stopped growing
sim_leaf_tnc_by_mass$DOY <- as.integer(sim_leaf_tnc_by_mass$time)
leaf.tnc.mean$DOY <- as.integer(leaf.tnc.mean$time)
ind.sampling.days <- which(sim_leaf_tnc_by_mass$DOY %in% leaf.tnc.all$DOY)
sim_leaf_tnc_sampling_days <- sim_leaf_tnc_by_mass[ind.sampling.days,]

date.list <- data.frame(
  DOY = c(186, 187, 209, 236, 258),
  DOY_date = c('DOY 186: 07/05', 'DOY 187: 07/06', 'DOY 209: 07/28', 'DOY 236: 08/24', 'DOY 258: 09/15'))

sim_leaf_tnc_sampling_days <- left_join(sim_leaf_tnc_sampling_days, date.list, by = 'DOY')

col.names <- c('DOY_date','hour','TNC', 'Source')
leaf.tnc.sampling.days <- rbind(sim_leaf_tnc_sampling_days[,col.names],
                                leaf.tnc.long[,col.names],
                                leaf.tnc.mean[,col.names])

sim_stem_tnc_by_mass$DOY <- as.integer(sim_stem_tnc_by_mass$time)
stem.tnc.mean$DOY <- as.integer(stem.tnc.mean$time)
ind.sampling.days <- which(sim_stem_tnc_by_mass$DOY %in% stem.tnc.all$DOY)
sim_stem_tnc_sampling_days <- sim_stem_tnc_by_mass[ind.sampling.days,]

date.list <- data.frame(
  DOY = c(186, 187, 209, 236, 258, 278),
  DOY_date = c('DOY 186: 07/05', 'DOY 187: 07/06',
               'DOY 209: 07/28', 'DOY 236: 08/24',
               'DOY 258: 09/15', 'DOY 278: 10/05'))

sim_stem_tnc_sampling_days <- left_join(sim_stem_tnc_sampling_days, date.list, by = 'DOY')

col.names <- c('DOY_date','hour','TNC', 'Source')
stem.tnc.sampling.days <- rbind(sim_stem_tnc_sampling_days[,col.names],
                                stem.tnc.long[,col.names],
                                stem.tnc.mean[,col.names])

leaf.tnc.sampling.days$Organ <- 'Leaf'
stem.tnc.sampling.days$Organ <- 'Stem'

TNC.sampling.days <- rbind(leaf.tnc.sampling.days, stem.tnc.sampling.days)

TNC.sampling.days <- TNC.sampling.days[-which(TNC.sampling.days$DOY_date=='DOY 278: 10/05'),]

# =============================================================================
# Figure 7 diurnal substrate C concentration
# =============================================================================
tnc_diurnal_plot <- ggplot(TNC.sampling.days, aes(hour, TNC, group = Source)) +
  geom_point(data=subset(TNC.sampling.days, Source != 'Simulated'), aes(shape=Source, color=Organ, size=Source, alpha = Source))+
  geom_line(data=subset(TNC.sampling.days, Source == 'Simulated' & Organ == 'Leaf'), aes(color=Organ))+
  geom_line(data=subset(TNC.sampling.days, Source == 'Simulated' & Organ == 'Stem'), aes(color=Organ))+
  facet_wrap("DOY_date") +
  scale_shape_manual(values=c(8, 18)) +
  scale_size_manual(values=c(2, 4)) +
  scale_color_manual(values = c('#117733', '#999933'))+
  scale_alpha_manual(values=c(0.5, 1)) +
  theme_classic() +
  theme(legend.position = c(0.84, 0.18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_x_continuous(breaks = seq(0,24,6))+
  labs(title=element_blank(),
       x='Hour',
       y='Substrate C (mol C / kg)')

print(tnc_diurnal_plot)
ggsave(file.path(FIGURE_DIR, 'Fig-substrate-diurnal.png'),
       plot = tnc_diurnal_plot,
       width = 9,
       height = 3.5,
       units = "in",
       dpi = 600
)


# plot assimilation by layers
times <- c(186.5, 209.5, 236.5, 258.5)
doys <- c(186, 209, 236, 258)
layer_assim <- data.frame(DOY = numeric(),
                          layer_number = numeric(),
                          layer_assimilation = numeric())
for (t in 1:length(times)){
  time <- times[t]
  doy <- doys[t]
  idx <- which(r$fractional_doy==time)
  for (i in 0:9){
    total_assim <- (r[idx, paste0('sunlit_Assim_layer_', i)]*
                r[idx, paste0('sunlit_fraction_layer_', i)] +
                r[idx, paste0('shaded_Assim_layer_', i)]*
                r[idx, paste0('shaded_fraction_layer_', i)]) *
      r[idx, 'lai'] /10
    new_row <- data.frame(
      DOY = doy,
      layer_number = i,
      layer_assimilation = total_assim)
    layer_assim <- rbind(layer_assim, new_row)
  }
}

# ==============================================================================
# Figure S4 Layer assimilation rate
# ==============================================================================
layer_assim_plot <- ggplot(layer_assim, aes(x = layer_number, y = layer_assimilation)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ DOY, nrow = 1, scales = "fixed") +
  labs(
    x = "Layer Number (Top:1 - Bottom: 9)",
    y = expression(paste("Layer Assimilation (", mu, "mol m"^{-2}, "s"^{-1},")"))) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(color = "grey90")
  ) +
  scale_x_continuous(breaks = unique(layer_assim$layer_number),
                     labels = as.integer(unique(layer_assim$layer_number)))

# Display the plot
print(layer_assim_plot)
ggsave(file.path(FIGURE_DIR, 'FigS-layer_assim.png'),
       plot = layer_assim_plot,
       width = 10,
       height = 3,
       units = "in",
       dpi = 600
)
