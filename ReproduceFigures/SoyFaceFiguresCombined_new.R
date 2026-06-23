library(BioCro)
library(UTRSoybeanBML)
library(BioCroWater)
library(ggplot2)
library(grid)
library(gridExtra)
library(lattice)
library(dplyr)
library(grid)
library(tidyr)
library(tidyverse)
library(patchwork)

# Clear workspace
rm(list=ls())

# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Parameters
source('../ParameterOptimization/soybean_parameter_expansion.R')
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')

# Source mse calculation and partitioning plotting functions
source('calculate_mse.R')
source('plot_partitioning.R')
source('plot_lai_comparison.R')
source('plot_biomass.R')
source('plot_sla.R')
RUN_SENSITIVITY <- FALSE # This takes a long time, so turn it off if not needed
PLOT_SLA <- TRUE

# years, sowing dates, and harvesting dates of growing seasons being fit to
years <- c('2002', '2004', '2005', '2006')
sow.date <- c(152, 149, 148, 148)
harv.date <- c(288, 289, 270, 270)
co2_opt = '_ambient_'
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
parameters$Catm      <- 372
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
RootVals <- list()
numrows <- vector()
figs <- list()
figs.elevCO2 <- list()
lai.figs <- list()
allocation.figs <- list()
allocation.elevCO2.figs <- list()

pod_sensitivity_lower <- NULL  # for -10%
pod_sensitivity_upper <- NULL  # for +10%

for (i in 1:length(years)) {   
  yr <- years[i]
  weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))
  sd.idx <- which(weather$doy == sow.date[i])[12]
  hd.ind <- which(weather$doy == harv.date[i])[24]
  
  weather.growingseason[[i]] <- weather[sd.idx: hd.ind,]
  
  ExpBiomass[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass.csv'))
  colnames(ExpBiomass[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  ExpBiomass.std[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
  colnames(ExpBiomass.std[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  RootVals[[i]] <- data.frame("DOY"=ExpBiomass[[i]]$DOY[3], "Root"=0.17*sum(ExpBiomass[[i]][5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  numrows[i] <- nrow(weather.growingseason[[i]])
  invwts <- ExpBiomass.std[[i]]
  
  # Load LAI data
  if(i>1){
    LAI[[i]]<-read.csv(file=paste0('../Data/SoyFACE_data/',yr,'_ambient_lai.csv'))
    LAI.elevCO2[[i]]<-read.csv(file=paste0('../Data/SoyFACE_data/',yr,'_elevated_lai.csv'))
  }
  
  soybean_optsolver[[i]] <- partial_run_biocro(initial_values,
                                               parameters,
                                               weather.growingseason[[i]],
                                               direct_modules,
                                               differential_modules,
                                               solver,
                                               arg_names,
                                               verbose = FALSE)
  
  result <- soybean_optsolver[[i]](optim_params_conversion(optim_params_short_SoyFACE)) 
  results[[i]] <- result
  
  # figure out the DVI when shell stops growing
  final_shell <- tail(ExpBiomass[[i]]$Pod,1) - tail(ExpBiomass[[i]]$Seed, 1) 
  print(final_shell)
  print(result$DVI[which.min(abs(result$Pod-final_shell))])
  
  # plot simulation vs observation biomass and calculate mse
  figs[[i]] <- plot_SoyFACE_biomass(result, ExpBiomass[[i]], ExpBiomass.std[[i]], co2_opt, years[i])
  
  # plot partitioning
  allocation.figs[[i]] <- plot_partitioning(result, years[i])
  
  # calculate sensitivity
  if (RUN_SENSITIVITY) {source('calculate_sensitivity.R')}
  
}

if (RUN_SENSITIVITY) {print(pod_sensitivity_lower)}

# initialize lists for figures
co2_opt = '_co2_'
parameters$Catm      <- 550

# Initialize lists
soybean_optsolver <- list()
RootVals <- list()
numrows <- vector()

for (i in 1:length(years)) {  
  yr <- years[i]
  
  ExpBiomass.elevCO2[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass.csv'))
  colnames(ExpBiomass.elevCO2[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  ExpBiomass.elevCO2.std[[i]] <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
  colnames(ExpBiomass.elevCO2.std[[i]])<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  RootVals[[i]] <- data.frame("DOY"=ExpBiomass.elevCO2[[i]]$DOY[3], "Root"=0.17*sum(ExpBiomass.elevCO2[[i]][5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  numrows[i] <- nrow(weather.growingseason[[i]])
  invwts <- ExpBiomass.std[[i]]

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
  
  # plot simulation vs observation biomass and calculate mse
  figs.elevCO2[[i]] <- plot_SoyFACE_biomass(result, ExpBiomass.elevCO2[[i]], ExpBiomass.elevCO2.std[[i]], co2_opt, years[i])
  
  # plot partitioning
  allocation.elevCO2.figs[[i]] <- plot_partitioning(result, years[i])
  
  # calculate sensitivity
  if (RUN_SENSITIVITY) {source('calculate_sensitivity.R')}
  
}

if (RUN_SENSITIVITY) {print(pod_sensitivity_lower)}


#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

common_legend <- g_legend(figs[[1]])
heights <- c(0.8, 1)
combined_graph_SoyFACE <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
                                  arrangeGrob(
                                    arrangeGrob(arrangeGrob(figs[[1]] + theme(axis.title.x = element_blank(),
                                                                              axis.text.x = element_blank(),
                                                                              legend.position="none"),
                                                            figs[[3]] + theme(legend.position="none"),
                                                            ncol = 2, left = '(A)'),
                                                top = 'Training'),
                                    
                                    arrangeGrob(arrangeGrob(figs[[2]] + theme(axis.title.x = element_blank(),
                                                                              axis.text.x = element_blank(),
                                                                              legend.position="none"),
                                                            figs[[4]] + theme(legend.position="none"),
                                                            nrow = 2, 
                                                            heights = heights, top = 'Ambient CO2'),
                                    
                                                arrangeGrob(figs.elevCO2[[1]] + theme(axis.title.x = element_blank(),
                                                                                      axis.text.x = element_blank(),
                                                                                      legend.position="none"),
                                                            figs.elevCO2[[2]] + theme(legend.position="none"),
                                                            nrow = 2, 
                                                            heights = heights, top = 'Elevated CO2'),
                                                arrangeGrob(figs.elevCO2[[3]] + theme(axis.title.x = element_blank(),
                                                                                      axis.text.x = element_blank(),
                                                                                      legend.position="none"),
                                                            figs.elevCO2[[4]] + theme(legend.position="none"),
                                                            nrow = 2,
                                                            heights = heights, top = 'Elevated CO2'),
                                                ncol = 3, left = '(B)'), 
                                    nrow = 2, heights = c(1.3,2)),
                                  common_legend,
                                  ncol=3, widths=c(0.3, 10, 1.5))

# Plot lai
for (i in 2:4){
  # results[[i]]$lai <- results[[i]]$lai * 1.12
  # results.elevCO2[[i]]$lai <- results.elevCO2[[i]]$lai * 1.12
  lai.figs[[i-1]] <- plot_amb_elev_lai(results[[i]], results.elevCO2[[i]], years[[i]], LAI[[i]], LAI.elevCO2[[i]])
  # print(lai.figs[[i-1]])
}

combined_graph.lai <- grid.arrange(arrangeGrob(textGrob(bquote("LAI"~(m^2~"/"~m^2)), rot = 90, gp=gpar(fontsize=12))),
                                   arrangeGrob(arrangeGrob(lai.figs[[1]] + theme(legend.position="none"),
                                                           lai.figs[[2]] + theme(legend.position="none"),
                                                           lai.figs[[3]] + theme(legend.position="none"),
                                                           ncol = 3),
                                               ncol = 1),
                                   ncol=2, widths=c(0.3, 5))

ggsave('FigS1-soyface-lai.png', 
       plot = combined_graph.lai, 
       width = 8,
       height = 2,
       units = "in",
       dpi = 600
)

######### Plot Root-Shoot Ratio #########
# Build combined data frame across all years
root_shoot_all <- do.call(rbind, lapply(1:length(years), function(i) {
  data.frame(
    year          = years[i],
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
root_shoot_ratio_plot <- ggplot(root_shoot_long, aes(x = fractional_doy, y = ratio, colour = treatment)) +
  geom_line() +
  facet_wrap(~ year, ncol = 2) +
  scale_colour_manual(values = c("Ambient CO2" = "steelblue", "Elevated CO2" = "tomato")) +
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
ggsave('Fig-root-shoot-ratio.png', 
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

ggsave('Fig-root-shoot-ratio-mean.png',
       plot   = root_shoot_mean_plot,
       width  = 4,
       height = 3,
       units  = "in",
       dpi    = 600)

######### Plot Harvest Index #########
# Build combined data frame across all years
root_shoot_all <- do.call(rbind, lapply(1:length(years), function(i) {
  data.frame(
    year           = years[i],
    fractional_doy = results[[i]]$fractional_doy,
    DVI            = results[[i]]$DVI,
    ambient_hi  = with(results[[i]],         0.761 * Pod / (Leaf + Stem + Pod)),
    elevated_hi = with(results.elevCO2[[i]], 0.761 * Pod / (Leaf + Stem + Pod))
  )
}))

# Pivot to long format for ggplot
root_shoot_long <- root_shoot_all %>%
  pivot_longer(
    cols      = c(ambient_hi, elevated_hi),
    names_to  = "co2_level",
    values_to = "HI"
  ) %>%
  mutate(co2_level = recode(co2_level,
                            ambient_hi  = "Ambient CO2",
                            elevated_hi = "Elevated CO2"
  ))

# Plot
hi_plot <- ggplot(root_shoot_long, aes(x = fractional_doy, y = HI, colour = co2_level)) +
  geom_line() +
  facet_wrap(~ year, ncol = 2) +
  scale_colour_manual(values = c("Ambient CO2" = "steelblue", "Elevated CO2" = "tomato")) +
  labs(
    x      = "Day of Year",
    y      = "Harvest Index",
    colour = "CO2 Level"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )
print(hi_plot)
ggsave('Fig-harvest-index.png', 
       plot = hi_plot, 
       width = 6,
       height = 4,
       units = "in",
       dpi = 600
)


####### Plot the shoot utilization ratio ########
# Build combined data frame across all years
utilization_all <- do.call(rbind, lapply(1:4, function(i) {
  
  # Ambient: daily sums then ratio
  ambient_daily <- results[[i]] %>%
    group_by(doy) %>%
    summarise(
      shoot_utilization   = sum(Leaf_utilization_rate + Stem_utilization_rate + Pod_utilization_rate),
      total_assimilation  = sum(canopy_assimilation_rate),
      .groups = "drop"
    ) %>%
    mutate(ratio = shoot_utilization / total_assimilation,
           treatment = "Ambient CO2")
  
  # Elevated: daily sums then ratio
  elevated_daily <- results.elevCO2[[i]] %>%
    group_by(doy) %>%
    summarise(
      shoot_utilization   = sum(Leaf_utilization_rate + Stem_utilization_rate + Pod_utilization_rate),
      total_assimilation  = sum(canopy_assimilation_rate),
      .groups = "drop"
    ) %>%
    mutate(ratio = shoot_utilization / total_assimilation,
           treatment = "Elevated CO2")
  
  rbind(ambient_daily, elevated_daily) %>%
    mutate(year = years[i])
}))

# Plot
shoot_utilization_assimilation_ratio_plot <- ggplot(utilization_all, aes(x = doy, y = ratio, colour = treatment)) +
  geom_line() +
  facet_wrap(~ year, ncol = 2) +
  scale_colour_manual(values = c("Ambient CO2" = "steelblue", "Elevated CO2" = "tomato")) +
  labs(
    x      = "Day of Year",
    y      = "Daily Shoot Utilization / Daily Canopy Assimilation",
    colour = "Treatment",
    title  = "Daily Shoot Utilization-to-Assimilation Ratio: Ambient vs Elevated CO2"
  ) +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  )

ggsave('FigS-shoot-utilization-assimilation-ratio.png', 
       plot = shoot_utilization_assimilation_ratio_plot, 
       width = 8,
       height = 6,
       units = "in",
       dpi = 600
)
# Build combined data frame across all years
utilization_all <- do.call(rbind, lapply(1:4, function(i) {
  
  # Ambient: daily sums then ratio
  ambient_daily <- results[[i]] %>%
    group_by(doy) %>%
    summarise(
      shoot_utilization   = sum(Leaf_utilization_rate + Stem_utilization_rate + Pod_utilization_rate),
      total_assimilation  = sum(canopy_assimilation_rate),
      .groups = "drop"
    ) %>%
    mutate(ratio = shoot_utilization / total_assimilation,
           treatment = "Ambient CO2")
  
  # Elevated: daily sums then ratio
  elevated_daily <- results.elevCO2[[i]] %>%
    group_by(doy) %>%
    summarise(
      shoot_utilization   = sum(Leaf_utilization_rate + Stem_utilization_rate + Pod_utilization_rate),
      total_assimilation  = sum(canopy_assimilation_rate),
      .groups = "drop"
    ) %>%
    mutate(ratio = shoot_utilization / total_assimilation,
           treatment = "Elevated CO2")
  
  rbind(ambient_daily, elevated_daily) %>%
    mutate(year = years[i])
}))

# Plot
shoot_utilization_assimilation_ratio_plot <- ggplot(utilization_all, aes(x = doy, y = ratio, colour = treatment)) +
  geom_line() +
  facet_wrap(~ year, ncol = 2) +
  scale_colour_manual(values = c("Ambient CO2" = "steelblue", "Elevated CO2" = "tomato")) +
  labs(
    x      = "Day of Year",
    y      = "Daily Shoot Utilization / Daily Canopy Assimilation",
    colour = "Treatment",
    title  = "Daily Shoot Utilization-to-Assimilation Ratio: Ambient vs Elevated CO2"
  ) +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  )

ggsave('FigS-shoot-utilization-assimilation-ratio.png', 
       plot = shoot_utilization_assimilation_ratio_plot, 
       width = 8,
       height = 6,
       units = "in",
       dpi = 600
)
# Build difference data frame across all years
ratio_diff_all <- do.call(rbind, lapply(1:4, function(i) {
  
  # Ambient: daily sums then ratio
  ambient_daily <- results[[i]] %>%
    group_by(doy) %>%
    summarise(
      ratio_ambient = sum(Leaf_utilization_rate + Stem_utilization_rate + Pod_utilization_rate) /
        sum(canopy_assimilation_rate),
      .groups = "drop"
    )
  
  # Elevated: daily sums then ratio
  elevated_daily <- results.elevCO2[[i]] %>%
    group_by(doy) %>%
    summarise(
      ratio_elevated = sum(Leaf_utilization_rate + Stem_utilization_rate + Pod_utilization_rate) /
        sum(canopy_assimilation_rate),
      .groups = "drop"
    )
  
  # Join and compute difference (ambient - elevated)
  inner_join(ambient_daily, elevated_daily, by = "doy") %>%
    mutate(ratio_diff = ratio_ambient - ratio_elevated,
           year = years[i])
}))

# Plot
ggplot(ratio_diff_all, aes(x = doy, y = ratio_diff)) +
  geom_line(colour = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  coord_cartesian(ylim = c(-1, 1)) +
  facet_wrap(~ year, ncol = 2) +
  labs(
    x     = "Day of Year",
    y     = "Ratio Difference (Ambient − Elevated)",
    title = "Difference in Daily Shoot Utilization-to-Assimilation Ratio: Ambient − Elevated CO2"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  )

###### Plot Transport Ratio #######
# Build combined data frame across all years
transportation_all <- do.call(rbind, lapply(1:4, function(i) {
  
  # Ambient: daily sums then ratio
  ambient_daily <- results[[i]] %>%
    group_by(doy) %>%
    summarise(
      Leaf_to_Stem   = sum(substrate_transport_Leaf_to_Stem),
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
      Leaf_to_Stem   = sum(substrate_transport_Leaf_to_Stem),
      Stem_to_Root  = sum(substrate_transport_Stem_to_Root),
      total_assimilation  = sum(canopy_assimilation_rate),
      .groups = "drop"
    ) %>%
    mutate(ratio = Stem_to_Root / total_assimilation,
           treatment = "Elevated CO2")
  
  rbind(ambient_daily, elevated_daily) %>%
    mutate(year = years[i])
}))

# Plot
transport_assimilation_ratio_plot <- ggplot(transportation_all, aes(x = doy, y = ratio, colour = treatment)) +
  geom_line() +
  facet_wrap(~ year, ncol = 2) +
  scale_colour_manual(values = c("Ambient CO2" = "steelblue", "Elevated CO2" = "tomato")) +
  labs(
    x      = "Day of Year",
    y      = "Daily Stem to Root Transport / Daily Canopy Assimilation",
    colour = "Treatment" #,
    # title  = "Daily Stem to Root Transport-to-Assimilation Ratio: Ambient vs Elevated CO2"
  ) +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold")
  )

ggsave('FigS-transport-assimilation-ratio.png', 
       plot = transport_assimilation_ratio_plot, 
       width = 8,
       height = 5,
       units = "in",
       dpi = 600
)

###### Plot Utilization Saturation Ratio #######
# Build combined data frame across all years
leaf_saturation_all <- do.call(rbind, lapply(1:4, function(i) {
  
  data.frame(
    fractional_doy = results[[i]]$fractional_doy,
    ambient        = with(results[[i]],
                          Leaf_utilization_rate / Leaf_structural_carbon / parameters$Leaf_utilization_rate_constant),
    elevated       = with(results.elevCO2[[i]],
                          Leaf_utilization_rate / Leaf_structural_carbon / parameters$Leaf_utilization_rate_constant),
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
utilization_saturation_plot <- ggplot(leaf_saturation_all, aes(x = fractional_doy, y = saturation,
                                colour = treatment, linewidth = treatment, alpha = treatment)) +
  geom_line() +
  facet_wrap(~ year, ncol = 2) +
  scale_colour_manual(values   = c("Ambient CO2" = "steelblue", "Elevated CO2" = "tomato")) +
  scale_linewidth_manual(values = c("Ambient CO2" = 0.8,         "Elevated CO2" = 0.4)) +
  scale_alpha_manual(values     = c("Ambient CO2" = 1,           "Elevated CO2" = 0.5)) +
  labs(x = "Day of Year", y = "Leaf Utilization Rate Saturation",
       colour = "Treatment", linewidth = "Treatment", alpha = "Treatment") +
  theme_bw() +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold"))

ggsave('FigS-utilization_saturation.png', 
       plot = utilization_saturation_plot, 
       width = 8,
       height = 6,
       units = "in",
       dpi = 600
)

# Compute difference before pivoting
leaf_saturation_diff <- do.call(rbind, lapply(1:4, function(i) {
  data.frame(
    fractional_doy = results[[i]]$fractional_doy,
    diff = with(results[[i]],
                Leaf_utilization_rate / Leaf_structural_carbon / parameters$Leaf_utilization_rate_constant) -
      with(results.elevCO2[[i]],
           Leaf_utilization_rate / Leaf_structural_carbon / parameters$Leaf_utilization_rate_constant),
    year = years[i]
  )
}))

ggplot(leaf_saturation_diff, aes(x = fractional_doy, y = diff)) +
  geom_line(colour = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ year, ncol = 2) +
  labs(x = "Day of Year", y = "Saturation Difference (Ambient − Elevated)",
       title = "Difference in Leaf Utilization Rate Saturation") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"))
################################# Plot shade treatment ##################################
source('plot_shade_R5.R')
################################# Plot pod removal treatment ############################
source('plot_depod_R5.R')
################################# Plot defoliation treatment ############################
source('plot_defoliation_parvej2025.R')
#################################### Plot hail event ####################################
source('plot_hail.R')
############################################# Plot LD11
years <- c('2021', '2022', '2023', '2024')
Catms <- c(414.7, 417.2, 419.3, 422.8) # from NOAA

updated_parameters <- read.csv('../Data/Soybean-BioCro_Parameters/updated-ld11-parameters-2022.csv')
updated_idx <- match(updated_parameters$New.variable.name, names(parameters))
parameters[updated_idx] <- updated_parameters$LD11

ExpBiomass <- list()
weather.afteremergence <- list()
ld11.results <- list()
ld11.figs <- list()
lai.ld11.figs <- list()
allocation.ld11.figs <- list()

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

ratio_df <- list()

co2_opt <- '_ambient_'
for (i in 1:length(years)){
  yr <- years[i]
  ExpBiomass[[i]] <- loadRData(paste0('../../energy-farm-biocro/soybean_ld11_biomass_', years[i],'/soybean_ld11_biomass_', years[i], '.RData'))
  sow.time <- min(ExpBiomass[[i]]$time)
  harv.time <- max(ExpBiomass[[i]]$time)
  weather <- loadRData(paste0('../../energy-farm-biocro/weather_', years[i], '/weather', years[i], '_hourly.RData'))
  parameters$Catm <- Catms[i]
  # update initial values
  sub_frac <- 0.1           # substrate_fraction
  str_frac <- 1 - sub_frac  # structural_fraction
  seed_mass <- ExpBiomass[[i]]$initial_seed[1]
  j <- 2
  mass_t <- sum(ExpBiomass[[i]][j, c('leaf', 'stem', 'root')])
  leaf_frac <- ExpBiomass[[i]]$leaf[j]/mass_t
  stem_frac <- ExpBiomass[[i]]$stem[j]/mass_t
  root_frac <- ExpBiomass[[i]]$root[j]/mass_t
  
  leaf_frac <- 0.45
  stem_frac <- 0.25
  root_frac <- 0.3
  
  print(paste0('lsr fractions:', round(leaf_frac,2), ', ', round(stem_frac, 2), ', ' , round(root_frac,2)))
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
  # parameters$par_energy_content <- 0.235
  
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
  if (RUN_SENSITIVITY) {source('calculate_sensitivity.R')}
}

if (RUN_SENSITIVITY) {print(pod_sensitivity_lower)}

### Plot LAI ###
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

ggsave('FigS1-ld11-lai.png', 
       plot = ld11_lai_comparison, 
       width = 8,
       height = 6,
       units = "in",
       dpi = 600
)

# extract the common legend
g_legend <-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

common_legend <- g_legend(ld11.figs[[1]])

combined_graph <- grid.arrange(arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90, gp=gpar(fontsize=12))),
                               arrangeGrob(arrangeGrob(ld11.figs[[1]] + theme(legend.position="none"), top = "Testing"),
                                           arrangeGrob(ld11.figs[[2]] + theme(legend.position="none"), top = "Testing"),
                                           arrangeGrob(ld11.figs[[3]] + theme(legend.position="none"), top = "Testing"),
                                           arrangeGrob(ld11.figs[[4]] + theme(legend.position="none"), top = "Testing"),
                                           ncol = 4, top = "LD11 at Energy Farm (Ambient CO2)"),
                               common_legend, 
                               ncol=3, widths=c(0.3, 5, 1.1))

common_legend <- g_legend(allocation.ld11.figs[[1]])

combined_graph.allocation <- grid.arrange(arrangeGrob(textGrob('Allocation %', rot = 90, gp=gpar(fontsize=12))),
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

ggsave('FigS-allocation-netimport.png', 
       plot = combined_graph.allocation, 
       width = 10,
       height = 7,
       units = "in",
       dpi = 600
)


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

combined_graph_v3 <- grid.arrange(arrangeGrob(
                                    arrangeGrob(textGrob('Biomass (Mg / ha)', rot = 90)),
                                    arrangeGrob(
                                      arrangeGrob(arrangeGrob(figs[[1]] + 
                                                                theme(legend.position="none"),
                                                              figs[[3]] + theme(legend.position="none"),
                                                              ncol = 2),
                                                  top = textGrob("Parameterization - Pioneer 93B15 at Ambient CO2", 
                                                                 gp=gpar(fontface="bold", fontsize=12))),
                                      textGrob(""), # spacer
                                      arrangeGrob(arrangeGrob(figs[[2]] + 
                                                                theme(axis.title.x = element_blank(),
                                                                                axis.text.x = element_blank(),
                                                                                legend.position="none"),
                                                              figs[[4]] + theme(legend.position="none"),
                                                              nrow = 2, 
                                                              heights = heights, top = 'Ambient CO2'),
                                                  
                                                  arrangeGrob(figs.elevCO2[[1]] + theme(axis.title.x = element_blank(),
                                                                                        axis.text.x = element_blank(),
                                                                                        legend.position="none"),
                                                              figs.elevCO2[[2]] + theme(legend.position="none"),
                                                              nrow = 2, 
                                                              heights = heights, top = 'Elevated CO2'),
                                                  arrangeGrob(figs.elevCO2[[3]] + theme(axis.title.x = element_blank(),
                                                                                        axis.text.x = element_blank(),
                                                                                        legend.position="none"),
                                                              figs.elevCO2[[4]] + theme(legend.position="none"),
                                                              nrow = 2,
                                                              heights = heights, top = 'Elevated CO2'),
                                                  ncol = 3, 
                                                  top = textGrob("Pioneer 93B15 at SoyFACE", 
                                                                 gp=gpar(fontface="bold", fontsize=12))), 
                                      textGrob(""), # spacer
                                      arrangeGrob(arrangeGrob(ld11.figs[[1]] + 
                                                                theme(axis.title.x = element_blank(),
                                                                                        axis.text.x = element_blank(),
                                                                                        legend.position="none"),
                                                              ld11.figs[[3]] + theme(legend.position="none"),
                                                              nrow = 2, 
                                                              heights = heights),
                                                  arrangeGrob(ld11.figs[[2]] + theme(axis.title.x = element_blank(),
                                                                                        axis.text.x = element_blank(),
                                                                                        legend.position="none"),
                                                              ld11.figs[[4]] + theme(legend.position="none"),
                                                              nrow = 2,
                                                              heights = heights),
                                                  ncol = 2, 
                                                  top = textGrob("LD11-2170 at Ambient CO2", 
                                                                 gp=gpar(fontface="bold", fontsize=12))),
                                      
                                      nrow = 5, heights = c(1.3, 0.1, 2, 0.1, 2.4)),
                                    ncol = 2, widths = c(0.3, 10)),
                                  common_legend_horizontal,
                                  nrow = 2, heights = c(11, 0.4)  # adjust 1 to give legend more/less space
                                  )


# Get the total height of the drawn object to calculate y positions
total_h <- 11 + 0.4  # matches heights above
plot_fraction <- 11 / total_h  # the fraction of the figure that is plots

sections_h <- 1.3 + 0.1 + 2 + 0.1 + 2.4 

y_A <- 0.996
y_B <- 0.996 - (1.3 + 0.1) / sections_h * plot_fraction
y_C <- 0.996 - (1.3 + 0.1 + 2 + 0.1) / sections_h * plot_fraction

png("Fig2-biomass.png", width = 5.5, height = 11, units = "in", res = 600)
grid.newpage()
grid.draw(combined_graph_v3)  # draw the plot FIRST

# overlay the labels on top
grid.text("(A)", x = 0.03, y = y_A, just = c("left", "top"),
          gp = gpar(fontface = "bold", fontsize = 12))
grid.text("(B)", x = 0.03, y = y_B, just = c("left", "top"),
          gp = gpar(fontface = "bold", fontsize = 12))
grid.text("(C)", x = 0.03, y = y_C, just = c("left", "top"),
          gp = gpar(fontface = "bold", fontsize = 12))

dev.off()

combined_graph_v4 <- grid.arrange(
  arrangeGrob(
    arrangeGrob(
      textGrob('Biomass (Mg / ha)', rot = 90), # Col1: y-axis Grob
      arrangeGrob( # Col 2: Biomass figures
        arrangeGrob( # Col 2, Row 1: Pioneer 93B15 Grob
          arrangeGrob( # ambient CO2 Grob
            arrangeGrob( # Parameterization Grob
              figs[[1]] + theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank(),
                                legend.position="none"),
              figs[[3]] + theme(legend.position="none", 
                                axis.title.x = element_blank()),
              nrow = 2),
            # top = textGrob("Parameterization", 
            #                gp=gpar(fontface="bold", fontsize=12))),
            arrangeGrob( # ambient CO2 2nd column Grob
              figs[[2]] + theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank(),
                                legend.position="none"),
              figs[[4]] + theme(legend.position="none",
                                axis.title.x = element_blank()),
              nrow = 2), 
            ncol = 2, 
            top = textGrob("(A) Ambient CO2", 
                           gp=gpar(fontface="bold", fontsize=12))
            ),
          arrangeGrob( # elevated CO2 Grob
            arrangeGrob( # eCO2 col 1
              figs.elevCO2[[1]] + theme(axis.title.x = element_blank(),
                                        axis.text.x = element_blank(),
                                        legend.position="none"),
              figs.elevCO2[[3]] + theme(axis.title.x = element_blank(),
                                        legend.position="none"),
              nrow = 2), 
            arrangeGrob(figs.elevCO2[[2]] + theme(axis.title.x = element_blank(),
                                                  axis.text.x = element_blank(),
                                                  legend.position="none"),
                        figs.elevCO2[[4]] + theme(axis.title.x = element_blank(),
                                                  legend.position="none"),
                        nrow = 2),
            ncol = 2, 
            top = textGrob("(B) Elevated CO2", 
                           gp=gpar(fontface="bold", fontsize=12))
            ),
          ncol = 2,
          top = textGrob("Pioneer 93B15", 
                         gp=gpar(fontface="bold", fontsize=12))
          ), # end of the Pioneer 93B15 Grob
        
        textGrob(""), # Row 2: spacer
        
        arrangeGrob( # Row 3: LD11 figures
          ld11.figs[[1]] + theme(legend.position="none"),
          ld11.figs[[2]] + theme(legend.position="none"),
          ld11.figs[[3]] + theme(legend.position="none"),
          ld11.figs[[4]] + theme(legend.position="none"),
          ncol = 4,
          top = textGrob("(C) LD11-2170 at Ambient CO2", 
                         gp=gpar(fontface="bold", fontsize=12))
        ),
        nrow = 3, heights = c(1.8, 0.1, 1)
      ),
      ncol = 2, widths = c(0.3, 10)
    ),
    common_legend_horizontal,
    nrow = 2, heights = c(11, 0.5))
)

combined_graph_v5 <- grid.arrange(
  arrangeGrob(
    arrangeGrob(
      textGrob('Biomass (Mg / ha)', rot = 90), # Col1: y-axis Grob
      arrangeGrob( # Col 2: Biomass figures
        arrangeGrob( # Col 2, Row 1: Pioneer 93B15 ambient CO2 Grob
              figs[[1]] + theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank(),
                                legend.position="none")
                        + labs(title = '(A) 2002'),
              figs[[2]] + theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank(),
                                legend.position="none")
                        + labs(title = '(B) 2004'),
              figs[[3]] + theme(legend.position="none", 
                                axis.title.x = element_blank(),
                                axis.text.x = element_blank())
                        + labs(title = '(C) 2005'),
              figs[[4]] + theme(legend.position="none",
                                axis.title.x = element_blank(),
                                axis.text.x = element_blank())
                        + labs(title = '(D) 2006'),
              ncol = 4,
              top = textGrob("Pioneer 93B15 at Ambient CO2", 
                           gp=gpar(fontface="bold", fontsize=12))
        ),
        textGrob(""), # Row 2: spacer,
        arrangeGrob( # Row 3: Pioneer 93B15 ambient CO2 Grob
          figs.elevCO2[[1]] + theme(axis.title.x = element_blank(),
                            axis.text.x = element_blank(),
                            legend.position="none")
                            + labs(title = '(E) 2002'),
          figs.elevCO2[[2]] + theme(axis.title.x = element_blank(),
                            axis.text.x = element_blank(),
                            legend.position="none")
                            + labs(title = '(F) 2004'),
          figs.elevCO2[[3]] + theme(legend.position="none", 
                            axis.title.x = element_blank(),
                            axis.text.x = element_blank())
                            + labs(title = '(G) 2005'),
          figs.elevCO2[[4]] + theme(legend.position="none",
                            axis.title.x = element_blank(),
                            axis.text.x = element_blank())
                            + labs(title = '(H) 2006'),
          ncol = 4,
          top = textGrob("Pioneer 93B15 at Elevated CO2", 
                       gp=gpar(fontface="bold", fontsize=12))
        ),
        textGrob(""), # Row 4: spacer,
        arrangeGrob( # Row 5: LD11 figures
          ld11.figs[[1]] + theme(legend.position="none") + labs(title = '(I) 2021'),
          ld11.figs[[2]] + theme(legend.position="none") + labs(title = '(J) 2022'),
          ld11.figs[[3]] + theme(legend.position="none") + labs(title = '(K) 2023'),
          ld11.figs[[4]] + theme(legend.position="none") + labs(title = '(L) 2024'),
          ncol = 4,
          top = textGrob("LD11-2170 at Ambient CO2", 
                         gp=gpar(fontface="bold", fontsize=12))
        ),
        nrow = 5, heights = c(1, 0.05, 1, 0.05, 1)
      ),
      ncol = 2, widths = c(0.3, 10)
    ),
    common_legend_horizontal,
    nrow = 2, heights = c(11, 0.5))
)


ggsave('Fig-Biomass.png', 
       plot = combined_graph_v5, 
       width = 10,
       height = 8,
       units = "in",
       dpi = 600)

total_precip <- sapply(1:4, function(i) {
  ld11.results[[i]] %>%
    filter(DVI >= 0 & DVI <= 1) %>%
    summarise(total = sum(precip, na.rm = TRUE)) %>%
    pull(total)
})

names(total_precip) <- years
print(total_precip)


# Plot substrate C concentration
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
# y='Leaf Substrate C (mol glucose eq./Mg)')

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

tnc_growingseasion_plot <- grid.arrange(arrangeGrob(textGrob('Substrate C (mol C / kg)', rot = 90, gp=gpar(fontsize=12))),
                                             arrangeGrob(arrangeGrob(leaf_tnc_plot, top = '(A) Leaf'),
                                                         arrangeGrob(stem_tnc_plot, top = '(B) Stem'),
                                                                     ncol = 2), #, top = 'LD11-2170 (2022) Leaf and Stem Substrate C Concentrations'),
                                             ncol=2, widths=c(0.2, 10))

ggsave('Fig-substrate-annual.png', 
       plot = tnc_growingseasion_plot, 
       width = 10,
       height = 3.5,
       units = "in",
       dpi = 600
)

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
  # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  scale_x_continuous(breaks = seq(0,24,6))+
  labs(title=element_blank(), 
       x='Hour',
       y='Substrate C (mol C / kg)')

print(tnc_diurnal_plot)
ggsave('Fig5-substrate-diurnal.png', 
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
                          layer_assimilation = numeric())# ,
                          # assim_type = character())
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

layer_assim_plot <- ggplot(layer_assim, aes(x = layer_number, y = layer_assimilation))+ #, group = assim_type)) +
  geom_point() + # aes(color=assim_type)) +
  geom_line() + #aes(color=assim_type)) +
  facet_wrap(~ DOY, nrow = 1, scales = "fixed") +
  labs(# title = "Layer Assimilation by Layer Number on Different DOYs",
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
ggsave('Fig6-layer_assim.png', 
       plot = layer_assim_plot, 
       width = 10,
       height = 3,
       units = "in",
       dpi = 600
)

if (RUN_SENSITIVITY) {source('generate_sensitivity_latex_table.R')}


max(r$Leaf_senescence_loss/(1-(1-optim_params_short_SoyFACE[22])))+
  max(r$Stem_senescence_loss/(1-(1-optim_params_short_SoyFACE[22])))+
  max(r$Root_senescence_loss/(1-(1-optim_params_short_SoyFACE[22])))

r_summary <- r %>%
  group_by(doy) %>%
  summarise(
    diff = sum(substrate_transport_Leaf_to_Stem, na.rm = TRUE) -
      sum(substrate_transport_Stem_to_Pod, na.rm = TRUE)
  )

transport_diff_fig <- ggplot(r_summary, aes(x = doy, y = diff)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Difference in Substrate Transport by Day of Year",
    x = "2022 Day of Year (DOY)",
    y = expression(paste("Leaf to Stem - Stem to Pod", "(mol m"^{-2}, "day"^{-1},")"))
  ) +
  theme_minimal()
ggsave('FigS-transport_diff.png', 
       plot = transport_diff_fig, 
       width = 4,
       height = 3,
       units = "in",
       dpi = 600
)




