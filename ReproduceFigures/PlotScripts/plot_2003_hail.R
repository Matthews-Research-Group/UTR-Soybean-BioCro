# =============================================================================
# 2003 Hail Event Simulation (Ambient vs. Elevated CO2)
# -----------------------------------------------------------------------------
# Purpose:
#   To simulate the 2003 SoyFACE growing season with a mid-season hail event
#   (modeled as an abrupt loss of leaf/stem tissue) under both ambient and
#   elevated CO2, using two C allocation models ("UTR" and "Partitioning"),
#   and produce three figures:
#     - Fig-hail-utr.png          (combined_with_hail_utr)
#     - FigS-hail-partitioning.png (combined_with_hail_partitioning)
#     - FigS-hail-utilization.png  (hail_utr_utilization_plot)
# -----------------------------------------------------------------------------
# Assumes the following objects already exist in the environment (created
# upstream in the analysis pipeline, e.g. by a BioCro model-setup script):
#   - initial_values, parameters, direct_modules, differential_modules, solver
#   - arg_names, optim_params_short_SoyFACE, optim_params_conversion()
#   - plot_SoyFACE_biomass() -- plots simulated vs. observed biomass
# =============================================================================

library(gridExtra)
library(grid)

# -----------------------------------------------------------------------------
# update_differential_quantities()
# -----------------------------------------------------------------------------
# Reduces the leaf and stem differential-state variables to remaining fractions
# (40% of leaf, 50% of stem) to represent tissue loss from hail.
#
# Arguments:
#   r                 - data frame of simulation output up to (and including)
#                       the time of the hail event
#   updated_values    - names of the differential-state variables to specify
#   model             - which model's state layout to use: 'utr' or
#                       'partitioning'
#
# Returns:
#   A named list of state values to use as the starting point for the
#   post-hail simulation
# -----------------------------------------------------------------------------
update_differential_quantities <- function(r, updated_values, model){
  
  leaf_remaining_percentage <- 0.4
  stem_remaining_percentage <- 0.5
  
  last_row <- nrow(r)
  differential_quantities_just_before_hail <-
    as.list(r[last_row, updated_values])
  
  differential_quantities_just_after_hail <-
    differential_quantities_just_before_hail
  
  if (model == 'utr'){
    differential_quantities_just_after_hail$Leaf_substrate_carbon <-
      differential_quantities_just_before_hail$Leaf_substrate_carbon * 
      leaf_remaining_percentage
    differential_quantities_just_after_hail$Leaf_structural_carbon <-
      differential_quantities_just_before_hail$Leaf_structural_carbon * 
      leaf_remaining_percentage
    
    differential_quantities_just_after_hail$Stem_substrate_carbon <-
      differential_quantities_just_before_hail$Stem_substrate_carbon * 
      stem_remaining_percentage
    differential_quantities_just_after_hail$Stem_structural_carbon <-
      differential_quantities_just_before_hail$Stem_structural_carbon * 
      stem_remaining_percentage
  }
  else if (model == 'partitioning'){
    differential_quantities_just_after_hail$Leaf <-
      differential_quantities_just_before_hail$Leaf * 
      leaf_remaining_percentage
    
    differential_quantities_just_after_hail$Stem <-
      differential_quantities_just_before_hail$Stem * 
      stem_remaining_percentage
  }
  
  return(differential_quantities_just_after_hail)
}

# -----------------------------------------------------------------------------
# g_legend()
# -----------------------------------------------------------------------------
# Extracts the legend ("guide-box" grob) from a ggplot object.
# -----------------------------------------------------------------------------
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# -----------------------------------------------------------------------------
# make_strip_label()
# -----------------------------------------------------------------------------
# Builds a panel title grob with a left-aligned panel letter (e.g. "(A)") and
# a centered title (e.g. "Ambient CO2"), for use as the `top=` argument of
# arrangeGrob() in the UTR hail figure.
# -----------------------------------------------------------------------------
make_strip_label <- function(letter, title, fontsize = 12, height_lines = 1.5) {
  arrangeGrob(
    grobTree(
      textGrob(letter, x = 0.01, hjust = 0,
               gp = gpar(fontsize = fontsize)),
      textGrob(title,  x = 0.5,  hjust = 0.5,
               gp = gpar(fontsize = fontsize))
    ),
    heights = unit(height_lines, "lines")  # reserves explicit vertical space
  )
}

# Load 2003 weather and observed biomass data, and define the hail timing
yr <- '2003'
weather <- read.csv(file = paste0('../Data/Weather_data/', yr, '_Bondville_IL_daylength.csv'))

co2_opt <- '_ambient_'
parameters$Catm <- 372
ExpBiomass <- read.csv(file = paste0('../Data/SoyFACE_data/biomasses/', yr, co2_opt, 'biomass.csv'))
colnames(ExpBiomass) <- c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

ExpBiomass.std <- read.csv(file = paste0('../Data/SoyFACE_data/biomasses/', yr, co2_opt, 'biomass_std.csv'))
colnames(ExpBiomass.std) <- c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

sd.idx <- which(weather$doy == 147)[12]              # planting/start day (Morgan et al. 2005)
hd.ind <- which(weather$doy == max(ExpBiomass$DOY))[24]  # harvest/end day
hail.ind <- which(weather$doy == 198)[14]      # day of the hail event

# Split the growing-season weather into "before hail" and "after hail"
# segments (used for both the ambient and elevated-CO2 runs below).
weather.growingseason_1 <- weather[sd.idx:hail.ind, ]
weather.growingseason_2 <- weather[hail.ind:hd.ind, ]

soybean$parameters$time_zone_offset <- -6

# =============================================================================
# Ambient CO2: simulate the hail event
# =============================================================================

#### Before hail ####
soybean_utr_optsolver_hail1 <- partial_run_biocro(initial_values,
                                                  parameters,
                                                  weather.growingseason_1,
                                                  direct_modules,
                                                  differential_modules,
                                                  solver,
                                                  arg_names,
                                                  verbose = FALSE)
result_utr_hail1 <- soybean_utr_optsolver_hail1(
  optim_params_conversion(optim_params_short_SoyFACE))
result_partitioning_hail1 <- with(soybean, 
                                  run_biocro(initial_values,
                                             parameters,
                                             weather.growingseason_1,
                                             direct_modules,
                                             differential_modules,
                                             solver,
                                             verbose = FALSE))

#### After hail ####
# Reduce leaf/stem state per update_differential_quantities(), then carry
# the adjusted state forward as the starting point for the "after hail" run.

## UTR model ##
differential_quantities_just_after_hail_utr <- update_differential_quantities(
  result_utr_hail1, names(initial_values), 'utr')

soybean_utr_optsolver_hail2 <- 
  partial_run_biocro(differential_quantities_just_after_hail_utr,
                     parameters,
                     weather.growingseason_2,
                     direct_modules,
                     differential_modules,
                     solver,
                     arg_names,
                     verbose = FALSE)

orignal_utr_params <- data.frame(
  optim_params_conversion(optim_params_short_SoyFACE))
rownames(orignal_utr_params) <- arg_names
colnames(orignal_utr_params) <- "Value"
parameters_after_hail <- orignal_utr_params # Could be changed to test other hypothesis

result_utr_hail2 <- soybean_utr_optsolver_hail2(
  parameters_after_hail$Value)

## Partitioning model ##
differential_quantities_just_after_hail_partitioning <- 
  update_differential_quantities(result_partitioning_hail1,
                                 names(soybean$initial_values),
                                 'partitioning')

result_partitioning_hail2 <- 
  run_biocro(differential_quantities_just_after_hail_partitioning,
             soybean$parameters,
             weather.growingseason_2,
             soybean$direct_modules,
             soybean$differential_modules,
             soybean$ode_solver,
             verbose = FALSE)

# Connect the "before" and "after" segments back into one data frame
result_utr_hail <- rbind(result_utr_hail1[seq_len(nrow(result_utr_hail1) - 1), ], 
                         result_utr_hail2)
result_partitioning_hail <- 
  rbind(result_partitioning_hail1[seq_len(nrow(result_partitioning_hail1) - 1), ], 
        result_partitioning_hail2)
result_partitioning_hail$Pod <- result_partitioning_hail$Grain + 
  result_partitioning_hail$Shell

# Simulated-vs-observed biomass at ambient CO2.
fig_2003_utr_with_hail <- plot_SoyFACE_biomass(
  result_utr_hail, ExpBiomass, ExpBiomass.std, co2_opt, yr)
fig_2003_utr_with_hail <- fig_2003_utr_with_hail +
  geom_vline(xintercept = 198, color = "red", alpha = 0.5, linewidth = 1)

fig_2003_partitioning_with_hail <- plot_SoyFACE_biomass(
  result_partitioning_hail, ExpBiomass, ExpBiomass.std, co2_opt, yr)
fig_2003_partitioning_with_hail <- fig_2003_partitioning_with_hail +
  geom_vline(xintercept = 198, color = "red", alpha = 0.5, linewidth = 1)


# =============================================================================
# Elevated CO2: repeat the hail simulation at Catm = 550 ppm
# =============================================================================
co2_opt <- '_CO2_'
parameters$Catm <- 550
soybean$parameters$Catm <- 550

ExpBiomass.elevCO2 <- read.csv(
  file = paste0('../Data/SoyFACE_data/biomasses/', yr, co2_opt, 'biomass.csv'))
colnames(ExpBiomass.elevCO2) <- c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

ExpBiomass.elevCO2.std <- read.csv(
  file = paste0('../Data/SoyFACE_data/biomasses/', yr, co2_opt, 'biomass_std.csv'))
colnames(ExpBiomass.elevCO2.std) <- c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")

#### Before hail ####
soybean_utr_optsolver_eCO2_hail1 <- partial_run_biocro(initial_values,
                                                       parameters,
                                                       weather.growingseason_1,
                                                       direct_modules,
                                                       differential_modules,
                                                       solver,
                                                       arg_names,
                                                       verbose = FALSE)
result_utr_eCO2_hail1 <- soybean_utr_optsolver_eCO2_hail1(
  optim_params_conversion(optim_params_short_SoyFACE))
result_partitioning_eCO2_hail1 <- with(soybean, run_biocro(initial_values,
                                                           parameters,
                                                           weather.growingseason_1,
                                                           direct_modules,
                                                           differential_modules,
                                                           solver,
                                                           verbose = FALSE))

#### After hail ####
## UTR model ##
differential_quantities_just_after_hail_utr_eCO2 <- 
  update_differential_quantities(result_utr_eCO2_hail1,
                                 names(initial_values),
                                 'utr')

soybean_utr_optsolver_eCO2_hail2 <- 
  partial_run_biocro(differential_quantities_just_after_hail_utr_eCO2,
                     parameters,
                     weather.growingseason_2,
                     direct_modules,
                     differential_modules,
                     solver,
                     arg_names,
                     verbose = FALSE)

result_utr_eCO2_hail2 <- soybean_utr_optsolver_eCO2_hail2(parameters_after_hail$Value)

## Partitioning model ##
differential_quantities_just_before_hail_partitioning_eCO2 <- 
  update_differential_quantities(result_partitioning_eCO2_hail1,
                                 names(soybean$initial_values),
                                 'partitioning')

result_partitioning_eCO2_hail2 <- 
  run_biocro(differential_quantities_just_before_hail_partitioning_eCO2,
             soybean$parameters,
             weather.growingseason_2,
             soybean$direct_modules,
             soybean$differential_modules,
             soybean$ode_solver,
             verbose = FALSE)

# Connect the "before" and "after" segments back into one data frame.
result_utr_eCO2_hail <- rbind(result_utr_eCO2_hail1[
  seq_len(nrow(result_utr_eCO2_hail1) - 1), ], 
  result_utr_eCO2_hail2)
result_partitioning_eCO2_hail <- rbind(result_partitioning_eCO2_hail1[
                                        seq_len(nrow(result_partitioning_eCO2_hail1) - 1), ], 
                                       result_partitioning_eCO2_hail2)
result_partitioning_eCO2_hail$Pod <- result_partitioning_eCO2_hail$Grain + 
  result_partitioning_eCO2_hail$Shell

# Simulated-vs-observed biomass at elevated CO2 plots
fig_2003_utr_eCO2_with_hail <- plot_SoyFACE_biomass(result_utr_eCO2_hail, 
                                                    ExpBiomass.elevCO2, 
                                                    ExpBiomass.elevCO2.std, 
                                                    co2_opt, yr)

fig_2003_utr_eCO2_with_hail <- fig_2003_utr_eCO2_with_hail +
  geom_vline(xintercept = 198, color = "red", alpha = 0.5, linewidth = 1)


fig_2003_partitioning_eCO2_with_hail <- 
  plot_SoyFACE_biomass(result_partitioning_eCO2_hail, ExpBiomass.elevCO2, 
                       ExpBiomass.elevCO2.std, co2_opt, yr)

fig_2003_partitioning_eCO2_with_hail <- fig_2003_partitioning_eCO2_with_hail +
  geom_vline(xintercept = 198, color = "red", alpha = 0.5, linewidth = 1)

# Shared legend for the combined figures below
common_legend <- g_legend(fig_2003_utr_eCO2_with_hail)

# =============================================================================
# Figure: Hail simulation, UTR model (ambient vs. elevated CO2)
# =============================================================================
combined_with_hail_utr <- grid.arrange(
  arrangeGrob(textGrob(expression('Biomass (Mg ha'^{-1}*')'), rot = 90)),
  arrangeGrob(
    arrangeGrob(
      arrangeGrob(fig_2003_utr_with_hail +
                    theme(axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          legend.position = "none"),
                  top = make_strip_label("(A)", "Ambient CO2")),
      arrangeGrob(fig_2003_utr_eCO2_with_hail +
                    theme(axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          legend.position = "none"),
                  top = make_strip_label("(B)", "Elevated CO2")),
      ncol = 2),
    arrangeGrob(textGrob('Day of Year (2003)')),
    nrow = 2, heights = c(4, 0.3)),
  common_legend,
  ncol = 3, widths = c(0.3, 5, 1))

ggsave(file.path(FIGURE_DIR, 'Fig-hail-utr.png'),
       plot = combined_with_hail_utr,
       width = 6,
       height = 2.5,
       units = "in",
       dpi = 600)

# =============================================================================
# Figure: Hail simulation, Partitioning model (ambient vs. elevated CO2)
# =============================================================================
combined_with_hail_partitioning <- grid.arrange(
  arrangeGrob(textGrob(expression('Biomass (Mg ha'^{-1}*')'), rot = 90)),
  arrangeGrob(
    arrangeGrob(
      arrangeGrob(fig_2003_partitioning_with_hail +
                    theme(axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          legend.position = "none"),
                  top = make_strip_label("(A)", "Ambient CO2")),
      arrangeGrob(fig_2003_partitioning_eCO2_with_hail +
                    theme(axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          legend.position = "none"),
                  top = make_strip_label("(B)", "Elevated CO2")),
      ncol = 2),
    arrangeGrob(textGrob('Day of Year (2003)')),
    nrow = 2, heights = c(4, 0.3)),
  common_legend,
  ncol = 3, widths = c(0.3, 5, 1))

ggsave(file.path(FIGURE_DIR, 'FigS-hail-partitioning.png'),
       plot = combined_with_hail_partitioning,
       width = 6,
       height = 2.5,
       units = "in",
       dpi = 600)

# =============================================================================
# Figure: Leaf and stem utilization rates (UTR model only)
# =============================================================================
# Compares daily leaf/stem carbon-utilization rates before and after the hail
# (marked with a red shaded vertical band).
utr_ambient_df <- result_utr_hail %>%
  group_by(doy) %>%
  summarise(
    Leaf = sum(`Leaf_utilization_rate`, na.rm = TRUE),
    Stem = sum(`Stem_utilization_rate`, na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(cols = c(Leaf, Stem), 
                      names_to = "type", values_to = "value") %>%
  mutate(scenario = "Ambient")

utr_elev_df <- result_utr_eCO2_hail %>%
  group_by(doy) %>%
  summarise(
    Leaf = sum(`Leaf_utilization_rate`, na.rm = TRUE),
    Stem = sum(`Stem_utilization_rate`, na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(cols = c(Leaf, Stem), 
                      names_to = "type", values_to = "value") %>%
  mutate(scenario = "Elevated")

hail_utr_utilization_plot <-
  bind_rows(utr_ambient_df, utr_elev_df) %>%
  ggplot(aes(x = doy, y = value, color = type, linetype = scenario)) +
  annotate("rect", xmin = 197.5, xmax = 198.5, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.5) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("Leaf" = "#117733", "Stem" = "#999933"),
    name = "Utilization Rate"
  ) +
  scale_linetype_manual(
    values = c("Ambient" = "solid", "Elevated" = "dotted"),
    name = "CO2 Level"
  ) +
  labs(
    x = "2003 Day of Year (DOY)",
    y = expression(paste("Utilization Rate", "(mol m"^{-2}, "day"^{-1},")"))
  ) +
  theme_minimal()

ggsave(file.path(FIGURE_DIR, 'FigS-hail-utilization.png'),
       plot = hail_utr_utilization_plot,
       width = 6,
       height = 2.5,
       units = "in",
       dpi = 600)