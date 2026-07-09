# =============================================================================
# Shade Yield Reduction Analysis
# -----------------------------------------------------------------------------
# Purpose:
#   Simulate the effect of mid-season canopy shading on soybean yield
#   using two C allocation models ("UTR" and "Partitioning"), compare the
#   simulated yield reduction against observed shading data from
#   Proulx and Naeve (2009), and produce a bar + point figure
#   ("Fig-shade_yield_reduction.png") summarizing the comparison.
#
# Assumes the following objects already exist in the environment (created
# upstream in the analysis pipeline in 'generate_all_figures.R':
#   - initial_values, parameters, direct_modules, differential_modules, solver
#   - arg_names, optim_params_short_SoyFACE, optim_params_conversion()
# =============================================================================

# -----------------------------------------------------------------------------
# calculate_pod_reduction()
# -----------------------------------------------------------------------------
# Simulates a growing season with and without a shading event, and returns the
# proportional reduction in final pod yield caused by shading.
#
# Arguments:
#   shade_dvi  - Developmental Vegetative Index (DVI) at which shading begins
#   yr         - character year label (e.g. "2002"), used only to tag output
#   shade_pct  - fraction (0-1) by which incoming solar radiation is reduced
#                during the shading period (e.g. 0.5 = 50% shade)
#
# Returns:
#   A list with: year, utr (yield reduction, UTR model), and
#   partitioning (yield reduction, Partitioning model).
#
# Notes on method:
#   Final yield reduction is computed as
#     1 - (final Pod mass, shaded) / (final Pod mass, unshaded)
#   i.e. it assumes a fixed seed:pod (shell) ratio between the shaded and
#   unshaded runs, so no separate "shell mass" correction is needed here.
# -----------------------------------------------------------------------------
calculate_pod_reduction <- function(shade_dvi, yr, shade_pct){
  
  weather.growingseason <- weather[sd.idx:hd.ind, ]
  
  soybean$parameters$time_zone_offset <- -6
  
  ###### No-shade (control) scenario #####
  soybean_utr_optsolver_no_shade <- partial_run_biocro(initial_values,
                                                       parameters,
                                                       weather.growingseason,
                                                       direct_modules,
                                                       differential_modules,
                                                       solver,
                                                       arg_names,
                                                       verbose = FALSE)
  
  result_utr_no_shade <- soybean_utr_optsolver_no_shade(
    optim_params_conversion(optim_params_short_SoyFACE)
  )
  
  result_partitioning_no_shade <- with(soybean, run_biocro(initial_values,
                                                           parameters,
                                                           weather.growingseason,
                                                           direct_modules,
                                                           differential_modules,
                                                           solver,
                                                           verbose = FALSE))
  
  ###### Shade scenario #####
  # Find the weather-row index closest to the target shading DVI, then split
  # the growing-season weather into a "before shading" and "during/after
  # shading" segment. Solar radiation in the second segment is reduced by
  # `shade_pct` to simulate the shading treatment.
  shade.ind <- which.min(abs(result_utr_no_shade$DVI - shade_dvi)) + sd.idx
  weather.growingseason_1 <- weather[sd.idx:shade.ind, ]
  weather.growingseason_2 <- weather[shade.ind:hd.ind, ]
  weather.growingseason_2$solar <- weather.growingseason_2$solar * (1 - shade_pct)
  
  #### Before shading ####
  soybean_utr_optsolver_shade1 <- partial_run_biocro(initial_values,
                                                     parameters,
                                                     weather.growingseason_1,
                                                     direct_modules,
                                                     differential_modules,
                                                     solver,
                                                     arg_names,
                                                     verbose = FALSE)
  result_utr_shade1 <- soybean_utr_optsolver_shade1(
    optim_params_conversion(optim_params_short_SoyFACE)
  )
  
  result_partitioning_shade1 <- with(soybean, run_biocro(initial_values,
                                                         parameters,
                                                         weather.growingseason_1,
                                                         direct_modules,
                                                         differential_modules,
                                                         solver,
                                                         verbose = FALSE))
  
  #### After shading ####
  last_row <- nrow(result_partitioning_shade1)
  
  ## UTR model ##
  differential_quantities_just_before_shade <- as.list(
    result_utr_shade1[last_row, names(initial_values)]
  )
  
  soybean_utr_optsolver_shade2 <- partial_run_biocro(differential_quantities_just_before_shade,
                                                     parameters,
                                                     weather.growingseason_2,
                                                     direct_modules,
                                                     differential_modules,
                                                     solver,
                                                     arg_names,
                                                     verbose = FALSE)
  
  orignal_utr_params <- data.frame(optim_params_conversion(optim_params_short_SoyFACE))
  rownames(orignal_utr_params) <- arg_names
  colnames(orignal_utr_params) <- "Value"
  parameters_after_shade <- orignal_utr_params
  
  result_utr_shade2 <- soybean_utr_optsolver_shade2(parameters_after_shade$Value)
  
  ## Partitioning model ##
  differential_quantities_just_before_shade <- as.list(
    result_partitioning_shade1[last_row, names(soybean$initial_values)]
  )
  result_partitioning_shade2 <- run_biocro(differential_quantities_just_before_shade,
                                           soybean$parameters,
                                           weather.growingseason_2,
                                           soybean$direct_modules,
                                           soybean$differential_modules,
                                           soybean$ode_solver,
                                           verbose = FALSE)
  
  # Connect the "before" and "after" segments back into one result.
  result_utr_shade <- rbind(
    result_utr_shade1[seq_len(nrow(result_utr_shade1) - 1), ],
    result_utr_shade2
  )
  result_partitioning_shade <- rbind(
    result_partitioning_shade1[seq_len(nrow(result_partitioning_shade1) - 1), ],
    result_partitioning_shade2
  )
  
  # Final pod-yield reduction, assuming a fixed seed:pod (shell) ratio
  # between shaded and unshaded runs.
  utr_yield_reduction <- 1 - tail(result_utr_shade$Pod, 1) / tail(result_utr_no_shade$Pod, 1)
  
  # Partitioning model reports Grain and Shell separately; Pod = Grain + Shell.
  result_partitioning_shade$Pod <- result_partitioning_shade$Grain + result_partitioning_shade$Shell
  result_partitioning_no_shade$Pod <- result_partitioning_no_shade$Grain + result_partitioning_no_shade$Shell
  partitioning_yield_reduction <- 1 - tail(result_partitioning_shade$Pod, 1) /
    tail(result_partitioning_no_shade$Pod, 1)
  
  return(list(year = yr,
              utr = utr_yield_reduction,
              partitioning = partitioning_yield_reduction))
}

# -----------------------------------------------------------------------------
# Run the simulation across years and shade levels
# -----------------------------------------------------------------------------
# For each of 4 growing-season years, load that year's daily weather record,
# identify the simulation start/end indices (`sd.idx`, `hd.ind`), then run
# `calculate_pod_reduction()` for shade levels of 50/60/70/80% (applied at
# DVI = 1.5). Results are collected into `reduction_df`.
# -----------------------------------------------------------------------------
reduction_df <- list()

for (yr in c('2002', '2004', '2005', '2006')) {
  
  weather <- read.csv(file = paste0('../Data/Weather_data/', yr, '_Bondville_IL_daylength.csv'))
  
  parameters$Catm <- 372  # ambient CO2 (ppm) scenario
  
  # Experimental biomass data is only needed here to locate the simulation
  # end date (the last day of year with an observed biomass sample).
  ExpBiomass <- read.csv(file = paste0('../Data/SoyFACE_data/biomasses/', yr, '_ambient_biomass.csv'))
  colnames(ExpBiomass) <- c("DOY", "Leaf", "Stem", "Pod", "Seed", "Litter", "CumLitter")
  
  sd.idx <- which(weather$doy == 147)[12]              # planting/start day (Morgan et al. 2005)
  hd.ind <- which(weather$doy == max(ExpBiomass$DOY))[24]  # harvest/end day
  
  for (shade_pct in c(0.5, 0.6, 0.7, 0.8)) {
    shade_dvi <- 1.5
    pod_reduction <- calculate_pod_reduction(shade_dvi, yr, shade_pct)
    
    reduction_df[[length(reduction_df) + 1]] <- list(
      year       = pod_reduction$year,
      Shade_PCT  = shade_pct,
      UTR        = pod_reduction$utr * 100,           # convert to %
      Partitioning = pod_reduction$partitioning * 100  # convert to %
    )
  }
}

# -----------------------------------------------------------------------------
# Observed shading yield-reduction data
# -----------------------------------------------------------------------------
shade_exp <- data.frame(
  Shade_PCT = as.factor(c(0.5, 0.6, 0.7, 0.8)),
  Source    = "Proulx and Naeve (2009)",
  Mean_YR   = c(18.1, 28.8, 32.3, 37.9)
)

# -----------------------------------------------------------------------------
# Reshape simulation output for plotting
# -----------------------------------------------------------------------------
# Flatten the list-of-lists into a data frame (one row per year x shade level).
reduction_df_flat <- bind_rows(lapply(reduction_df, as.data.frame))
reduction_df_flat$Shade_PCT <- as.factor(reduction_df_flat$Shade_PCT)

# Long format: one row per (year, shade level, model), used to draw the
# individual per-year points on top of the bars.
df_points <- reduction_df_flat %>%
  tidyr::pivot_longer(cols = c("UTR", "Partitioning"),
                      names_to = "Source", values_to = "Yield_R")

# Mean across years, per shade level and model, used to draw the bars.
df_bars <- df_points %>%
  group_by(Shade_PCT, Source) %>%
  summarise(Mean_YR = mean(Yield_R), .groups = "drop")

# Combine simulated bar means with the observed bar values.
all_shade_bars <- bind_rows(df_bars, shade_exp)

# Fix factor level order so bars/points are colored and legended consistently:
# observed data first, then the two simulated models.
source_levels <- c("Proulx and Naeve (2009)", "UTR", "Partitioning")
all_shade_bars$Source <- factor(all_shade_bars$Source, levels = source_levels)
df_points$Source      <- factor(df_points$Source, levels = source_levels)

# -----------------------------------------------------------------------------
# Manual dodge positions
# -----------------------------------------------------------------------------
dodge_offset <- data.frame(
  Source = c("Proulx and Naeve (2009)", "UTR", "Partitioning"),
  offset = c(-0.267, 0, 0.267)
)

shade_levels <- c(0.5, 0.6, 0.7, 0.8)  # order used for the x-axis

all_shade_bars_pos <- all_shade_bars %>%
  left_join(dodge_offset, by = "Source") %>%
  mutate(x_idx = match(as.numeric(as.character(Shade_PCT)), shade_levels),
         x_pos = x_idx + offset)

df_points_pos <- df_points %>%
  left_join(dodge_offset, by = "Source") %>%
  mutate(x_idx = match(as.numeric(as.character(Shade_PCT)), shade_levels),
         x_pos = x_idx + offset)

# -----------------------------------------------------------------------------
# Build and save the figure
# -----------------------------------------------------------------------------
# Bars show the mean yield reduction (observed, UTR, Partitioning) for each
# shade level; open circles show individual per-year simulated values.
shade_reduction_fig <- ggplot() +
  geom_bar(data = all_shade_bars_pos,
           aes(x = x_pos, y = Mean_YR, fill = Source),
           stat = "identity", width = 0.25) +
  geom_jitter(data = df_points_pos,
              aes(x = x_pos, y = Yield_R, fill = Source),
              width = 0.04, height = 0,
              shape = 21, size = 2.5, color = "black") +
  scale_fill_manual(name = "Source",
                    values = c("Proulx and Naeve (2009)" =  "#009E73",
                               "UTR"                      = "#E65F00",
                               "Partitioning"             = "#0072B2"),
                    breaks = source_levels) +
  scale_x_continuous(breaks = 1:4, labels = c("50%", "60%", "70%", "80%")) +
  labs(x = "Shade (%)",
       y = "Yield Reduction (%)") +
  theme_bw() +
  theme(legend.position = "bottom")

print(shade_reduction_fig)

ggsave(file.path(FIGURE_DIR, 'Fig-shade_yield_reduction.png'),
       plot = shade_reduction_fig,
       width = 6,
       height = 4,
       units = "in",
       dpi = 600)