# =============================================================================
# Pod Removal Yield Reduction Analysis
# -----------------------------------------------------------------------------
# Purpose:
#   Simulate the effect of pod removal at R5 on soybean final yield using
#   two C allocation models ("UTR" and "Partitioning"), compare the
#   simulated yield reduction against observed yield reduction data
#   from Proulx and Naeve (2009).
#
# Assumes the following objects already exist in the environment (created
# upstream in the analysis pipeline in 'generate_all_figures.R':
#   - initial_values, parameters, direct_modules, differential_modules, solver
#   - arg_names, optim_params_short_SoyFACE, optim_params_conversion()
# =============================================================================

# -----------------------------------------------------------------------------
# update_differential_quantities()
# -----------------------------------------------------------------------------
# Reduces the pod-related differential-state variable(s) by the fraction of
# pods removed.
# Used to set the model's starting state for the "after pod removal"
# simulation segment.
#
# Arguments:
#   r                - data frame of simulation output up to (and including)
#                      the time of pod removal
#   updated_values   - names of the differential-state variables to be specified
#   UPDATE_PARAMETERS - if TRUE, also step DVI back by 0.2 to represent a
#                      developmental reset caused by pod removal. Currently
#                      always called with FALSE (no DVI adjustment), but the
#                      option is kept here since the calling code toggles it.
#   model            - which model's state layout to use: 'utr' or
#                      'partitioning'
#   pod_removal_frac  - fraction (0-1) of pod mass removed
#
# Returns:
#   A named list of state values to use as the starting point for the
#   post-pod-removal simulation segment.
# -----------------------------------------------------------------------------
update_differential_quantities <- 
  function(r, updated_values, UPDATE_PARAMETERS, model, pod_removal_frac){
  pod_remaining_percentage <- 1 - pod_removal_frac
  
  last_row <- nrow(r)
  differential_quantities_just_before_podrm <-
    as.list(r[last_row, updated_values])
  
  differential_quantities_just_after_podrm <-
    differential_quantities_just_before_podrm
  
  if (model == 'utr'){
    differential_quantities_just_after_podrm$Pod_structural_carbon <-
      differential_quantities_just_before_podrm$Pod_structural_carbon * pod_remaining_percentage
  }
  else if (model == 'partitioning'){
    differential_quantities_just_after_podrm$Shell <-
      differential_quantities_just_before_podrm$Shell * pod_remaining_percentage
    
    differential_quantities_just_after_podrm$Grain <-
      differential_quantities_just_before_podrm$Grain * pod_remaining_percentage
  }
  
  if (UPDATE_PARAMETERS){
    differential_quantities_just_after_podrm$DVI <-
      differential_quantities_just_before_podrm$DVI - 0.2
  }
  
  return(differential_quantities_just_after_podrm)
}

# -----------------------------------------------------------------------------
# calculate_pod_reduction()
# -----------------------------------------------------------------------------
# Simulates a growing season with and without pod-removal at R5,
# and returns the reduction in final pod yield caused by pod removal,
# as estimated by both the "UTR" and "Partitioning" models.
#
# Arguments:
#   podrm_dvi        - Developmental Index (DVI) at which pod removal occurs
#   yr               - year of the simulation
#   pod_removal_frac - fraction (0-1) of pod mass removed
#
# Returns:
#   A list with: year, utr (yield reduction predicted by the UTR model), and
#   partitioning (yield reduction predicted by the partitioning model).
#
# Notes on method:
#   Final yield reduction is computed as
#     1 - (final pod mass, pod-removed) / (final pod mass,  control)
#   i.e. it assumes a fixed seed:pod (shell) ratio between the pod-removed
#   and control runs.
# -----------------------------------------------------------------------------
calculate_pod_reduction <- function(podrm_dvi, yr, pod_removal_frac){
  
  weather.growingseason <- weather[sd.idx: hd.ind,]
  soybean$parameters$time_zone_offset <- -6
  
  ###### No-pod-removal (control) scenario #####
  soybean_utr_optsolver_no_podrm <- partial_run_biocro(initial_values,
                                                       parameters,
                                                       weather.growingseason,
                                                       direct_modules,
                                                       differential_modules,
                                                       solver,
                                                       arg_names,
                                                       verbose = FALSE)
  
  result_utr_no_podrm <- soybean_utr_optsolver_no_podrm(optim_params_conversion(optim_params_short_SoyFACE))
  result_partitioning_no_podrm <- with(soybean, run_biocro(initial_values,
                                                           parameters,
                                                           weather.growingseason,
                                                           direct_modules,
                                                           differential_modules,
                                                           solver,
                                                           verbose = FALSE))
  
  ###### Pod removal scenario #####
  # Find the weather-row index closest to the target pod-removal DVI, then
  # split the growing-season weather into "before pod removal" and
  # "during/after pod removal" parts
  pod.rm.ind <- which.min(abs(result_utr_no_podrm$DVI - podrm_dvi)) + sd.idx
  weather.growingseason_1 <- weather[sd.idx:pod.rm.ind,]
  weather.growingseason_2 <- weather[pod.rm.ind:hd.ind,]
  
  #### Before pod removal ####
  soybean_utr_optsolver_podrm1 <- partial_run_biocro(initial_values,
                                                     parameters,
                                                     weather.growingseason_1,
                                                     direct_modules,
                                                     differential_modules,
                                                     solver,
                                                     arg_names,
                                                     verbose = FALSE)
  result_utr_podrm1 <- soybean_utr_optsolver_podrm1(optim_params_conversion(optim_params_short_SoyFACE))
  result_partitioning_podrm1 <- with(soybean, run_biocro(initial_values,
                                                         parameters,
                                                         weather.growingseason_1,
                                                         direct_modules,
                                                         differential_modules,
                                                         solver,
                                                         verbose = FALSE))
  
  #### After pod removal ####
  # Reduce the pod-related state variables by the removed fraction,
  # then carry the adjusted state forward as the starting point for 
  # the "after pod removal" run.
  
  ## UTR model ##
  UPDATE_PARAMETERS <- FALSE
  differential_quantities_just_after_podrm_utr <- 
    update_differential_quantities(result_utr_podrm1,
                                   names(initial_values),
                                   UPDATE_PARAMETERS, 'utr',
                                   pod_removal_frac)
  
  soybean_utr_optsolver_podrm2 <- 
    partial_run_biocro(differential_quantities_just_after_podrm_utr,
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
  parameters_after_podrm <- orignal_utr_params
  
  result_utr_podrm2 <- soybean_utr_optsolver_podrm2(
    parameters_after_podrm$Value)
  
  ## Partitioning model ##
  differential_quantities_just_after_podrm_partitioning <- 
    update_differential_quantities(result_partitioning_podrm1,
                                   names(soybean$initial_values),
                                   UPDATE_PARAMETERS, 'partitioning',
                                   pod_removal_frac)
  
  result_partitioning_podrm2 <- run_biocro(differential_quantities_just_after_podrm_partitioning,
                                           soybean$parameters,
                                           weather.growingseason_2,
                                           soybean$direct_modules,
                                           soybean$differential_modules,
                                           soybean$ode_solver,
                                           verbose = FALSE)
  
  # Connect the "before" and "after" segments back into one data frame
  result_utr_podrm <- rbind(result_utr_podrm1[seq_len(nrow(result_utr_podrm1) - 1), ], 
                            result_utr_podrm2)
  result_partitioning_podrm <- rbind(result_partitioning_podrm1[seq_len(nrow(result_partitioning_podrm1) - 1), ], 
                                     result_partitioning_podrm2)
  
  # Final pod-yield reduction, assuming a fixed seed:pod (shell) ratio
  # between the pod-removed and control runs.
  utr_yield_reduction <- 1 - tail(result_utr_podrm$Pod, 1) / tail(result_utr_no_podrm$Pod, 1)
  
  # Partitioning model reports Grain and Shell separately; Pod = Grain + Shell.
  result_partitioning_podrm$Pod <- result_partitioning_podrm$Grain + 
    result_partitioning_podrm$Shell
  result_partitioning_no_podrm$Pod <- result_partitioning_no_podrm$Grain + 
    result_partitioning_no_podrm$Shell
  partitioning_yield_reduction <- 1 - tail(result_partitioning_podrm$Pod, 1) / 
    tail(result_partitioning_no_podrm$Pod, 1)
  
  return(list(year = yr,
              utr = utr_yield_reduction,
              partitioning = partitioning_yield_reduction))
}

# -----------------------------------------------------------------------------
# Run the simulation across years and pod-removal levels
# -----------------------------------------------------------------------------
# For each of 4 growing-season years, load that year's daily weather record,
# identify the simulation start/end indices (`sd.idx`, `hd.ind`), then run
# `calculate_pod_reduction()` for pod-removal levels of 20/40/60/70% (applied
# at DVI = 1.6). Results are collected into `reduction_df`.
# -----------------------------------------------------------------------------
reduction_df <- list()

for (yr in c('2002', '2004', '2005', '2006')){
  
  weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))
  parameters$Catm <- 372  # ambient CO2 (ppm) scenario
  
  # Experimental biomass data is only needed here to locate the simulation
  # end date (the last day of year with an observed biomass sample).
  ExpBiomass <- read.csv(file = paste0('../Data/SoyFACE_data/biomasses/', yr, '_ambient_biomass.csv'))
  colnames(ExpBiomass) <- c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  sd.idx <- which(weather$doy == 147)[12]              # planting/start day (Morgan et al. 2005)
  hd.ind <- which(weather$doy == max(ExpBiomass$DOY))[24]  # harvest/end day
  
  for(pod_removal_frac in c(0.2, 0.4, 0.6, 0.7)){
    pod_removal_dvi <- 1.6
    pod_reduction <- calculate_pod_reduction(pod_removal_dvi, yr, pod_removal_frac)
    
    reduction_df[[length(reduction_df) + 1]] <- list(
      year         = pod_reduction$year,
      Pod_RM       = pod_removal_frac,
      UTR          = pod_reduction$utr * 100,          # convert to %
      Partitioning = pod_reduction$partitioning * 100  # convert to %
    )
  }
}

# Reshape simulation output for plotting
# Flatten the list-of-lists into a data frame (one row per year x pod-removal
# level).
df <- bind_rows(lapply(reduction_df, as.data.frame))
df$Pod_RM <- as.factor(df$Pod_RM)

# Observed pod-removal yield-reduction data
obs_df <- data.frame(
  Pod_RM  = as.factor(c(0.2, 0.4, 0.6, 0.7)),
  Source  = "Proulx and Naeve (2009)",
  Mean_YR = c(11.7, 33.9, 50.1, 52.3)
)

# Long format: one row per (year, pod-removal level, model), used to draw the
# individual per-year points on top of the bars for the means.
df_points <- df %>%
  tidyr::pivot_longer(cols = c("UTR", "Partitioning"),
                      names_to = "Source", values_to = "Yield_R")

# Mean across years, per pod-removal level and model, used to draw the bars.
df_bars <- df_points %>%
  group_by(Pod_RM, Source) %>%
  summarise(Mean_YR = mean(Yield_R), .groups = "drop")

# Combine simulated means with the observed means.
all_depod_bars <- bind_rows(df_bars, obs_df)

# Fix factor level order so bars/points have consistent legends
source_levels <- c("Proulx and Naeve (2009)", "UTR", "Partitioning")
all_depod_bars$Source <- factor(all_depod_bars$Source, levels = source_levels)
df_points$Source <- factor(df_points$Source, levels = source_levels)


# Manual dodge positions
dodge_offset <- data.frame(
  Source = c("Proulx and Naeve (2009)", "UTR", "Partitioning"),
  offset = c(-0.267, 0, 0.267)
)

pod_levels <- c(0.2, 0.4, 0.6, 0.7)  # order used for the x-axis

all_depod_bars_pos <- all_depod_bars %>%
  left_join(dodge_offset, by = "Source") %>%
  mutate(x_idx = match(as.numeric(as.character(Pod_RM)), pod_levels),
         x_pos = x_idx + offset)

df_points_pos <- df_points %>%
  left_join(dodge_offset, by = "Source") %>%
  mutate(x_idx = match(as.numeric(as.character(Pod_RM)), pod_levels),
         x_pos = x_idx + offset)

# Build and save the figure
yield_reduction_fig <- ggplot() +
  geom_bar(data = all_depod_bars_pos,
           aes(x = x_pos, y = Mean_YR, fill = Source),
           stat = "identity", width = 0.25) +
  geom_jitter(data = df_points_pos,
              aes(x = x_pos, y = Yield_R, fill = Source),
              width = 0.04, height = 0,
              shape = 21, size = 2.5, color = "black") +
  scale_fill_manual(name = "Source",
                    values = c("Proulx and Naeve (2009)" = "#009E73",
                               "UTR"                     = "#E65F00",
                               "Partitioning"            = "#0072B2"),
                    breaks = source_levels) +
  scale_x_continuous(breaks = 1:4, labels = c("20%", "40%", "60%", "70%")) +
  labs(x = "Pod Removal (%)",
       y = "Yield Reduction (%)") +
  theme_bw() +
  theme(legend.position = "bottom")

print(yield_reduction_fig)

ggsave(file.path(FIGURE_DIR, 'Fig-pod_removal_yield_reduction.png'),
       plot = yield_reduction_fig,
       width = 6,
       height = 4,
       units = "in",
       dpi = 600
)