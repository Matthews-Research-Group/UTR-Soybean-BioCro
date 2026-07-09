library(ggtext)

# =============================================================================
# PURPOSE OF THIS SCRIPT
# -----------------------------------------------------------------------------
# For a set of historical growing seasons, simulate soybean yield loss
# caused by defoliation (leaf removal) at two crop stages
# (R4, DVI = 1.35 and R5, DVI = 1.5) and a range of defoliation intensities
# (0%, 25%, 50%, 75%, ~100%).
#
# Two C allocation crop models are compared, the UTR and the partitioning 
# models, and their simulated yield reductions are compared to field
# experimental results from Parvej et al. (2025).
#
# Each simulated defoliation event is run in two stages:
#   1) run the model normally up to the time of defoliation,
#   2) remove `defoliate_frac` of the leaf biomass from the model state, 
#   then resume the simulation for the rest of the season.
#
# The resulting simulated yield reductions (relative to a no-defoliation
# simulation) are then plotted alongside yield reduction data from 
# Parvej et al. (2025).
# =============================================================================


# -----------------------------------------------------------------------------
# update_differential_quantities()
#
# After a defoliation event, this function takes the model's differential
# (state) variables at the last simulated time step and reduces the leaf
# biomass or carbon pools by `defoliate_frac` (the fraction of leaf removed),
# leaving every other state variable unchanged. This updated state is then used 
# as the starting point ("initial_values") for the second part of the
# simulation (i.e., after defoliation).
#
# Arguments:
#   r              - data frame of simulated results up to defoliation; the
#                   last row holds the state just before defoliation.
#   updated_values - names of the state variables to specify
#   model          - which crop model is being used ("utr" or "partitioning"),
#                   since the two models name their carbon pools differently.
#   defoliate_frac - fraction of leaf biomass removed (e.g., 0.5 = 50%).
# -----------------------------------------------------------------------------
update_differential_quantities <- 
  function(r, updated_values, model, defoliate_frac){
  # Fraction of leaf biomass that remains after defoliation
  leaf_remaining_percentage <- 1 - defoliate_frac
  
  last_row <- nrow(r)
  differential_quantities_just_before_defoliation <-
    as.list(r[last_row, updated_values])
  
  differential_quantities_just_after_defoliation <-
    differential_quantities_just_before_defoliation
  
  if (model == 'utr'){
    # UTR model tracks leaf carbon in two separate pools: substrate and
    # structural carbon. Both are reduced by the same defoliation fraction.
    differential_quantities_just_after_defoliation$Leaf_substrate_carbon <-
      differential_quantities_just_before_defoliation$Leaf_substrate_carbon *
      leaf_remaining_percentage
    differential_quantities_just_after_defoliation$Leaf_structural_carbon <-
      differential_quantities_just_before_defoliation$Leaf_structural_carbon *
      leaf_remaining_percentage
  }
  else if (model == 'partitioning'){
    # Partitioning model tracks leaf biomass as a single pool.
    differential_quantities_just_after_defoliation$Leaf <-
      differential_quantities_just_before_defoliation$Leaf *
      leaf_remaining_percentage
  }
  
  return(differential_quantities_just_after_defoliation)
}


# -----------------------------------------------------------------------------
# calculate_pod_reduction()
#
# Runs the full "no defoliation vs. defoliation" comparison for a single
# year, a single defoliation timing (`defoliation_dvi`), and a single
# defoliation intensity (`defoliate_frac`), for both the UTR and partitioning
# models, and returns the resulting fractional yield (pod biomass) reduction
# for each model.
#
# Arguments:
#   defoliation_dvi - developmental stage (DVI) at which defoliation occurs
#                      (e.g., 1.35 for R4, 1.5 for R5).
#   yr               - growing season year
#   defoliate_frac   - fraction of leaf biomass removed at defoliation.
#
# Returns a list with:
#   year         - the growing season year
#   utr          - fractional yield reduction from the UTR model
#   partitioning - fractional yield reduction from the partitioning model
# -----------------------------------------------------------------------------
calculate_pod_reduction <- function(defoliation_dvi, yr, defoliate_frac){
  weather.growingseason <- weather[sd.idx: hd.ind,]
  soybean$parameters$time_zone_offset <- -6
  
  ###### No Defoliation scenario (control run, full season) #####
  soybean_utr_optsolver_no_defoliation <- partial_run_biocro(initial_values,
                                                             parameters,
                                                             weather.growingseason,
                                                             direct_modules,
                                                             differential_modules,
                                                             solver,
                                                             arg_names,
                                                             verbose = FALSE)
  
  result_utr_no_defoliation <- soybean_utr_optsolver_no_defoliation(optim_params_conversion(optim_params_short_SoyFACE))
  result_partitioning_no_defoliation <- with(soybean, run_biocro(initial_values,
                                                                 parameters,
                                                                 weather.growingseason,
                                                                 direct_modules,
                                                                 differential_modules,
                                                                 solver,
                                                                 verbose = FALSE))
  
  ###### Defoliation scenario #####
  # Find the time index in the no-defoliation run whose DVI is closest to
  # the requested defoliation stage; this splits the weather data into
  # "before defoliation" and "after defoliation" parts
  defoliate.ind <- which.min(abs(result_utr_no_defoliation$DVI - defoliation_dvi)) + sd.idx
  weather.growingseason_1 <- weather[sd.idx:defoliate.ind,]
  weather.growingseason_2 <- weather[defoliate.ind:hd.ind,]
  
  #### Stage 1: run normally up to the moment of defoliation ####
  soybean_utr_optsolver_defoliation1 <- partial_run_biocro(initial_values,
                                                           parameters,
                                                           weather.growingseason_1,
                                                           direct_modules,
                                                           differential_modules,
                                                           solver,
                                                           arg_names,
                                                           verbose = FALSE)
  result_utr_defoliation1 <- soybean_utr_optsolver_defoliation1(optim_params_conversion(optim_params_short_SoyFACE))
  result_partitioning_defoliation1 <- with(soybean, run_biocro(initial_values,
                                                               parameters,
                                                               weather.growingseason_1,
                                                               direct_modules,
                                                               differential_modules,
                                                               solver,
                                                               verbose = FALSE))
  
  #### Stage 2: remove leaf biomass, then resume for the rest of the season ####
  
  # -- UTR model --
  differential_quantities_just_after_defoliation_utr <- update_differential_quantities(
    result_utr_defoliation1,
    names(initial_values),
    'utr',
    defoliate_frac)
  
  soybean_utr_optsolver_defoliation2 <- partial_run_biocro(differential_quantities_just_after_defoliation_utr,
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
  parameters_after_defoliation <- orignal_utr_params
  
  result_utr_defoliation2 <- soybean_utr_optsolver_defoliation2(parameters_after_defoliation$Value)
  
  # -- Partitioning model --
  differential_quantities_just_after_defoliation_partitioning <- update_differential_quantities(
    result_partitioning_defoliation1,
    names(soybean$initial_values),
    'partitioning',
    defoliate_frac)
  
  result_partitioning_defoliation2 <- run_biocro(differential_quantities_just_after_defoliation_partitioning,
                                                 soybean$parameters,
                                                 weather.growingseason_2,
                                                 soybean$direct_modules,
                                                 soybean$differential_modules,
                                                 soybean$ode_solver,
                                                 verbose = FALSE)
  
  # Connect the "before" and "after" defoliation runs back together into a
  # single full growing season result for each model.
  result_utr_defoliation <- rbind(result_utr_defoliation1
                                  [seq_len(nrow(result_utr_defoliation1) - 1), ],
                                  result_utr_defoliation2)
  result_partitioning_defoliation <- rbind(result_partitioning_defoliation1
                                           [seq_len(nrow(result_partitioning_defoliation1) - 1), ],
                                           result_partitioning_defoliation2)
  
  # Fractional yield reduction
  utr_yield_reduction <- 1 - tail(result_utr_defoliation$Pod, 1) / 
    tail(result_utr_no_defoliation$Pod, 1)
  
  # Partitioning model reports pod biomass as separate Grain + Shell pools,
  # so these are summed to get total pod biomass before comparing.
  result_partitioning_defoliation$Pod <- result_partitioning_defoliation$Grain + 
    result_partitioning_defoliation$Shell
  result_partitioning_no_defoliation$Pod <- result_partitioning_no_defoliation$Grain + 
    result_partitioning_no_defoliation$Shell
  partitioning_yield_reduction <- 1 - tail(result_partitioning_defoliation$Pod, 1) / 
    tail(result_partitioning_no_defoliation$Pod, 1)
  
  return(list(year = yr,
              utr = utr_yield_reduction,
              partitioning = partitioning_yield_reduction))
}


# =============================================================================
# RUN THE SIMULATION
# -----------------------------------------------------------------------------
# Loop over 4 growing seasons, 2 defoliation stages (R4/R5), and
# 5 defoliation intensities, collecting the yield-reduction percentages.
# =============================================================================
reduction_df <- list()
defoliate_fracs <- c(0, 0.25, 0.5, 0.75, 0.999)
defoliate_DVIs <- c(1.35, 1.5)

for (yr in c('2002', '2004', '2005', '2006')){
  # Load the daily weather record for this growing season.
  weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))
  co2_opt <- '_ambient_'
  parameters$Catm <- 372
  
  # Load observed biomass data, used only to define the start/end of the
  # simulated growing season window (planting -> harvest DOY).
  ExpBiomass <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt , 'biomass.csv'))
  colnames(ExpBiomass)<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  # sd.idx / hd.ind mark the simulated planting (sowing) and harvest dates,
  # matched to the corresponding day-of-year in the weather file.
  sd.idx <- which(weather$doy == 147)[12] # Morgan et al. 2005
  hd.ind <- which(weather$doy == max(ExpBiomass$DOY))[24]
  
  for(defoliate_dvi in defoliate_DVIs){
    for(defoliate_frac in defoliate_fracs){
      pod_reduction <- calculate_pod_reduction(defoliate_dvi, yr, defoliate_frac)
      
      # Store just the fields needed for the final comparison plot.
      reduction_df[[length(reduction_df) + 1]] <- list(
        year             = pod_reduction$year,
        Defoliation_DVI  = defoliate_dvi,
        Defoliation_frac  = defoliate_frac,
        UTR              = pod_reduction$utr * 100,
        Partitioning     = pod_reduction$partitioning * 100
      )
    }
  }
}

# Flatten the `reduction_df` list into a single data frame.
df <- bind_rows(lapply(reduction_df, function(r) {
  if (any(sapply(r, length) == 0)) return(NULL)
  as.data.frame(r)
}))

df$Defoliation_DVI <- as.factor(df$Defoliation_DVI)


# =============================================================================
# OBSERVED FIELD DATA FOR COMPARISON
# -----------------------------------------------------------------------------
# Yield-reduction values from Parvej et al. (2025), for defoliation
# applied at R4 (DVI = 1.35) and R5 (DVI = 1.5) in Iowa and Indiana 
# field trials.
# =============================================================================
exp_yield_reduction <- data.frame(
  Defoliation_DVI = c(
    rep(1.35, 5), rep(1.5, 5),   # Iowa
    rep(1.35, 5), rep(1.5, 5)    # Indiana
  ),
  Defoliation_frac = rep(c(0, 25, 50, 75, 99.9), 4),
  Yield_R = 100-c(
    # Iowa - R4 (DVI 1.35)
    95.1, 93.0, 82.0, 62.1, 33.3,
    # Iowa - R5 (DVI 1.5)
    95.1, 88.8, 73.6, 49.5, 16.4,
    # Indiana - R4 (DVI 1.35)
    95.7, 93.0, 79.5, 55.2, 20.3,
    # Indiana - R5 (DVI 1.5)
    95.7, 93.0, 79.5, 55.2, 20.3
  ),
  Source = "Parvej et al. (2025)"
)


# =============================================================================
# PREPARE DATA FOR PLOTTING
# =============================================================================

# Reshape simulated results to long format: one row per (year, DVI stage,
# defoliation %, model) combination, with a single Yield_R column.
# Defoliation_frac is rescaled from a 0-1 fraction to a 0-100 percentage to
# match the field data.
df_long <- df %>%
  mutate(Defoliation_frac = Defoliation_frac * 100) %>%
  pivot_longer(cols = c("UTR", "Partitioning"),
               names_to = "Source",
               values_to = "Yield_R")

# Average across years to get one bar per (DVI stage, defoliation %,
# model).
df_bars_models <- df_long %>%
  group_by(Defoliation_DVI, Defoliation_frac, Source) %>%
  summarise(Mean_YR = mean(Yield_R), .groups = "drop")

# Average the field data across locations (Iowa/Indiana) to get one mean
# per DVI stage - defoliation % scenario.
exp_data_formatted <- exp_yield_reduction %>%
  mutate(Defoliation_DVI = as.factor(Defoliation_DVI))

df_bars_parvej <- exp_data_formatted %>%
  group_by(Defoliation_DVI, Defoliation_frac, Source) %>%
  summarise(Mean_YR = mean(Yield_R), .groups = "drop")

# Combine simulated and observed data: bars show the mean, 
# points show the individual underlying observations/replicates.
all_bars <- bind_rows(df_bars_models, df_bars_parvej)
all_points <- bind_rows(
  df_long %>% select(Defoliation_DVI, Defoliation_frac, Source, Yield_R),
  exp_data_formatted %>% select(Defoliation_DVI, Defoliation_frac, Source, Yield_R)
)

# Fix factor level order so bars/points are consistently colored and ordered
# (field data first, then the two simulated models).
source_levels <- c("Parvej et al. (2025)", "UTR", "Partitioning")
all_bars$Source <- factor(all_bars$Source, levels = source_levels)
all_points$Source <- factor(all_points$Source, levels = source_levels)


# =============================================================================
# BUILD AND SAVE THE FINAL FIGURE
# =============================================================================
facet_labels <- c(
  "1.35" = "(A)                                                             R4 (DVI = 1.35)",
  "1.5"  = "(B)                                                             R5 (DVI = 1.5)"
)

yield_reduction_fig <- ggplot() +
  geom_bar(data = all_bars,
           aes(x = as.factor(Defoliation_frac), y = Mean_YR, fill = Source),
           stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(data = all_points,
             aes(x = as.factor(Defoliation_frac), y = Yield_R,
                 group = Source, fill = Source),
             position = position_dodge(width = 0.8),
             shape = 21, size = 1.5, color = "black", alpha = 0.7) +
  facet_wrap(~Defoliation_DVI, nrow = 2,
             labeller = as_labeller(facet_labels)) +
  scale_fill_manual(values = c("Parvej et al. (2025)" = "#009E73",
                               "UTR" =                  "#E65F00",
                               "Partitioning" =         "#0072B2")) +
  labs(x = "Defoliation Percentage (%)",
       y = "Yield Reduction (%)",
       fill = "Source") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(hjust = 0))

print(yield_reduction_fig)

ggsave(file.path(FIGURE_DIR, 'Fig-defoliation_comparison.png'),
       plot = yield_reduction_fig,
       width = 6, height = 6, dpi = 600)