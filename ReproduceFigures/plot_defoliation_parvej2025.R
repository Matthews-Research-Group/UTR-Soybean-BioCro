# define a function to update differential values
update_differential_quantities <- function(r, updated_values, UPDATE_PARAMETERS, model, defoliate_ptc){
  ## Define remaining percentages
  leaf_remaining_percentage <- 1 - defoliate_ptc
  
  last_row <- nrow(r)
  differential_quantities_just_before_defoliation <-
    as.list(r[last_row, updated_values])
  
  differential_quantities_just_after_defoliation <-
    differential_quantities_just_before_defoliation
  
  if (model == 'utr'){
    differential_quantities_just_after_defoliation$Leaf_substrate_carbon <-
      differential_quantities_just_before_defoliation$Leaf_substrate_carbon * 
      leaf_remaining_percentage
    differential_quantities_just_after_defoliation$Leaf_structural_carbon <-
      differential_quantities_just_before_defoliation$Leaf_structural_carbon * 
      leaf_remaining_percentage # Could be changed to a different percentage
  } 
  else if (model == 'partitioning'){
    differential_quantities_just_after_defoliation$Leaf <-
      differential_quantities_just_before_defoliation$Leaf * 
      leaf_remaining_percentage
  }
  
  return(differential_quantities_just_after_defoliation)
}

# define a function to calculate yield reduction from pod removal
calculate_pod_reduction <- function(defoliation_dvi, yr, defoliate_pct){
  weather.growingseason <- weather[sd.idx: hd.ind,]
  RootVals <- data.frame("DOY"=ExpBiomass$DOY[3], "Root"=0.17*sum(ExpBiomass[5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  soybean$parameters$time_zone_offset <- -6
  ###### No Defoliation scenario #####
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
  defoliate.ind <- which.min(abs(result_utr_no_defoliation$DVI - defoliation_dvi)) + sd.idx
  weather.growingseason_1 <- weather[sd.idx:defoliate.ind,]
  weather.growingseason_2 <- weather[defoliate.ind:hd.ind,]
  
  #### Before Defoliation #### 
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
  #### After Defoliation #### 
  ### Without Parameter Update ####
  ## UTR ##
  UPDATE_PARAMETERS <- FALSE
  differential_quantities_just_after_defoliation_utr <- update_differential_quantities(result_utr_defoliation1, 
                                                                                 names(initial_values), 
                                                                                 UPDATE_PARAMETERS, 'utr',
                                                                                 defoliate_pct)
  
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
  
  ## Partitioning model##
  differential_quantities_just_after_defoliation_partitioning <- update_differential_quantities(result_partitioning_defoliation1, 
                                                                                          names(soybean$initial_values),
                                                                                          UPDATE_PARAMETERS, 'partitioning',
                                                                                          defoliate_pct)
  
  result_partitioning_defoliation2 <- run_biocro(differential_quantities_just_after_defoliation_partitioning,
                                           soybean$parameters,
                                           weather.growingseason_2,
                                           soybean$direct_modules,
                                           soybean$differential_modules,
                                           soybean$ode_solver,
                                           verbose = FALSE)
  
  result_utr_defoliation <- rbind(result_utr_defoliation1
                                  [seq_len(nrow(result_utr_defoliation1) - 1), ], 
                                  result_utr_defoliation2)
  result_partitioning_defoliation <- rbind(result_partitioning_defoliation1
                                           [seq_len(nrow(result_partitioning_defoliation1) - 1), ], 
                                           result_partitioning_defoliation2)
  
  # Method 3: Assume seed:pod ratio to be fixed
  utr_yield_reduction3 <- 1 - tail(result_utr_defoliation$Pod, 1)/tail(result_utr_no_defoliation$Pod, 1)
  
  # Partitioning Model Yield Change
  result_partitioning_defoliation$Pod <- result_partitioning_defoliation$Grain  + result_partitioning_defoliation$Shell
  result_partitioning_no_defoliation$Pod <- result_partitioning_no_defoliation$Grain  + result_partitioning_no_defoliation$Shell
  partitioning_yield_reduction <- 1 - tail(result_partitioning_defoliation$Pod, 1) / tail(result_partitioning_no_defoliation$Pod, 1)
  
  return(list(year=yr, 
              utr=utr_yield_reduction3, 
              partitioning=partitioning_yield_reduction))
}


reduction_df <- list()
defolicate_pcts <- c(0, 0.25, 0.5, 0.75, 0.999)
defolicate_DVIs <- c(1.4, 1.5)

for (yr in c('2002', '2004', '2005', '2006')){
  weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength.csv'))
  co2_opt <- '_ambient_'
  parameters$Catm <- 372
  ExpBiomass <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt , 'biomass.csv'))
  colnames(ExpBiomass)<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  ExpBiomass.std <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, co2_opt, 'biomass_std.csv'))
  colnames(ExpBiomass.std)<-c("DOY","Leaf","Stem","Pod", "Seed", "Litter", "CumLitter")
  
  sd.idx <- which(weather$doy == 147)[12] # Morgan et al. 2005
  hd.ind <- which(weather$doy == max(ExpBiomass$DOY))[24]
  
  for(defoliate_dvi in defolicate_DVIs){
    for(defoliate_pct in defolicate_pcts){
      pod_reduction <- calculate_pod_reduction(defoliate_dvi, yr, defoliate_pct)
      print(pod_reduction)
      # Append each result to the list
      reduction_df[[length(reduction_df) + 1]] <- list(
        year             = pod_reduction$year,
        Defolication_DVI = defoliate_dvi,
        Defolication_PCT = defoliate_pct,
        UTR              = pod_reduction$utr * 100,
        Partitioning     = pod_reduction$partitioning * 100
      )
    }
  } 
}

# Flatten reduction_df
df <- bind_rows(lapply(reduction_df, function(r) {
  if (any(sapply(r, length) == 0)) return(NULL)
  as.data.frame(r)
}))
df$Defolication_DVI <- as.factor(df$Defolication_DVI)

# Creating the digitized dataframe from Parvej et al. (2025)
exp_yield_reduction <- data.frame(
  location = c(
    rep("Iowa", 10),    # 5 points for R4, 5 points for R5
    rep("Indiana", 10)  # 5 points for R4, 5 points for R5
  ),
  Defoliation_DVI = c(
    # Iowa
    rep(1.4, 5), rep(1.5, 5),
    # Indiana
    rep(1.4, 5), rep(1.5, 5)
  ),
  Defoliation_PCT = rep(c(0, 25, 50, 75, 99.9), 4),
  Mean_YR = 100-c(
    # Iowa - R4 (DVI 1.4)
    95.1, 93.0, 82.0, 62.1, 33.3,
    # Iowa - R5 (DVI 1.5)
    95.1, 88.8, 73.6, 49.5, 16.4,
    # Indiana - R4 (DVI 1.4)
    95.7, 93.0, 79.5, 55.2, 20.3,
    # Indiana - R5 (DVI 1.5) 
    95.7, 93.0, 79.5, 55.2, 20.3
  ),
  Source = "Parvej et al. (2025)"
)

# --- 1. Data Preparation ---

# Process your model data (df)
# Ensure Defolication_PCT is on the same scale as exp_yield_reduction (0-100)
df_long <- df %>%
  mutate(Defolication_PCT = Defolication_PCT * 100) %>%
  pivot_longer(cols = c("UTR", "Partitioning"), 
               names_to = "Source", 
               values_to = "Yield_R")

# Calculate means for the bars (Models)
df_bars_models <- df_long %>%
  group_by(Defolication_DVI, Defolication_PCT, Source) %>%
  summarise(Mean_YR = mean(Yield_R), .groups = "drop")

# Process Parvej et al. (2025) data (already contains means and location points)
# Rename columns to match model data for binding
exp_data_formatted <- exp_yield_reduction %>%
  rename(Defolication_DVI = Defoliation_DVI, 
         Defolication_PCT = Defoliation_PCT,
         Yield_R = Mean_YR) %>%
  mutate(Defolication_DVI = as.factor(Defolication_DVI))

# Calculate means for the bars (Parvej)
df_bars_parvej <- exp_data_formatted %>%
  group_by(Defolication_DVI, Defolication_PCT, Source) %>%
  summarise(Mean_YR = mean(Yield_R), .groups = "drop")

# Combine all bars and all points
all_bars <- bind_rows(df_bars_models, df_bars_parvej)
all_points <- bind_rows(
  df_long %>% select(Defolication_DVI, Defolication_PCT, Source, Yield_R),
  exp_data_formatted %>% select(Defolication_DVI, Defolication_PCT, Source, Yield_R)
)

# Set Factor levels for consistent coloring/ordering
source_levels <- c("UTR", "Partitioning", "Parvej et al. (2025)")
all_bars$Source <- factor(all_bars$Source, levels = source_levels)
all_points$Source <- factor(all_points$Source, levels = source_levels)

# --- 2. Plotting ---

yield_reduction_fig <- ggplot() +
  geom_bar(data = all_bars, 
           aes(x = as.factor(Defolication_PCT), y = Mean_YR, fill = Source),
           stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(data = all_points, 
             aes(x = as.factor(Defolication_PCT), y = Yield_R, 
                 group = Source, fill = Source),   # <-- fill inside aes()
             position = position_dodge(width = 0.8), 
             shape = 21, size = 1.5, color = "black", alpha = 0.7) +
  facet_wrap(~Defolication_DVI, nrow = 2, 
             labeller = as_labeller(c("1.4" = "R4 (DVI = 1.4)", 
                                      "1.5" = "R5 (DVI = 1.5)"))) +
  scale_fill_manual(values = c("UTR" = "#E65F00", 
                               "Partitioning" = "#0072B2", 
                               "Parvej et al. (2025)" = "#009E73")) +
  labs(x = "Defoliation Percentage (%)",
       y = "Yield Reduction (%)",
       fill = "Source") +          # <-- overrides the legend title
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"))

print(yield_reduction_fig)

# --- 3. Save ---
ggsave('Fig-defoliation_comparison.png', 
       plot = yield_reduction_fig, 
       width = 6, height = 6, dpi = 600)















