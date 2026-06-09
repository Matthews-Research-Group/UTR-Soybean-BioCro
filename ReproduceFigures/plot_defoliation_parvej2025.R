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
  
  # if(defoliate_dvi==1.5){
  #   if(defoliate_pct==0.5){utr_yield_reduction3 <- 1- (1-utr_yield_reduction3) * 0.95 }
  #   if(defoliate_pct==0.75){utr_yield_reduction3 <- 1- (1-utr_yield_reduction3) * 0.86}
  #   if(defoliate_pct==0.999){utr_yield_reduction3 <- 1- (1-utr_yield_reduction3) * 0.8 }
  # }


  
  # Method 1: Calculate Seed by taking out the final shell
  final_shell <- (tail(ExpBiomass$Pod, 1) - tail(ExpBiomass$Seed, 1))
  no_defoliation_seed <- tail(result_utr_no_defoliation$Pod, 1) - final_shell
  shade_seed <- tail(result_utr_defoliation$Pod, 1) - final_shell
  utr_yield_reduction1 <- 1 - shade_seed/no_defoliation_seed

  # Method 2: Calculate Seed by taking out the shell when shading treatment started
  estimated_shell_mass <- result_utr_defoliation$Pod[which.min(abs(result_utr_defoliation$DVI - 1.35))]
  no_defoliation_seed <- tail(result_utr_no_defoliation$Pod, 1) - estimated_shell_mass
  defoliation_seed <- tail(result_utr_defoliation$Pod, 1) - estimated_shell_mass
  utr_yield_reduction2 <- 1 - defoliation_seed/no_defoliation_seed
  
  # print('Estimated shell mass:')
  # print(estimated_shell_mass)
  # 
  # full_shell_doy <- result_utr_defoliation$fractional_doy[which.min(abs(result_utr_defoliation1$DVI - 1.45))]     
  # print(paste0('full shell DOY: ', full_shell_doy))
  # 
  # 
  ExpBiomass <- read.csv(file=paste0('../Data/SoyFACE_data/biomasses/',yr, '_ambient_biomass.csv'))
  # closest_exp_idx <- which.min(abs(ExpBiomass$DOY-full_shell_doy))
  # print('Closest experimental Pod mass')
  # print(ExpBiomass$Rep_Mg_per_ha[closest_exp_idx])
  
  
  # Partitioning Model Yield Change
  result_partitioning_defoliation$Pod <- result_partitioning_defoliation$Grain + result_partitioning_defoliation$Shell
  result_partitioning_no_defoliation$Pod <- result_partitioning_no_defoliation$Grain + result_partitioning_no_defoliation$Shell
  partitioning_yield_reduction <- 1 - tail(result_partitioning_defoliation$Pod, 1) / tail(result_partitioning_no_defoliation$Pod, 1)
  
  return(list(year=yr, 
              utr=utr_yield_reduction3, 
              partitioning=partitioning_yield_reduction,
              utr_shell_mass_1 = tail(result_utr_defoliation$Pod, 1) * (1-0.76),
              utr_shell_mass_2 = estimated_shell_mass,
              partitioning_shell = tail(result_partitioning_defoliation$Shell, 1),
              exp_final_shell = tail(ExpBiomass$Rep_Mg_per_ha,1)-tail(ExpBiomass$Seed_Mg_per_ha,1)))
}


reduction_df <- list()
defoliate_pcts <- c(0, 0.25, 0.5, 0.75, 0.999)
defoliate_DVIs <- c(1.35, 1.5)

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
  
  for(defoliate_dvi in defoliate_DVIs){
    for(defoliate_pct in defoliate_pcts){
      pod_reduction <- calculate_pod_reduction(defoliate_dvi, yr, defoliate_pct)
      print(pod_reduction)
      # Append each result to the list
      reduction_df[[length(reduction_df) + 1]] <- list(
        year             = pod_reduction$year,
        Defoliation_DVI  = defoliate_dvi,
        Defoliation_PCT  = defoliate_pct,
        UTR              = pod_reduction$utr * 100,
        Partitioning     = pod_reduction$partitioning * 100,
        utr_shell_mass_1 = pod_reduction$utr_shell_mass_1,
        utr_shell_mass_2 = pod_reduction$utr_shell_mass_2,
        partitioning_shell = pod_reduction$partitioning_shell,
        exp_final_shell  = pod_reduction$exp_final_shell
      )
    }
  } 
}

# Flatten reduction_df
df <- bind_rows(lapply(reduction_df, function(r) {
  if (any(sapply(r, length) == 0)) return(NULL)
  as.data.frame(r)
}))
df$Defoliation_DVI <- as.factor(df$Defoliation_DVI)

# Creating the digitized dataframe from Parvej et al. (2025)
exp_yield_reduction <- data.frame(
  location = c(
    rep("Iowa", 10),    # 5 points for R4, 5 points for R5
    rep("Indiana", 10)  # 5 points for R4, 5 points for R5
  ),
  Defoliation_DVI = c(
    # Iowa
    rep(1.35, 5), rep(1.5, 5),
    # Indiana
    rep(1.35, 5), rep(1.5, 5)
  ),
  Defoliation_PCT = rep(c(0, 25, 50, 75, 99.9), 4),
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

# --- 1. Data Preparation ---

# Process the simulated data frame
# Ensure Defoliation_PCT is on the same scale as exp_yield_reduction (0-100)
df_long <- df %>%
  mutate(Defoliation_PCT = Defoliation_PCT * 100) %>%
  pivot_longer(cols = c("UTR", "Partitioning"), 
               names_to = "Source", 
               values_to = "Yield_R")

# Calculate means for the bars (Models)
df_bars_models <- df_long %>%
  group_by(Defoliation_DVI, Defoliation_PCT, Source) %>%
  summarise(Mean_YR = mean(Yield_R), .groups = "drop")


# Calculate means for the bars (Parvej)
exp_data_formatted <- exp_yield_reduction %>%
  mutate(Defoliation_DVI = as.factor(Defoliation_DVI))

df_bars_parvej <- exp_data_formatted %>%
  group_by(Defoliation_DVI, Defoliation_PCT, Source) %>%
  summarise(Mean_YR = mean(Yield_R), .groups = "drop")

# Combine all bars and all points
all_bars <- bind_rows(df_bars_models, df_bars_parvej)
all_points <- bind_rows(
  df_long %>% select(Defoliation_DVI, Defoliation_PCT, Source, Yield_R),
  exp_data_formatted %>% select(Defoliation_DVI, Defoliation_PCT, Source, Yield_R)
)

# Set Factor levels for consistent coloring/ordering
source_levels <- c("Parvej et al. (2025)", "UTR", "Partitioning")
all_bars$Source <- factor(all_bars$Source, levels = source_levels)
all_points$Source <- factor(all_points$Source, levels = source_levels)

# --- 2. Plotting ---

yield_reduction_fig <- ggplot() +
  geom_bar(data = all_bars, 
           aes(x = as.factor(Defoliation_PCT), y = Mean_YR, fill = Source),
           stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(data = all_points, 
             aes(x = as.factor(Defoliation_PCT), y = Yield_R, 
                 group = Source, fill = Source),   # <-- fill inside aes()
             position = position_dodge(width = 0.8), 
             shape = 21, size = 1.5, color = "black", alpha = 0.7) +
  facet_wrap(~Defoliation_DVI, nrow = 2, 
             labeller = as_labeller(c("1.35" = "R4 (DVI = 1.35)", 
                                      "1.5" = "R5 (DVI = 1.5)"))) +
  scale_fill_manual(values = c("Parvej et al. (2025)" = "#009E73",
                               "UTR" =                  "#E65F00", 
                               "Partitioning" =         "#0072B2")) +
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




shell_table <- do.call(rbind, lapply(reduction_df, function(x) {
  data.frame(
    year                    = x$year,
    Defoliation_DVI         = x$Defoliation_DVI,
    Defoliation_PCT         = x$Defoliation_PCT,
    utr_shell_mass_1        = x$utr_shell_mass_1,
    utr_shell_mass_2        = x$utr_shell_mass_2,
    partitioning_shell      = x$partitioning_shell,
    exp_final_shell         = x$exp_final_shell,
    stringsAsFactors = FALSE
  )
}))

write.csv(shell_table, file = 'shell_table.csv', row.names = FALSE)








