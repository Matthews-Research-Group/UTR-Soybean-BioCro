# define a function to update differential values
update_differential_quantities <- function(r, updated_values, UPDATE_PARAMETERS, model, pod_removal_pct){
  ## Define remaining percentages
  leaf_remaining_percentage <- 1.0
  stem_remaining_percentage <- 1.0
  pod_remaining_percentage <- 1-pod_removal_pct
  
  last_row <- nrow(r)
  differential_quantities_just_before_podrm <-
    as.list(r[last_row, updated_values])
  
  differential_quantities_just_after_podrm <-
    differential_quantities_just_before_podrm
  
  if (model == 'utr'){
    differential_quantities_just_after_podrm$Leaf_substrate_carbon <-
      differential_quantities_just_before_podrm$Leaf_substrate_carbon * leaf_remaining_percentage
    differential_quantities_just_after_podrm$Leaf_structural_carbon <-
      differential_quantities_just_before_podrm$Leaf_structural_carbon * leaf_remaining_percentage # Could be changed to a different percentage
    differential_quantities_just_after_podrm$Pod_structural_carbon <-
      differential_quantities_just_before_podrm$Pod_structural_carbon * pod_remaining_percentage # Could be changed to a different percentage
    
    differential_quantities_just_after_podrm$Stem_substrate_carbon <-
      differential_quantities_just_before_podrm$Stem_substrate_carbon * stem_remaining_percentage
    differential_quantities_just_after_podrm$Stem_structural_carbon <-
      differential_quantities_just_before_podrm$Stem_structural_carbon * stem_remaining_percentage
    differential_quantities_just_after_podrm$Pod_structural_carbon <-
      differential_quantities_just_before_podrm$Pod_structural_carbon * pod_remaining_percentage
  } 
  else if (model == 'partitioning'){
    differential_quantities_just_after_podrm$Leaf <-
      differential_quantities_just_before_podrm$Leaf * leaf_remaining_percentage
    
    differential_quantities_just_after_podrm$Stem <-
      differential_quantities_just_before_podrm$Stem * stem_remaining_percentage
    
    differential_quantities_just_after_podrm$Shell <-
      differential_quantities_just_before_podrm$Shell * pod_remaining_percentage
    
    differential_quantities_just_after_podrm$Grain <-
      differential_quantities_just_before_podrm$Grain * pod_remaining_percentage
  }
  
  if(UPDATE_PARAMETERS){
    differential_quantities_just_after_podrm$DVI <-
      differential_quantities_just_before_podrm$DVI - 0.2
  }
  return(differential_quantities_just_after_podrm)
}

# define a function to calculate yield reduction from pod removal
calculate_pod_reduction <- function(podrm_dvi, yr, pod_removal_pct){
  weather.growingseason <- weather[sd.idx: hd.ind,]
  RootVals <- data.frame("DOY"=ExpBiomass$DOY[3], "Root"=0.17*sum(ExpBiomass[5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  soybean$parameters$time_zone_offset <- -6
  ###### No Pod Removal scenario #####
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
  
  # print the DVI or DOY of seed growth
  print(paste0("Pod Growth Start DVI: ", result_partitioning_no_podrm$DVI[which(result_partitioning_no_podrm$Grain>0.001)[1]]))
  print(paste0("Pod Growth Start DOY: ", result_partitioning_no_podrm$doy[which(result_partitioning_no_podrm$Grain>0.001)[1]]))
  
  ###### Pod Removal scenario #####
  pod.rm.ind <- which.min(abs(result_utr_no_podrm$DVI - podrm_dvi)) + sd.idx
  weather.growingseason_1 <- weather[sd.idx:pod.rm.ind,]
  weather.growingseason_2 <- weather[pod.rm.ind:hd.ind,]
  
  #### Before Pod Removal #### 
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
  #### After Pod Removal #### 
  ### Without Parameter Update ####
  ## UTR ##
  UPDATE_PARAMETERS <- FALSE
  differential_quantities_just_after_podrm_utr <- update_differential_quantities(result_utr_podrm1, 
                                                                                 names(initial_values), 
                                                                                 UPDATE_PARAMETERS, 'utr',
                                                                                 pod_removal_pct)
  
  soybean_utr_optsolver_podrm2 <- partial_run_biocro(differential_quantities_just_after_podrm_utr,
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
  parameters_after_podrm <- orignal_utr_params
  
  result_utr_podrm2 <- soybean_utr_optsolver_podrm2(parameters_after_podrm$Value)
  
  ## Partitioning model##
  differential_quantities_just_after_podrm_partitioning <- update_differential_quantities(result_partitioning_podrm1, 
                                                                                          names(soybean$initial_values),
                                                                                          UPDATE_PARAMETERS, 'partitioning',
                                                                                          pod_removal_pct)
  
  result_partitioning_podrm2 <- run_biocro(differential_quantities_just_after_podrm_partitioning,
                                           soybean$parameters,
                                           weather.growingseason_2,
                                           soybean$direct_modules,
                                           soybean$differential_modules,
                                           soybean$ode_solver,
                                           verbose = FALSE)
  
  result_utr_podrm <- rbind(result_utr_podrm1[seq_len(nrow(result_utr_podrm1) - 1), ], result_utr_podrm2)
  result_partitioning_podrm <- rbind(result_partitioning_podrm1[seq_len(nrow(result_partitioning_podrm1) - 1), ], result_partitioning_podrm2)
  
  # Check the DOY of the assumed Pod Removal DVI
  print(paste0('DVI = ', podrm_dvi, " is on DOY = ", 
               result_utr_no_podrm$doy[which.min(abs(result_utr_no_podrm$DVI - podrm_dvi))]))
  
  
  exp_first_seed_idx <- which(ExpBiomass$Seed>0)[1]
  exp_first_seed_doy <- ExpBiomass$DOY[exp_first_seed_idx]
  exp_before_seed_doy <- ExpBiomass$DOY[exp_first_seed_idx-1]
  exp_first_seed_dvi <- result_utr_no_podrm$DVI[which(result_utr_no_podrm$doy==exp_first_seed_doy)[1]]
  exp_before_seed_dvi <- result_utr_no_podrm$DVI[which(result_utr_no_podrm$doy==exp_before_seed_doy)[1]]
  
  # Method 1: Calculate Seed by taking out the final shell
  final_shell <- (tail(ExpBiomass$Pod, 1) - tail(ExpBiomass$Seed, 1))
  no_podrm_seed <- tail(result_utr_no_podrm$Pod, 1) - final_shell
  podrm_seed <- tail(result_utr_podrm$Pod, 1) - final_shell * (1-pod_removal_pct)
  utr_yield_reduction1 <- 1 - podrm_seed/no_podrm_seed

  # Method 2: Calculate Seed by taking out the shell at Pod Removal
  pod_at_removal_time <- tail(result_utr_podrm1$Pod, 1)
  no_podrm_seed <- tail(result_utr_no_podrm$Pod, 1) - pod_at_removal_time
  podrm_seed <- tail(result_utr_podrm$Pod, 1) - pod_at_removal_time * (1-pod_removal_pct)
  utr_yield_reduction2 <- 1 - podrm_seed/no_podrm_seed
  
  # Method 2_2: Calculate Seed by taking out the shell when shading treatment started
  if(podrm_dvi>1.45){
    shell_at_podrm_start_time <- result_utr_podrm1$Pod[which.min(abs(result_utr_podrm1$DVI - 1.45))]
    no_podrm_seed <- tail(result_utr_no_podrm$Pod, 1) - shell_at_podrm_start_time
    podrm_seed <- tail(result_utr_podrm$Pod, 1) - shell_at_podrm_start_time
    utr_yield_reduction2 <- 1 - podrm_seed/no_podrm_seed
  }
  
  # Method 3: Assume seed:pod ratio to be fixed
  utr_yield_reduction3 <- 1 - tail(result_utr_podrm$Pod, 1)/tail(result_utr_no_podrm$Pod, 1)
  
  # Partitioning Model Yield Change
  result_partitioning_podrm$Pod <- result_partitioning_podrm$Grain  + result_partitioning_podrm$Shell
  result_partitioning_no_podrm$Pod <- result_partitioning_no_podrm$Grain  + result_partitioning_no_podrm$Shell
  partitioning_yield_reduction <- 1 - tail(result_partitioning_podrm$Pod, 1) / tail(result_partitioning_no_podrm$Pod, 1)
  
  return(list(year=yr, 
              exp_before_seed_dvi = exp_before_seed_dvi, 
              exp_first_seed_dvi = exp_first_seed_dvi,
              utr=utr_yield_reduction3, 
              partitioning=partitioning_yield_reduction))
}


reduction_df <- list()

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
  
  for(pod_removal_pct in c(0.2, 0.4, 0.6, 0.7)){
    pod_removal_dvi <- 1.6
    pod_reduction <- calculate_pod_reduction(pod_removal_dvi, yr, pod_removal_pct)
    print(pod_reduction)
    # Append each result to the list
    reduction_df[[length(reduction_df) + 1]] <- list(
      year     = pod_reduction$year,
      Pod_RM   = pod_removal_pct,
      DVI_before_seed = pod_reduction$exp_before_seed_dvi,
      DVI_first_seed = pod_reduction$exp_first_seed_dvi,
      # UTR1      = pod_reduction$utr1 * 100,
      # UTR2      = pod_reduction$utr2 * 100,
      UTR      = pod_reduction$utr * 100,
      Partitioning  = pod_reduction$partitioning * 100
    )
  } 
}

# Flatten reduction_df
df <- bind_rows(lapply(reduction_df, function(r) {
  if (any(sapply(r, length) == 0)) return(NULL)
  as.data.frame(r)
}))
df$Pod_RM <- as.factor(df$Pod_RM)

# Observed data from Proulx and Naeve (2009)
obs_df <- data.frame(
  Pod_RM  = as.factor(c(0.2, 0.4, 0.6, 0.7)),
  Source  = "Proulx and Naeve (2009)",
  Mean_YR = c(11.7, 33.9, 50.1, 52.3)
)

# Compute per-year long format for points (Source 1 & 2)
df_points <- df %>%
  tidyr::pivot_longer(cols = c("UTR", "Partitioning"),
                      names_to = "Source", values_to = "Yield_R")

# Compute mean across years for bars (Source 1 & 2)
df_bars <- df_points %>%
  group_by(Pod_RM, Source) %>%
  summarise(Mean_YR = mean(Yield_R), .groups = "drop")

# Combine bars from all 3 sources
all_depod_bars <- bind_rows(df_bars, obs_df)
all_depod_bars$Source <- factor(all_depod_bars$Source, 
                          levels = c("Proulx and Naeve (2009)",
                            "UTR", "Partitioning"))

# Plot
df_points$Source <- factor(df_points$Source, 
                           levels = c("Proulx and Naeve (2009)",
                                      "UTR", "Partitioning"))


# Manually compute x positions to match dodged bars exactly
# With 3 groups and dodge width 0.8, offsets are: -0.267, 0, +0.267
# Convert Pod_RM to numeric in all_depod_bars too
dodge_offset <- data.frame(
  Source = c("Proulx and Naeve (2009)", "UTR", "Partitioning"),
  offset = c(-0.267, 
             #-0.1335, 
             0, 
             #0.1335, 
             0.267)
)

all_depod_bars_pos <- all_depod_bars %>%
  left_join(dodge_offset, by = "Source") %>%
  mutate(x_pos = as.numeric(as.character(Pod_RM)) * 10 + offset)
# Scale Pod_RM to 1,2,3,4 index instead
pod_levels <- c(0.2, 0.4, 0.6, 0.7)

all_depod_bars_pos <- all_depod_bars %>%
  left_join(dodge_offset, by = "Source") %>%
  mutate(x_idx = match(as.numeric(as.character(Pod_RM)), pod_levels),
         x_pos = x_idx + offset)

df_points_pos <- df_points %>%
  left_join(dodge_offset, by = "Source") %>%
  mutate(x_idx = match(as.numeric(as.character(Pod_RM)), pod_levels),
         x_pos = x_idx + offset)

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
                    breaks = c("Proulx and Naeve (2009)", "UTR", "Partitioning")) +
  scale_x_continuous(breaks = 1:4, labels = c("20%", "40%", "60%", "70%")) +
  labs(x = "Pod Removal (%)",
       y = "Yield Reduction (%)") +
  theme_bw() +
  theme(legend.position = "bottom")

print(yield_reduction_fig)

ggsave('Fig-pod_removal_yield_reduction.png', 
       plot = yield_reduction_fig, 
       width = 6,
       height = 4,
       units = "in",
       dpi = 600
)

dvi_ranges <- reduction_df %>%
  bind_rows() %>%
  select(year, DVI_before_seed, DVI_first_seed) %>%
  group_by(year) %>%
  slice(1) %>%
  ungroup()

