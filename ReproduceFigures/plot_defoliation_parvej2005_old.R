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
defolicate_pcts <- c(0.95)
defolicate_DVIs <- c(1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9)

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
    defoliate_pct <- 0.99
    pod_reduction <- calculate_pod_reduction(defoliate_dvi, yr, defoliate_pct)
    print(pod_reduction)
    # Append each result to the list
    reduction_df[[length(reduction_df) + 1]] <- list(
      year           = pod_reduction$year,
      Defolication_DVI = defoliate_dvi,
      # Defolication   = defoliate_pct,
      UTR            = pod_reduction$utr * 100,
      Partitioning   = pod_reduction$partitioning * 100
    )
  } 
}

# Flatten reduction_df
df <- bind_rows(lapply(reduction_df, function(r) {
  if (any(sapply(r, length) == 0)) return(NULL)
  as.data.frame(r)
}))
df$Defolication_DVI <- as.factor(df$Defolication_DVI)

# Observed data from Fehr, 1977
obs_df <- data.frame(
  Defolication_DVI  = as.factor(c(1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9)),
  Source  = "Fehr (1977)",
  Mean_YR = c(21, 27, 55, 82, 53, 3, 0)
)

# Compute per-year long format for points (Source 1 & 2)
df_points <- df %>%
  tidyr::pivot_longer(cols = c(# "UTR1", "UTR2", 
    "UTR", "Partitioning"),
    names_to = "Source", values_to = "Yield_R")

# Compute mean across years for bars (Source 1 & 2)
df_bars <- df_points %>%
  group_by(Defolication_DVI, Source) %>%
  summarise(Mean_YR = mean(Yield_R), .groups = "drop")

# Combine bars from all 3 sources
all_bars <- bind_rows(df_bars, obs_df)
all_bars$Source <- factor(all_bars$Source, 
                          levels = c(
                            "UTR", "Partitioning", "Fehr (1977)"))

# Plot
df_points$Source <- factor(df_points$Source, 
                           levels = c(
                             "UTR", "Partitioning", "Fehr (1977)"))


# Manually compute x positions to match dodged bars exactly
# With 3 groups and dodge width 0.8, offsets are: -0.267, 0, +0.267
# Convert Defolication to numeric in all_bars too
dodge_offset <- data.frame(
  Source = c("UTR", "Partitioning", "Fehr (1977)"),
  offset = c(-0.267, 
             0, 
             0.267)
)

all_bars_pos <- all_bars %>%
  left_join(dodge_offset, by = "Source") %>%
  mutate(x_pos = as.numeric(as.character(Defolication_DVI)) * 10 + offset)
# Scale Defolication to 1,2,3,4 index instead
Defolication_DVI <- c(1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9)

all_bars_pos <- all_bars %>%
  left_join(dodge_offset, by = "Source") %>%
  mutate(x_idx = match(as.numeric(as.character(Defolication_DVI)), defolicate_DVIs),
         x_pos = x_idx + offset)

df_points_pos <- df_points %>%
  left_join(dodge_offset, by = "Source") %>%
  mutate(x_idx = match(as.numeric(as.character(Defolication_DVI)), defolicate_DVIs),
         x_pos = x_idx + offset)

yield_reduction_fig <- ggplot() +
  geom_bar(data = all_bars_pos,
           aes(x = x_pos, y = Mean_YR, fill = Source),
           stat = "identity", width = 0.25) +
  geom_jitter(data = df_points_pos,
              aes(x = x_pos, y = Yield_R, fill = Source),
              width = 0.04, height = 0,
              shape = 21, size = 2.5, color = "black") +
  scale_fill_manual(name = "Source",
                    values = c(# "UTR1"                    = "#E69F00",
                      # "UTR2"                    = "#E67F00",
                      "UTR"                       = "#E65F00",
                      "Partitioning"              = "#0072B2",
                      "Proulx and Naeve (2009)"   = "#009E73"),
                    breaks = c(# "UTR1", "UTR2", 
                      "UTR", "Partitioning", "Proulx and Naeve (2009)")) +
  scale_x_continuous(breaks = 1:7, labels = as.character(Defolication_DVI)) +
  labs(title = "Yield Reduction by Defoliation Percentage",
       x = "100% Defoliation DVI",
       y = "Yield Reduction (%)") +
  theme_bw() +
  theme(legend.position = "bottom")

print(yield_reduction_fig)

ggsave('Fig-defoliate_yield_reduction.png', 
       plot = yield_reduction_fig, 
       width = 6,
       height = 4,
       units = "in",
       dpi = 600
)


