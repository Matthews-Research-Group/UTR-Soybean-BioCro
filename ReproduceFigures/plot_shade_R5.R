# define a function to calculate yield reduction from Shade
calculate_pod_reduction <- function(shade_dvi, yr, shade_pct){
  weather.growingseason <- weather[sd.idx: hd.ind,]
  RootVals <- data.frame("DOY"=ExpBiomass$DOY[3], "Root"=0.17*sum(ExpBiomass[5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  soybean$parameters$time_zone_offset <- -6
  ###### No Shade scenario #####
  soybean_utr_optsolver_no_shade <- partial_run_biocro(initial_values,
                                                       parameters,
                                                       weather.growingseason,
                                                       direct_modules,
                                                       differential_modules,
                                                       solver,
                                                       arg_names,
                                                       verbose = FALSE)
  
  result_utr_no_shade <- soybean_utr_optsolver_no_shade(optim_params_conversion(optim_params_short_SoyFACE))
  result_partitioning_no_shade <- with(soybean, run_biocro(initial_values,
                                                           parameters,
                                                           weather.growingseason,
                                                           direct_modules,
                                                           differential_modules,
                                                           solver,
                                                           verbose = FALSE))
  
  ###### Shade scenario #####
  shade.ind <- which.min(abs(result_utr_no_shade$DVI - shade_dvi)) + sd.idx
  weather.growingseason_1 <- weather[sd.idx:shade.ind,]
  weather.growingseason_2 <- weather[shade.ind:hd.ind,]
  weather.growingseason_2$solar <- weather.growingseason_2$solar * (1-shade_pct)
  #### Before Shade #### 
  soybean_utr_optsolver_shade1 <- partial_run_biocro(initial_values,
                                                     parameters,
                                                     weather.growingseason_1,
                                                     direct_modules,
                                                     differential_modules,
                                                     solver,
                                                     arg_names,
                                                     verbose = FALSE)
  result_utr_shade1 <- soybean_utr_optsolver_shade1(optim_params_conversion(optim_params_short_SoyFACE))
  result_partitioning_shade1 <- with(soybean, run_biocro(initial_values,
                                                         parameters,
                                                         weather.growingseason_1,
                                                         direct_modules,
                                                         differential_modules,
                                                         solver,
                                                         verbose = FALSE))
  #### After Shade#### 
  ## UTR ##
  last_row <- nrow(result_partitioning_shade1)
  differential_quantities_just_before_shade <- as.list(result_utr_shade1
                                                       [last_row, names(initial_values)])
  
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
  
  ## Partitioning model##
  differential_quantities_just_before_shade <- as.list(result_partitioning_shade1
                                                       [last_row, names(soybean$initial_values)])
  result_partitioning_shade2 <- run_biocro(differential_quantities_just_before_shade,
                                           soybean$parameters,
                                           weather.growingseason_2,
                                           soybean$direct_modules,
                                           soybean$differential_modules,
                                           soybean$ode_solver,
                                           verbose = FALSE)
  
  result_utr_shade <- rbind(result_utr_shade1[seq_len(nrow(result_utr_shade1) - 1), ], result_utr_shade2)
  result_partitioning_shade <- rbind(result_partitioning_shade1[seq_len(nrow(result_partitioning_shade1) - 1), ], result_partitioning_shade2)
  # fig_utr_with_shade <- plot_SoyFACE_biomass(result_utr_shade, ExpBiomass, ExpBiomass.std, co2_opt, yr)
  # fig_partitioning_with_shade <- plot_SoyFACE_biomass(result_partitioning_shade, ExpBiomass, ExpBiomass.std, co2_opt, yr)
  
  # # Method 1: Calculate Seed by taking out the final shell
  # final_shell <- (tail(ExpBiomass$Pod, 1) - tail(ExpBiomass$Seed, 1))
  # no_shade_seed <- tail(result_utr_no_shade$Pod, 1) - final_shell
  # shade_seed <- tail(result_utr_shade$Pod, 1) - final_shell
  # utr_yield_reduction1 <- 1 - shade_seed/no_shade_seed
  # 
  # # Method 2: Calculate Seed by taking out the shell when shading treatment started
  # pod_at_shade_start_time <- tail(result_utr_shade1$Pod, 1)
  # no_shade_seed <- tail(result_utr_no_shade$Pod, 1) - pod_at_shade_start_time
  # shade_seed <- tail(result_utr_shade$Pod, 1) - pod_at_shade_start_time
  # utr_yield_reduction2 <- 1 - shade_seed/no_shade_seed
  
  # Method 3: Assume seed:pod ratio to be fixed
  utr_yield_reduction3 <- 1 - tail(result_utr_shade$Pod, 1)/tail(result_utr_no_shade$Pod, 1)
  
  # Partitioning Model Yield Change
  result_partitioning_shade$Pod <- result_partitioning_shade$Grain + result_partitioning_shade$Shell
  result_partitioning_no_shade$Pod <- result_partitioning_no_shade$Grain + result_partitioning_no_shade$Shell
  partitioning_yield_reduction <- 1 - tail(result_partitioning_shade$Pod, 1) / tail(result_partitioning_no_shade$Pod, 1)
  
  return(list(year=yr, 
              #utr1=utr_yield_reduction1,
              #utr2=utr_yield_reduction2,
              utr=utr_yield_reduction3, 
              partitioning=partitioning_yield_reduction))
}


reduction_df <- list()
ratio_df <- list()
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
  
  for(i in (min(which(ExpBiomass$Seed > 0)) - 1):max(which(ExpBiomass$Seed > 0))){
    ratio_df[[length(ratio_df) + 1]] <- list(
      year  = yr,
      DOY   = ExpBiomass[i, 'DOY'],
      Shell = ExpBiomass[i, 'Pod'] - ExpBiomass[i, 'Seed'], 
      Ratio = ExpBiomass[i, 'Seed'] / ExpBiomass[i, 'Pod']
    )
  }
  
  for(shade_pct in c(0.5, 0.6, 0.7, 0.8)){
    shade_dvi <- 1.5
    pod_reduction <- calculate_pod_reduction(shade_dvi, yr, shade_pct)
    print(pod_reduction)
    # Append each result to the list
    reduction_df[[length(reduction_df) + 1]] <- list(
      year     = pod_reduction$year,
      Shade_PCT   = shade_pct,
      # UTR1      = pod_reduction$utr1 * 100,
      # UTR2      = pod_reduction$utr2 * 100,
      UTR       = pod_reduction$utr * 100,
      Partitioning  = pod_reduction$partitioning * 100
    )
  }
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

for(yr in 2021:2024){
  ExpBiomass <- loadRData(paste0('../../energy-farm-biocro/soybean_ld11_biomass_', yr,'/soybean_ld11_biomass_', yr, '.RData'))
  for(j in (min(which(ExpBiomass$seed > 0)) - 1):max(which(ExpBiomass$seed > 0))){
    ratio_df[[length(ratio_df) + 1]] <- list(
      year  = as.character(yr),
      DOY   = ExpBiomass[j, 'doy'],
      Shell = ExpBiomass[j, 'pod'] - ExpBiomass[j, 'seed'], 
      Ratio = ExpBiomass[j, 'seed'] / ExpBiomass[j, 'pod']
    )
  }
}

# Flatten and plot
ratio_df <- bind_rows(lapply(ratio_df, as.data.frame))
ratio_df$year <- as.factor(ratio_df$year)
year_colors <- c("2002" = "#E69F00",
                 "2004" = "#0072B2",
                 "2005" = "#009E73",
                 "2006" = "#D55E00",
                 "2021" = "#CC79A7",
                 "2022" = "#F0E442",
                 "2023" = "#56B4E9",
                 "2024" = "#999999")

ggplot(ratio_df, aes(x = DOY, y = Ratio, color = year)) +
  geom_line() +
  geom_point(size = 2.5) +
  scale_color_manual(name = "Year",
                     values = year_colors) +
  labs(title = "Seed : Pod Ratio Over Growing Season",
       x = "Day of Year (DOY)",
       y = "Seed : Pod Ratio") +
  theme_bw() +
  theme(legend.position = "bottom")

ggplot(ratio_df, aes(x = DOY, y = Shell, color = year)) +
  geom_line() +
  geom_point(size = 2.5) +
  scale_color_manual(name = "Year",
                     values = year_colors) +
  labs(title = "Shell Dry Mass",
       x = "Day of Year (DOY)",
       y = "Shell Dry Mass (Mg / ha)") +
  theme_bw() +
  theme(legend.position = "bottom")

# Observed shading yield reduction data
shade_exp <- data.frame(
  Shade_PCT = as.factor(c(0.5, 0.6, 0.7, 0.8)),
  Source  = "Proulx and Naeve (2009)",
  Mean_YR   = c(18.1, 28.8, 32.3, 37.9)
)


# Flatten reduction_df
reduction_df_flat <- bind_rows(lapply(reduction_df, function(r) {
  r <- r[sapply(r, length) > 0]  # remove NULL elements, keep the row
  as.data.frame(r)
}))
reduction_df_flat$Shade_PCT <- as.factor(reduction_df_flat$Shade_PCT)

# Compute per-year long format for points 
df_points <- reduction_df_flat %>%
  tidyr::pivot_longer(cols = c(# "UTR1", "UTR2", 
                               "UTR", "Partitioning"),
                      names_to = "Source", values_to = "Yield_R")


# Compute mean across years for bars
df_bars <- df_points %>%
  group_by(Shade_PCT, Source) %>%
  summarise(Mean_YR = mean(Yield_R), .groups = "drop")

# Combine bars from all 3 sources
all_bars <- bind_rows(df_bars, shade_exp)
all_bars$Source <- factor(all_bars$Source, 
                          levels = c(#"UTR1", "UTR2", 
                            "UTR", 
                            "Partitioning", 
                            "Proulx and Naeve (2009)"))

# Plot
df_points$Source <- factor(df_points$Source, 
                           levels = c(# "UTR1", "UTR2", 
                                      "UTR", "Partitioning", "Proulx and Naeve (2009)"))

# Manually compute x positions to match dodged bars exactly
# With 3 groups and dodge width 0.8, offsets are: -0.267, 0, +0.267
# Convert Shade_PCT to numeric in all_bars too
dodge_offset <- data.frame(
  Source = c(# "UTR1", "UTR2", 
    "UTR", "Partitioning", 
    "Proulx and Naeve (2009)"),
  offset = c(#-0.267, -0.1335, 0, 0.1335, 0.267
    -0.267,  0, 0.267)
)

all_bars_pos <- all_bars %>%
  left_join(dodge_offset, by = "Source") %>%
  mutate(x_pos = as.numeric(as.character(Shade_PCT)) * 10 + offset)
# Scale Shade_PCT to 1,2,3,4 index instead
shade_levels <- c(0.5, 0.6, 0.7, 0.8)

all_bars_pos <- all_bars %>%
  left_join(dodge_offset, by = "Source") %>%
  mutate(x_idx = match(as.numeric(as.character(Shade_PCT)), shade_levels),
         x_pos = x_idx + offset)

df_points_pos <- df_points %>%
  left_join(dodge_offset, by = "Source") %>%
  mutate(x_idx = match(as.numeric(as.character(Shade_PCT)), shade_levels),
         x_pos = x_idx + offset)

shade_reduction_fig <- ggplot() +
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
                                "UTR"                     = "#E65F00",
                                "Partitioning"                = "#0072B2",
                                "Proulx and Naeve (2009)" = "#009E73"),
                     breaks = c(# "UTR1", "UTR2", 
                                "UTR",
                                "Partitioning", 
                                "Proulx and Naeve (2009)")) +
  scale_x_continuous(breaks = 1:4, labels = c("50%", "60%", "70%", "80%")) +
  labs(x = "Shade (%)",
       y = "Yield Reduction (%)") +
  theme_bw() +
  theme(legend.position = "bottom")

print(shade_reduction_fig)

ggsave('Fig-shade_yield_reduction.png', 
       plot = shade_reduction_fig, 
       width = 6,
       height = 4,
       units = "in",
       dpi = 600
)

ratio_df %>%
  filter(!is.na(Ratio)) %>%
  group_by(year) %>%
  slice_max(DOY) %>%
  ungroup() %>%
  summarise(mean_ratio = mean(Ratio))

ratio_df %>%
  filter(!is.na(Ratio)) %>%
  group_by(year) %>%
  slice_max(DOY) %>%
  ungroup() %>%
  summarise(sd_ratio = sd(Ratio))

