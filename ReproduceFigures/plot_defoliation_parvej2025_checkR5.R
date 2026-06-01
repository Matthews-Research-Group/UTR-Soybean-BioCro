library(ggplot2)
library(tidyr)
library(dplyr)

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
              utr_nd <- result_utr_no_defoliation,
              partitioning_nd <- result_partitioning_no_defoliation,
              utr=result_utr_defoliation, 
              partitioning=result_partitioning_defoliation))
}


reduction_df <- list()
defoliate_pcts <- c(0, 0.25, 0.5, 0.75, 0.999)

utr_rs <- list()
partitioning_rs <- list()

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
    
  for(i in 1:length(defoliate_pcts)){
      defoliate_dvi <- 1.5
      pod_reduction <- calculate_pod_reduction(defoliate_dvi, yr, defoliate_pcts[i])
      # Append each result to the list
      utr_rs[[i]]             <- pod_reduction$utr 
      partitioning_rs[[i]]     <- pod_reduction$partitioning 
  }
  
  # Helper to reshape one dataframe to long format
  to_long <- function(df, line_group) {
    df %>%
      pivot_longer(cols = c(Leaf, Stem, Root, Pod),
                   names_to = "Organ",
                   values_to = "Biomass") %>%
      mutate(LineGroup = line_group)
  }
  
  # Helper: reshape one df to long, tagging line group, source, and percentage
  to_long <- function(df, line_group, source, pct_label) {
    df %>%
      pivot_longer(cols = c(Leaf, Stem, Root, Pod),
                   names_to = "Organ",
                   values_to = "Biomass") %>%
      mutate(LineGroup = line_group,
             Source = source,
             Pct = pct_label)
  }
  
  # Build one big combined dataframe across all comparisons
  all_data <- bind_rows(lapply(2:length(defoliate_pcts), function(i) {
    pct_label <- paste0(defoliate_pcts[i]*100, "%")
    bind_rows(
      to_long(utr_rs[[1]],          "Baseline",   "UTR",          pct_label),
      to_long(utr_rs[[i]],          "Defoliated", "UTR",          pct_label),
      to_long(partitioning_rs[[1]], "Baseline",   "Partitioning", pct_label),
      to_long(partitioning_rs[[i]], "Defoliated", "Partitioning", pct_label)
    )
  }))
  
  # Keep facet rows in numeric order, columns in the order you want
  all_data <- all_data %>%
    mutate(
      Pct    = factor(Pct, levels = paste0(defoliate_pcts[-1]*100, "%")),
      Source = factor(Source, levels = c("UTR", "Partitioning")),
      Organ  = factor(Organ, levels = c("Leaf", "Stem", "Root", "Pod"))
    )
  
  # Plot
  final_plot <- ggplot(all_data,
                       aes(x = fractional_doy, y = Biomass,
                           color = Organ, linetype = LineGroup)) +
    geom_line(linewidth = 1) +
    facet_grid(rows = vars(Pct), cols = vars(Source)) +   # shared axes by default
    scale_color_manual(values = c("Leaf" = "#117733",
                                  "Stem" = "#999933",
                                  "Root" = "#332288",
                                  "Pod"  = "#882255")) +
    scale_linetype_manual(values = c("Baseline" = "dotted", "Defoliated" = "solid")) +
    labs(x = "Day of Year", y = "Biomass",
         color = NULL, linetype = NULL,
         title = yr) +
    theme_minimal() +
    theme(strip.text.y = element_text(angle = 0))  # horizontal left-side labels
  
  print(final_plot)
  
  # Helper: find the fractional_doy where DVI first reaches a target value
  dvi_to_doy <- function(df, target) {
    idx <- which.min(abs(df$DVI - target))
    df$fractional_doy[idx]
  }
  
  # Build a dataframe of vline positions, one set per facet (Pct x Source)
  vline_df <- bind_rows(lapply(2:length(defoliate_pcts), function(i) {
    pct_label <- paste0(defoliate_pcts[i]*100, "%")
    bind_rows(
      data.frame(Pct = pct_label, Source = "UTR",
                 DVI = c(1.35, 1.5),
                 xint = c(dvi_to_doy(utr_rs[[1]], 1.35),
                          dvi_to_doy(utr_rs[[1]], 1.5))),
      data.frame(Pct = pct_label, Source = "Partitioning",
                 DVI = c(1.35, 1.5),
                 xint = c(dvi_to_doy(partitioning_rs[[1]], 1.35),
                          dvi_to_doy(partitioning_rs[[1]], 1.5)))
    )
  }))
  
  vline_df <- vline_df %>%
    mutate(
      Pct    = factor(Pct, levels = paste0(defoliate_pcts[-1]*100, "%")),
      Source = factor(Source, levels = c("UTR", "Partitioning")),
      DVI    = factor(DVI, levels = c(1.35, 1.5))
    )
  
  # Plot
  final_plot <- ggplot(all_data,
                       aes(x = fractional_doy, y = Biomass,
                           color = Organ, linetype = LineGroup)) +
    geom_vline(data = vline_df,
               aes(xintercept = xint, group = DVI),
               color = ifelse(vline_df$DVI == "1.35", "grey75", "grey45"),
               linewidth = 0.6, inherit.aes = FALSE) +
    geom_line(linewidth = 1) +
    facet_grid(rows = vars(Pct), cols = vars(Source)) +
    scale_color_manual(values = c("Leaf" = "#117733",
                                  "Stem" = "#999933",
                                  "Root" = "#332288",
                                  "Pod"  = "#882255")) +
    scale_linetype_manual(values = c("Baseline" = "dotted", "Defoliated" = "solid")) +
    labs(x = "Day of Year", y = "Biomass",
         color = NULL, linetype = NULL,
         title = yr) +
    theme_minimal() +
    theme(strip.text.y = element_text(angle = 0))
  
  print(final_plot)
}
