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

calculate_pod_reduction <- function(podrm_doy, yr, pod_removal_pct){
  pod.rm.ind <- which(weather$doy == podrm_doy)[12] # Specify pod removal time
  weather.growingseason <- weather[sd.idx: hd.ind,]
  weather.growingseason_1 <- weather[sd.idx:pod.rm.ind,]
  weather.growingseason_2 <- weather[pod.rm.ind:hd.ind,]
  # print(weather.growingseason_2$time)
  
  RootVals <- data.frame("DOY"=ExpBiomass$DOY[3], "Root"=0.17*sum(ExpBiomass[5,2:4])) # See Ordonez et al. 2020, https://doi.org/10.1016/j.eja.2020.126130
  
  soybean$parameters$time_zone_offset <- -6
  # ##### No podrm scenario #####
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
  
  
  #### Before podrm #### 
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
  #### After podrm #### 
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
  
  # result_partitioning_podrm$Pod <- result_partitioning_podrm$Grain + result_partitioning_podrm$Shell
  # fig_utr_with_podrm <- plot_SoyFACE_biomass(result_utr_podrm, ExpBiomass, ExpBiomass.std, co2_opt, yr)
  # fig_partitioning_with_podrm <- plot_SoyFACE_biomass(result_partitioning_podrm, ExpBiomass, ExpBiomass.std, co2_opt, yr)
  
  # utr_yield_reduction <- (tail(result_utr_no_podrm$Pod, 1)-tail(result_utr_podrm$Pod, 1))/tail(result_utr_no_podrm$Pod, 1)
  # result_partitioning_no_podrm$Pod <- result_partitioning_no_podrm$Grain + result_partitioning_no_podrm$Shell
  # partitioning_yield_reduction <- (tail(result_partitioning_no_podrm$Pod, 1)-tail(result_partitioning_podrm$Pod, 1))/
  #   tail(result_partitioning_no_podrm$Pod, 1)
  
  no_podrm_seed <- tail(result_utr_no_podrm$Pod, 1)
  podrm_seed <- tail(result_utr_podrm$Pod, 1)
  utr_yield_reduction <- 1 - podrm_seed/no_podrm_seed
  
  partitioning_yield_reduction <- 1 - tail(result_partitioning_podrm$Grain, 1) / tail(result_partitioning_no_podrm$Grain, 1)
  
  podrm.dvi <- result_utr_no_podrm$DVI[which(result_utr_no_podrm$doy==podrm_doy)][1]
  # print(paste0(yr, "-", podrm_doy, "DVI: ", podrm.dvi,"  UTR: ", utr_yield_reduction, "; original: ", partitioning_yield_reduction))
  return(list(year=yr, doy=podrm_doy, dvi=podrm.dvi,  utr=utr_yield_reduction, original=partitioning_yield_reduction))
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
    for(podrm_doy in seq(200, max(ExpBiomass$DOY), by = 10)){
      pod_reduction <- calculate_pod_reduction(podrm_doy, yr, pod_removal_pct)
      # print(pod_reduction)
      # Append each result to the list
      reduction_df[[length(reduction_df) + 1]] <- list(
        year     = pod_reduction$year,
        doy      = pod_reduction$doy,
        DVI      = pod_reduction$dvi,
        Pod_RM   = pod_removal_pct,
        UTR      = pod_reduction$utr * 100,
        Original = pod_reduction$original * 100
      )
    } 
  }
}

# Fix and flatten the list
df <- bind_rows(lapply(reduction_df, function(r) {
  if (any(sapply(r, length) == 0)) return(NULL)
  as.data.frame(r)
}))

df$Pod_RM <- as.factor(df$Pod_RM)
df$year   <- as.factor(df$year)

p1 <- ggplot(df) +
  geom_line(aes(x = doy, y = Original, color = Pod_RM, linetype = "Original"), alpha = 0.5) +
  geom_line(aes(x = doy, y = UTR, color = Pod_RM, linetype = "UTR")) +
  facet_wrap(~ year, ncol = 2) +
  scale_color_brewer(palette = "Set1", name = "Pod Removal %") +
  scale_linetype_manual(name = "Type",
                        values = c("Original" = "dashed", "UTR" = "solid")) +
  guides(linetype = guide_legend(override.aes = list(alpha = c(0.5, 1)))) +
  labs(title = "Pod Reduction by DOY",
       x = "Day of Year (DOY)", y = "Final Pod Reduction") +
  theme_bw() +
  theme(legend.position = "bottom")
print(p1)

# Observed yield reduction data
obs_df <- data.frame(
  DVI     = 1.4,
  Pod_RM  = as.factor(c(0.2, 0.4, 0.6, 0.7)),
  Yield_R = c(11.7, 33.9, 50.1, 52.3)
)

p2 <- ggplot(df) +
  geom_line(aes(x = DVI, y = Original, color = Pod_RM, linetype = "Original"), alpha = 0.5) +
  geom_line(aes(x = DVI, y = UTR, color = Pod_RM, linetype = "UTR")) +
  geom_point(data = obs_df, aes(x = DVI, y = Yield_R, color = Pod_RM), 
             size = 3, shape = 16) +
  facet_wrap(~ year, ncol = 2) +
  scale_color_brewer(palette = "Set1", name = "Pod Removal %") +
  scale_linetype_manual(name = "Type",
                        values = c("Original" = "dashed", "UTR" = "solid")) +
  guides(linetype = guide_legend(override.aes = list(alpha = c(0.5, 1)))) +
  labs(title = "Pod Reduction by DVI",
       x = "Development Index (DVI) of the Pod Removal Day", y = "Yield Reduction (%)") +
  coord_cartesian(xlim = c(1, 2)) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p2)
