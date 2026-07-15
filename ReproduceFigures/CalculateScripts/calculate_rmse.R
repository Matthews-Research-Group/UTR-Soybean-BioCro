# A function to format decimal places
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
# A function to calculate mse
calculate_rmse <- function(year, biocro_organ_biomass_tall, field_organ_biomass_tall){
  sampling.times <- unique(field_organ_biomass_tall$time)
  # delete the first data because the error is 0 for both models
  sampling.times <- sampling.times[-1]
  field_organ_biomass_tall <- field_organ_biomass_tall[which(field_organ_biomass_tall$time %in% sampling.times),]
  simulation.results <- biocro_organ_biomass_tall[which(biocro_organ_biomass_tall$time %in% sampling.times), ]
  
  merged.t <- merge(field_organ_biomass_tall, simulation.results, by = c("time", "Organ"), all = T)
  merged.t$diff = merged.t$biomass.x - merged.t$biomass.y
  rmse <- sqrt(mean((merged.t$diff)^2))
  print(paste0(year,' rmse:' , specify_decimal(rmse, 2)))
}

# A function to calculate msqe
calculate_rmse_soyFACE <- function(year, biocro_organ_biomass_tall, field_organ_biomass_tall){
  sampling.times <- unique(field_organ_biomass_tall$time)
  # delete the first data because the error is 0 for both models
  sampling.times <- sampling.times[-1]
  field_organ_biomass_tall <- field_organ_biomass_tall[which(field_organ_biomass_tall$time %in% sampling.times),]
  simulation.results <- biocro_organ_biomass_tall[which(biocro_organ_biomass_tall$time %in% sampling.times & 
                                                          biocro_organ_biomass_tall$Organ %in% c('Leaf', 'Stem', 'Pod')),] # SoyFACE does not have root data
  
  merged.t <- merge(field_organ_biomass_tall, simulation.results, by = c("time", "Organ"), all = T)
  merged.t$diff = merged.t$biomass.x - merged.t$biomass.y
  rmse <- sqrt(mean((merged.t$diff)^2))
  print(paste0(year,' rmse:' , specify_decimal(rmse, 2)))
  # print(merged.t)
  # print(paste0(year,' rmsqe:' , specify_decimal(sqrt(msqe),2)))
  # print(merged.t[ ,c('time','Organ', 'diff')])
}
