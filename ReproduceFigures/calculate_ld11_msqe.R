# Clear workspace
rm(list=ls())
# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Load Data
years <- c('2021','2022','2023')
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
organ_biomass <- list()
measured_biomass <- list()
figs.UTR <- list()

# A function to format decimal places
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

calculate_msqe <- function(year, biocro_organ_biomass_tall, field_organ_biomass_tall){
  sampling.times <- unique(field_organ_biomass_tall$time)
  # delete the first data because the error is 0 for both models
  sampling.times <- sampling.times[-1]
  field_organ_biomass_tall <- field_organ_biomass_tall[-which(field_organ_biomass_tall$time==field_organ_biomass_tall$time[1]),]
  simulation.results <- biocro_organ_biomass_tall[which(biocro_organ_biomass_tall$time %in% sampling.times), ]
 
  merged.t <- merge(field_organ_biomass_tall, simulation.results, by = c("time", "Organ"), all = T)
  # print(merged.t)
  merged.t$diff = merged.t$biomass.x - merged.t$biomass.y
  msqe <- mean((merged.t$diff)^2)
  print(paste0(year,' msqe:' , specify_decimal(msqe, 2)))
  print(paste0(year,' rmsqe:' , specify_decimal(sqrt(msqe),2)))
  # print(merged.t[ ,c('time','Organ', 'diff')])
  }
for (i in 1:length(years)){
  organ_biomass[[i]] <- loadRData(paste0('organ_biomass_sim_', years[i],'_ld11.RData'))
  measured_biomass[[i]] <- loadRData(paste0('organ_biomass_mea_', years[i],'_ld11.RData'))
  calculate_msqe(years[[i]],organ_biomass[[i]], measured_biomass[[i]])
}


