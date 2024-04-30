# Clear workspace
rm(list=ls())
# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Load Data
load('organ_biomass_plot_2021_ld11.RData')
load('organ_biomass_plot_2022_ld11.RData')
load('organ_biomass_plot_2022_ld10.RData')
load('organ_biomass_plot_2021_ld11_partitioning.RData')
load('organ_biomass_plot_2022_ld11_partitioning.RData')
load('organ_biomass_plot_2022_ld10_partitioning.RData')

years <- c('2021', '2022')
figs.UTR <- list()
figs.partitioning <- list()
# A function to format decimal places
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

calculate_msqe <- function(year, biocro_organ_biomass_tall, field_organ_biomass_tall, model = "partitioning"){
  sampling.times <- unique(field_organ_biomass_tall$time)
  # delete the first data because the error is 0 for both models
  sampling.times <- sampling.times[-1]
  field_organ_biomass_tall <- field_organ_biomass_tall[-which(field_organ_biomass_tall$time==field_organ_biomass_tall$time[1]),]
  print(sampling.times)
  simulation.results <- biocro_organ_biomass_tall[which(biocro_organ_biomass_tall$time %in% sampling.times), ]
 
  merged.t <- merge(field_organ_biomass_tall, simulation.results, by = c("time", "Organ"), all = T)
  merged.t$diff = merged.t$biomass.x - merged.t$biomass.y
  print(paste0(year,':' ,specify_decimal(sum((merged.t$diff)^2 / length(merged.t$diff)),2)))
  # print(merged.t[ ,c('time','Organ', 'diff')])
  }

calculate_msqe(2021,biocro_organ_biomass_tall_2021_ld11, field_organ_biomass_tall_2021_ld11, 'UTR')
calculate_msqe(2022,biocro_organ_biomass_tall_2022_ld11, field_organ_biomass_tall_2022_ld11, 'UTR')
calculate_msqe(2022,biocro_organ_biomass_tall_2022_ld10, field_organ_biomass_tall_2022_ld10, 'UTR')

calculate_msqe(2021, biocro_organ_biomass_tall_2021_ld11_partitioning, field_organ_biomass_tall_2021_ld11_partitioning, 'partitioning')
calculate_msqe(2022, biocro_organ_biomass_tall_2022_ld11_partitioning, field_organ_biomass_tall_2022_ld11_partitioning, 'partitioning')
calculate_msqe(2022, biocro_organ_biomass_tall_2022_ld10_partitioning, field_organ_biomass_tall_2022_ld10_partitioning, 'partitioning')
