# Clear workspace
rm(list=ls())
# Set working directory to location of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Load Data
load('SoyFACE_results_and_measurements.RData')

# A function to format decimal places
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

# years
years <- c('2002','2004','2005','2006')

# A function to calculate mean squared error
calculate_msqe <- function(year, bio, mea, model = 'Thornley'){ # bio: biocro results, mea: measurement
# year <- '2002'
# bio <- results[[1]]
# mea <- ExpBiomass[[1]]
  # unify the colume names
  mea$DOY <- mea$DOY + 0.5
  names(mea)[names(mea) == "DOY"] <- "time"
  mea <- mea[,c('time', 'Leaf', 'Stem', 'Pod')]
  if (model == 'partitioning'){
    names(bio)[names(bio) == "Grain"] <- "Pod"
  }
  bio <- bio[which(bio$time %in% mea$time), c('time', 'Leaf', 'Stem', 'Pod')]
  diff <- mea[,-1] - bio[,-1]
  rmse <- sum(diff^2) / (dim(diff)[1] * (dim(diff)[2]))
  print(paste0(year ,': ' , specify_decimal(rmse,2)))
}

for (i in 1:length(years)){
  calculate_msqe(years[i], results[[i]], ExpBiomass[[i]])
  calculate_msqe(years[i], results.elevCO2[[i]], ExpBiomass.elevCO2[[i]])
}
