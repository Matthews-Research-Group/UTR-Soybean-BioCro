# Clear workspace
rm(list=ls())
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')
utr_params <- data.frame(param_names = arg_names_short,
                         SoyFACE = optim_params_short_SoyFACE, 
                         EnergyFarm = optim_params_short)
utr_params$change <- as.integer(100 * (utr_params$EnergyFarm - utr_params$SoyFACE) / utr_params$SoyFACE)
write.csv(utr_params,"params_comparison.csv", row.names = FALSE)