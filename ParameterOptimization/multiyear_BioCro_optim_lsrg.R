## Optimization problem to solve for partitioning parameters for soybean
library(BioCro)
# library(profvis)
multiyear_BioCro_optim <- function(optim_params_short, biocro.fun, ExpData, num_rows, weights, wts2, RootVals){
  # print(optim_params_short)
  cost.avg <- 0
  optim_params <- optim_params_conversion(optim_params_short)
  
  for (i in 1:length(biocro.fun)) {
    gc()
    gro.opt <- match.fun(biocro.fun[[i]])
    result <- biocro.fun[[i]](optim_params)
    # print(paste0('num_rows_result:', nrow(result), ' num_rows_weather', num_rows[i]))
    if (nrow(result) < num_rows[i]) {
      # if simulation does not complete, assign a very high cost and exit
      # print(paste0("stop at DVI: ", max(result$DVI)))
      # l <- which.max(result$time)# last index
      # print(paste0("Leaf_substrate_carbon: ", result$Leaf_substrate_carbon[l]))
      # print(paste0("Stem_substrate_carbon: ", result$Stem_substrate_carbon[l]))
      # print(paste0("Root_substrate_carbon: ", result$Root_substrate_carbon[l]))
      # print(paste0("Pod_substrate_carbon: ", result$Pod_substrate_carbon[l]))
      # print(paste0("Pod mass: ", result$Pod[l]))
      # print(paste0("Pod C concentration: ", result$Pod_substrate_carbon[l]/result$Pod[l]))
      # print(paste0("Stem C concentration: ", result$Stem_substrate_carbon[l]/result$Stem[l]))
      # print(paste0("Previous time step Pod_substrate_carbon: ", result$Pod_substrate_carbon[l-1]))
      # print(paste0("Previous time step Pod_utilization: ", result$Pod_utilization_rate[l-1]))
      # print(paste0("Previous time step Stem_to_Pod_transport: ", result$Stem_to_Pod_transport_rate[l-1]))
      # print(paste0("Current time step Pod_utilization: ", result$Pod_utilization_rate[l]))
      # print(paste0("Current time step Stem_to_Pod_transport: ", result$Stem_to_Pod_transport_rate[l]))
      # print('unfinished simulation.')
      cost.avg <- 1e10
      break
      
    } else{
      # print('finished simulation.')
      TrueValues <- ExpData[[i]]
      RootValues <- RootVals[[i]]
      # Predicted values at DOY equal the DOYs in experimental data
      Pred <- data.frame("DOY"=TrueValues$DOY)
      doy_inds <- which(result$time %in% TrueValues$DOY)
      Pred$Stem <- result$Stem[doy_inds]
      Pred$Leaf <- result$Leaf[doy_inds]
      Pred$Pod <- result$Pod[doy_inds]
      
      Pred.Root <- data.frame("DOY"=RootValues$DOY)
      doy_inds.Root <- which(result$time %in% RootValues$DOY)
      Pred.Root$Vals <- result$Root[doy_inds.Root]
      
      cf <- optim_params_short[1]
      Pred$CumLitter <- cf * (result$Leaf_senescence_loss[doy_inds] + result$Stem_senescence_loss[doy_inds])
      
      # factor to scale experimental and simulated results between 0 and ~1 for all components
      scale.leaf <- max(TrueValues$Leaf)
      scale.stem <- max(TrueValues$Stem)
      scale.pod <- max(TrueValues$Pod)
      scale.root <- max(RootValues$Root)
      scale.CumLitter <- max(TrueValues$CumLitter)
      
      # weights
      wts <- weights[[i]]
      
      # weighted rmses
      err.leaf <- sum(wts$Leaf*(((Pred$Leaf-TrueValues$Leaf)/scale.leaf)^2))/length(Pred$Leaf)
      err.stem <- sum(wts$Stem*(((Pred$Stem-TrueValues$Stem)/scale.stem)^2))/length(Pred$Stem)
      err.pod <- sum(wts$Pod*(((Pred$Pod-TrueValues$Pod)/scale.pod)^2))/length(Pred$Pod)
      err.root <- sum(((Pred.Root$Vals-RootValues$Root)/scale.root)^2)/length(Pred.Root$Vals)
      
      cost <- wts2$Leaf*err.leaf + wts2$Stem*err.stem + wts2$Pod*err.pod + wts2$Root*err.root
      
      # add litter to the cost function
      err.litter <- sum(wts$CumLitter*(((Pred$CumLitter-TrueValues$CumLitter)/scale.CumLitter)^2))/length(Pred$CumLitter)
      
      cost <- cost + wts2$CumLitter * err.litter
      
      if(is.nan(cost)){
        cost.avg <- 1e10
        break
      }
      
      cost <- round(100 * cost,2) / length(biocro.fun)
      cost.avg <- cost.avg + cost
    }
  }
  # print(cost.avg)
  return(cost.avg)
  
}