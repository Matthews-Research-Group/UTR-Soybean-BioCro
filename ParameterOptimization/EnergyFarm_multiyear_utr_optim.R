## Optimization problem to solve for partitioning parameters for soybean
library(BioCro)
library(lattice)
EF_utr_optim <- function(optim_params_short, 
                                biocro.fun, 
                                ExpData.BM, 
                                ExpData.TNC, 
                                num_rows, 
                                wts){
  # print(optim_params_short)
  cost.avg <- 0
  optim_params <- optim_params_conversion(optim_params_short)
  
  for (i in 1:length(biocro.fun)) {
    gc()
    gro.opt <- match.fun(biocro.fun[[i]])
    result <- biocro.fun[[i]](optim_params)
    if (nrow(result) < num_rows[i]) {
      # print('did not finished simulation.')
      # print(nrow(result))
      # print(num_rows[i])
      cost.avg <- 1e10
      break
    } else{
      # print('finished simulation.')
      TrueValues.BM <- ExpData.BM[[i]] 
      # Predicted values at DOY equal the time in experimental data
      Pred.BM <- data.frame("time" = TrueValues.BM$DOY + 0.5) # Assume the biomass is collected at noon everyday.
      time_inds.BM <- which(result$time %in% Pred.BM$time)
      Pred.BM$Leaf <- result$Leaf[time_inds.BM]
      Pred.BM$Stem <- result$Stem[time_inds.BM]
      Pred.BM$Pod <- result$Pod[time_inds.BM]
      Pred.BM$Root <- result$Root[time_inds.BM]
      cf <- optim_params_short[1]
      Pred.BM$LeafLitter <- cf * result$Leaf_senescence_loss[time_inds.BM]
      Pred.BM$StemLitter <- cf * result$Stem_senescence_loss[time_inds.BM]
      
      # factor to scale experimental and simulated results between 0 and ~1 for all components
      scale.leaf <- max(TrueValues.BM$Leaf)
      scale.stem <- max(TrueValues.BM$Stem)
      scale.pod <- max(TrueValues.BM$Pod)
      scale.Root <- max(TrueValues.BM$Root)
      scale.LeafLitter <- max(TrueValues.BM$LeafLitter)
      scale.StemLitter <- max(TrueValues.BM$StemLitter)
      
      # weighted rmses
      err.leaf <- sum(((Pred.BM$Leaf-TrueValues.BM$Leaf)/scale.leaf)^2)/length(Pred.BM$Leaf)
      err.stem <- sum(((Pred.BM$Stem-TrueValues.BM$Stem)/scale.stem)^2)/length(Pred.BM$Stem)
      err.pod <- sum(((Pred.BM$Pod-TrueValues.BM$Pod)/scale.pod)^2)/length(Pred.BM$Pod)
      err.root <- sum(((Pred.BM$Root-TrueValues.BM$Root)/scale.Root)^2)/length(Pred.BM$Root)
      
      # print(paste0('err.leaf:', err.leaf))
      # print(paste0('Pred.BM$Leaf:', Pred.BM$Leaf))
      # print(paste0('err.stem:', err.stem))
      # print(paste0('err.pod:', err.pod))
      # print(paste0('err.root:', err.root))
      cost <- wts$Leaf*err.leaf + wts$Stem*err.stem + wts$Pod*err.pod + wts$Root*err.root
      
      # add additional weight to the 3 measurements closest to the pod_start_dvi
      pod_start_dvi <- optim_params_short[length(optim_params_short)-1]
      pod_start_time <- result$time[which.min(abs(result$DVI-pod_start_dvi))]
      pod_start_closest_true_value_ind <- which.min(abs(TrueValues.BM$DOY - pod_start_time))
      if(pod_start_closest_true_value_ind>1){
        pod_start_closest_true_value_inds <- c(pod_start_closest_true_value_ind-1,
                                               pod_start_closest_true_value_ind,
                                               pod_start_closest_true_value_ind+1)
      }else{
        pod_start_closest_true_value_inds <- c(pod_start_closest_true_value_ind,
                                               pod_start_closest_true_value_ind+1)
      }
      
      err.around.pod.start <- sum(((Pred.BM$Pod[pod_start_closest_true_value_inds]-
                                      TrueValues.BM$Pod[pod_start_closest_true_value_inds])/scale.pod)^2)/3
      
      cost <- cost + wts$Pod_start * err.around.pod.start
      
      # add litter to the cost function
      err.LeafLitter <- sum(((Pred.BM$LeafLitter-TrueValues.BM$LeafLitter)/scale.LeafLitter)^2)/length(Pred.BM$LeafLitter)
      err.StemLitter <- sum(((Pred.BM$StemLitter-TrueValues.BM$StemLitter)/scale.StemLitter)^2)/length(Pred.BM$StemLitter)
      
      cost <- cost + wts$Litter * (err.LeafLitter + err.StemLitter)
      
      # 2022 TNC (hard coded i to 2 for now)
      # 1 mol C/m^2 * ha/Mg 
      # = (10^9 nmol C) / (m^2) * (10^4 m^2) / (10^9 mg Leaf)
      # = 10^4 nmol/ mg Leaf
      if (i == 2){
        TrueValues.TNC <-ExpData.TNC
        
        # Calculate the time from DOY and hour
        # since "time" in the CSV file don't match with "time" from the simulation result
        TrueValues.TNC$time <- TrueValues.TNC$DOY + TrueValues.TNC$hour/24
        Pred.TNC <- data.frame("time" = TrueValues.TNC$time) # Assume the biomass is collected at noon everyday.
        time_inds.TNC <- which(result$time %in% TrueValues.TNC$time)
        
        Pred.TNC$Leaf <- 10^(4) * # nmol C/ mg Leaf = mol C/ Mg Leaf
          result$Leaf_substrate_carbon[time_inds.TNC] / # mol C/m^2
          result$Leaf[time_inds.TNC] # Mg / ha
        
        Pred.TNC$Stem <- 10^(4) * # nmol C/ mg Stem = mol C/ Mg Stem
          result$Stem_substrate_carbon[time_inds.TNC] / # mol C/m^2
          result$Stem[time_inds.TNC] # Mg / ha
        
        scale.leaf.TNC <- max(TrueValues.TNC$Leaf)
        scale.stem.TNC <- max(TrueValues.TNC$Stem)
        
        # add substrate C data to the cost function
        err.Leaf.TNC <- sum(((Pred.TNC$Leaf-TrueValues.TNC$Leaf)/scale.leaf.TNC)^2)/length(Pred.TNC$Leaf)
        err.Stem.TNC <- sum(((Pred.TNC$Stem-TrueValues.TNC$Stem)/scale.stem.TNC)^2)/length(Pred.TNC$Stem)
        
        # print(paste0('err.Leaf.TNC:', err.Leaf.TNC))
        # print(paste0('err.Stem.TNC:', err.Stem.TNC))
        
        cost <- cost + wts$TNC * (err.Leaf.TNC + err.Stem.TNC)
      }
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

