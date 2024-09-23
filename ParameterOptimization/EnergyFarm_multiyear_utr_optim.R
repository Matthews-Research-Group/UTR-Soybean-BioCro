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
      
      # factor to give heavier weights to max values
      weights.leaf <- c(rep(1, length(TrueValues.BM$Leaf)))
      weights.leaf[which.max(TrueValues.BM$Leaf)] <- 2
      weights.stem <- c(rep(1, length(TrueValues.BM$Stem)))
      weights.stem[which.max(TrueValues.BM$Stem)] <- 2
      weights.pod <- c(rep(1, length(TrueValues.BM$Pod)))
      weights.pod[which.max(TrueValues.BM$Pod)] <- 2
      weights.root <- c(rep(1, length(TrueValues.BM$Root)))
      weights.root[which.max(TrueValues.BM$Root)] <- 2

      scale.LeafLitter <- max(TrueValues.BM$LeafLitter)
      scale.StemLitter <- max(TrueValues.BM$StemLitter)
      weights.leaf.litter <- c(rep(1, length(TrueValues.BM$LeafLitter)))
      weights.leaf.litter[which.max(TrueValues.BM$LeafLitter)] <- 2
      weights.stem.litter <- c(rep(1, length(TrueValues.BM$StemLitter)))
      weights.stem.litter[which.max(TrueValues.BM$StemLitter)] <- 2
      
      # weighted rmses
      err.leaf <- sum((weights.leaf*(Pred.BM$Leaf-TrueValues.BM$Leaf)^2)/(TrueValues.BM$Leaf + 0.1))/length(Pred.BM$Leaf)
      err.stem <- sum((weights.stem*(Pred.BM$Stem-TrueValues.BM$Stem)^2)/(TrueValues.BM$Stem + 0.1))/length(Pred.BM$Stem)
      err.pod <- sum((weights.pod*(Pred.BM$Pod-TrueValues.BM$Pod)^2)/(TrueValues.BM$Pod + 0.1))/length(Pred.BM$Pod)
      err.root <- sum((weights.root*(Pred.BM$Root-TrueValues.BM$Root)^2)/(TrueValues.BM$Root+ 0.1))/length(Pred.BM$Root)
      
      # print(paste0('err.leaf:', err.leaf))
      # print(paste0('err.stem:', err.stem))
      # print(paste0('err.pod:', err.pod))
      # print(paste0('err.root:', err.root))
      
      cost <- wts$Leaf*err.leaf + wts$Stem*err.stem + wts$Pod*err.pod + wts$Root*err.root
      
      # add litter to the cost function
      err.LeafLitter <- sum((weights.leaf.litter*(Pred.BM$LeafLitter-TrueValues.BM$LeafLitter)^2)/
                              (TrueValues.BM$LeafLitter + 0.1))/length(Pred.BM$LeafLitter)
      err.StemLitter <- sum((weights.stem.litter*(Pred.BM$StemLitter-TrueValues.BM$StemLitter)^2)/
                              (TrueValues.BM$StemLitter + 0.1))/length(Pred.BM$StemLitter)
      
      # print(paste0('err.LeafLitter:', err.LeafLitter))
      # print(paste0('err.StemLitter:', err.StemLitter))
      
      cost <- cost + wts$Litter * (err.LeafLitter + err.StemLitter)
      
      # 2022 TNC (hard coded i to 2 for now)
      # 1 mol C/m^2 * ha/Mg 
      # = (10^9 nmol C) / (m^2) * (10^4 m^2) / (10^9 mg Leaf)
      # = 10^4 nmol/ mg Leaf
      if (i == 1){
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
        err.Leaf.TNC <- sum(((Pred.TNC$Leaf-TrueValues.TNC$Leaf)^2)/(TrueValues.TNC$Leaf + 10))/length(Pred.TNC$Leaf)
        err.Stem.TNC <- sum(((Pred.TNC$Stem-TrueValues.TNC$Stem)^2)/(TrueValues.TNC$Stem + 10))/length(Pred.TNC$Stem)
        # print(paste0('err.Leaf.TNC: ', err.Leaf.TNC))
        # print(paste0('err.Stem.TNC: ', err.Stem.TNC))
        
        cost <- cost + wts$TNC * (err.Leaf.TNC + err.Stem.TNC)
      }
      if(is.nan(cost)){
        cost.avg <- 1e10
        break
      }
      cost <- round(cost,2) / length(biocro.fun)
      cost.avg <- cost.avg + cost
    }
  }
  # print(cost.avg)
  return(cost.avg)
}

