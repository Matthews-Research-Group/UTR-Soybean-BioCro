## Optimization problem to solve for partitioning parameters for soybean
library(BioCro)
library(lattice)
thornley_2022_optim <- function(optim_params_short, 
                                   biocro.fun, 
                                   ExpData.BM, 
                                   ExpData.TNC, 
                                   num_rows, 
                                   wts2){
  # print(optim_params_short)
  cost.avg <- 0
  optim_params <- optim_params_conversion(optim_params_short)
  
  gc()
  gro.opt <- match.fun(biocro.fun)
  result <- biocro.fun(optim_params)
  
  if (nrow(result) < num_rows) {
    # if simulation does not complete, assign a very high cost and exit
    # print(nrow(result))
    # to test whether the simulation ended at Pod start DVI
    # print(max(result$DVI))
    # print(result$Stem_substrate_carbon[nrow(result)])
    # print(result$Pod_substrate_carbon[nrow(result)])
    # print(result$substrate_transport_Stem_to_Pod[nrow(result)])
    # print(optim_params_short)
    cost <- 1e10
  
  } else{
    TrueValues.BM <- ExpData.BM 
    TrueValues.TNC <-ExpData.TNC

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
    
    # Calculate the time from DOY and hour
    # since "time" in the CSV file don't match with "time" from the simulation result
    TrueValues.TNC$time <- TrueValues.TNC$DOY + TrueValues.TNC$hour/24
    Pred.TNC <- data.frame("time" = TrueValues.TNC$time) # Assume the biomass is collected at noon everyday.
    time_inds.TNC <- which(result$time %in% TrueValues.TNC$time)
    
    # 1 mol C/m^2 * ha/Mg 
    # = (10^9 nmol C) / (m^2) * (10^4 m^2) / (10^9 mg Leaf)
    # = 10^4 nmol/ mg Leaf
    Pred.TNC$Leaf <- 10^(4) * # nmol C/ mg Leaf = mol C/ Mg Leaf
      result$Leaf_substrate_carbon[time_inds.TNC] / # mol C/m^2
      result$Leaf[time_inds.TNC] # Mg / ha
    
    Pred.TNC$Stem <- 10^(4) * # nmol C/ mg Stem = mol C/ Mg Stem
      result$Stem_substrate_carbon[time_inds.TNC] / # mol C/m^2
      result$Stem[time_inds.TNC] # Mg / ha

    # factor to scale experimental and simulated results between 0 and ~1 for all components
    scale.leaf <- max(TrueValues.BM$Leaf)
    scale.stem <- max(TrueValues.BM$Stem)
    scale.pod <- max(TrueValues.BM$Pod)
    scale.Root <- max(TrueValues.BM$Root)
    scale.LeafLitter <- max(TrueValues.BM$LeafLitter)
    scale.StemLitter <- max(TrueValues.BM$StemLitter)
    
    scale.leaf.TNC <- max(TrueValues.TNC$Leaf)
    scale.stem.TNC <- max(TrueValues.TNC$Stem)
    
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
    
    cost <- wts2$Leaf*err.leaf + wts2$Stem*err.stem + wts2$Pod*err.pod + wts2$Root*err.root
    
    # add litter to the cost function
    err.LeafLitter <- sum(((Pred.BM$LeafLitter-TrueValues.BM$LeafLitter)/scale.LeafLitter)^2)/length(Pred.BM$LeafLitter)
    err.StemLitter <- sum(((Pred.BM$StemLitter-TrueValues.BM$StemLitter)/scale.StemLitter)^2)/length(Pred.BM$StemLitter)
    
    cost <- cost + wts2$Litter * (err.LeafLitter + err.StemLitter)

    # add substrate C data to the cost function
    err.Leaf.TNC <- sum(((Pred.TNC$Leaf-TrueValues.TNC$Leaf)/scale.leaf.TNC)^2)/length(Pred.TNC$Leaf)
    err.Stem.TNC <- sum(((Pred.TNC$Stem-TrueValues.TNC$Stem)/scale.stem.TNC)^2)/length(Pred.TNC$Stem)
    
    # print(paste0('err.Leaf.TNC:', err.Leaf.TNC))
    # print(paste0('err.Stem.TNC:', err.Stem.TNC))
    
    cost <- cost + wts2$TNC * (err.Leaf.TNC + err.Stem.TNC)
    
    if(is.nan(cost)){
      cost <- 1e10
    }else{
      cost <- round(100 * cost,2) 
    }
    
    # print(paste0('cost:', cost))
  }
  return(cost)
  
}

