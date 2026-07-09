## Optimization problem to solve for partitioning parameters for soybean

multiyear_BioCro_optim_obj <- function(optim_params, biocro.fun, ExpData, num_rows, weights, wts2, RootVals, obs_lai){
  log_file <- "bad_parameters_log.txt"

  tryresult <- tryCatch({
    cost.avg <- calc_cost(optim_params, biocro.fun, ExpData, num_rows, weights, wts2, RootVals,obs_lai) 
    cost.avg # Return the real cost
  }, error = function(e) {
      # Log the bad parameter set if an error occurs
      cat("Error caught with parameters:\n", file = log_file, append = TRUE)
      cat(paste0("Parameters: ", paste(format(optim_params, digits = 16), collapse = ", "), "\n"), file = log_file, append = TRUE)
      cat(paste0("Error message: ", e$message, "\n\n"), file = log_file, append = TRUE)
      return(1e10)  # Penalize this set heavily
  })

  # Check if the result is invalid (NaN, Inf)
  if (is.nan(tryresult) || is.infinite(tryresult)) {
    # Log the bad parameter set
    cat("Invalid result (NaN or Inf) with parameters:\n", file = log_file, append = TRUE)
    cat(paste0("Parameters: ", paste(format(optim_params, digits = 16), collapse = ", "), "\n\n"), file = log_file, append = TRUE)
    
    return(1e10)  # Penalize this set heavily
  }
  return(tryresult)
}

calc_cost<-function(optim_params, biocro.fun, ExpData, num_rows, weights, wts2, RootVals,obs_lai){
  cost.avg <- 0

  for (i in 1:length(biocro.fun)) {

    gro.opt <- match.fun(biocro.fun[[i]])

    result <- gro.opt(optim_params)
    result$time = result$doy + result$hour/24
#add a penalty on kLeaf, kStem, kShell & kGrain
    doy = result$doy
    lai = result$lai
    kLeaf   = result$kLeaf
    kStem   = result$kStem
    kShell  = result$kShell
    kGrain  = result$kGrain
#get the first doys when these terms become non-zero
    doy_leaf  = doy[kLeaf>0.01][1]
    doy_stem  = doy[kStem>0.01][1]
    doy_shell = doy[kShell>0.01][1]
    doy_grain = doy[kGrain>0.01][1]
# Find the last DOY where the obs seed value is 0
# i.e., the DOY before seed is non-zero
    TrueValues <- ExpData[[i]]
    last_zero_doy <- max(TrueValues$DOY[TrueValues$Seed < 0.1])
#add a penalty if 
#1. doy_leaf and doy_stem are separated by more than 5 days
#2. doy_shell or doy_grain occurs before doy_leaf
    if(is.na(doy_leaf)|is.na(doy_stem)|is.na(doy_shell)|is.na(doy_grain)){
       penalty = 9999
    }else if(abs(doy_leaf - doy_stem) > 5){
       penalty = 9999  #now i'm simply not allowing such case 
#    }else if((doy_leaf - doy[1]) > 20 | (doy_leaf - doy[1]) < 10){
#       penalty = 9999  #now i'm simply not allowing such case 
#    }else if(doy_shell < doy_leaf | doy_grain < doy_leaf){
#       penalty = 9999
#    }else if((doy_leaf-doy[1]) < 7 | (doy_stem-doy[1]) < 7){ #make sure leaf does not start very early
#       penalty = 9999
    }else if(doy_grain < last_zero_doy-14){ # doy_grain should not occur before seed is present
       penalty = 9999
    }else{
       penalty = 0
    }

    if (nrow(result) < num_rows[i]) {
      # if simulation does not complete, assign a very high cost and exit
      cost.avg <- 1e10
      break

    } else{

      TrueValues <- ExpData[[i]]

      # Predicted values at DOY equal the DOYs in experimental data
      Pred <- data.frame("DOY"=TrueValues$DOY)
      doy_inds <- which(result$time %in% TrueValues$DOY)
      Pred$Stem <- result$Stem[doy_inds]
      Pred$Leaf <- result$Leaf[doy_inds]
      Pred$Shell  <- result$Shell[doy_inds]
      Pred$Seed <- result$Grain[doy_inds]
      Pred$TotalLitter <- result$LeafLitter[doy_inds] + result$StemLitter[doy_inds]

      doy_inds.Root <- which(result$time %in% RootVals[[i]]$DOY)
      Pred.Root <- result$Root[doy_inds.Root]

      # factor to scale experimental and simulated results between 0 and ~1 for all components
      scale.leaf <- max(TrueValues$Leaf)
      scale.stem <- max(TrueValues$Stem)
      scale.shell  <- max(TrueValues$Shell)
      scale.seed <- max(TrueValues$Seed)
      scale.litter <- max(TrueValues$CumLitter)
      scale.Root <- max(RootVals[[i]]$Root)

      # weights
      wts <- weights[[i]]

      # weighted rmses
      err.stem <- sum(wts$Stem*(((Pred$Stem-TrueValues$Stem)/scale.stem)^2))/length(Pred$Stem)
      err.leaf <- sum(wts$Leaf*(((Pred$Leaf-TrueValues$Leaf)/scale.leaf)^2))/length(Pred$Leaf)
      err.shell <- sum(wts$Shell*(((Pred$Shell-TrueValues$Shell)/scale.shell)^2))/length(Pred$Shell)
      err.seed <- sum(wts$Seed*(((Pred$Seed-TrueValues$Seed)/scale.seed)^2))/length(Pred$Seed)
      err.Root <- sum(((Pred.Root-RootVals[[i]]$Root)/scale.Root)^2)/length(Pred.Root)

      err.litter <- sum(wts$CumLitter*(((Pred$TotalLitter-TrueValues$CumLitter)/scale.litter)^2))/length(Pred$TotalLitter)

      cost <- (wts2$Stem*err.stem + wts2$Shell*err.shell + 
               wts2$Leaf*err.leaf + wts2$Root*err.Root +
               wts2$Seed*err.seed + wts2$TotalLitter* err.litter)

#include errors on LAI
      true_lai  <- obs_lai[[i]]
      doy_inds  <- which(result$time %in% true_lai$Group.1)
      Pred_lai  <- result$lai[doy_inds]
      scale.lai <- max(true_lai$LAI_mean)
      err.lai   <- sum(((Pred_lai - true_lai$LAI_mean)/scale.lai)^2)/length(Pred_lai)
      cost      <- cost + wts2$LAI * err.lai
#

      cost <- round(100 * cost,2) / length(biocro.fun)
#double weight on the second calibration year (2005)
#It seems that Soybean-BioCro's orignal fitting is better in 2005, 
#so I put more weight on this year as well
#      cost <- cost * i
#
      cost.avg <- cost.avg + cost + penalty
      
     #print(c(i,err.stem,err.shell,err.leaf,err.Root,err.seed,err.litter))

    }

  }#end for (i in 1:length(biocro.fun)) 
  return(cost.avg)
}
