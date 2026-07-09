## Optimization problem to solve for partitioning parameters for soybean
library(BioCro)
library(BioCroWater)
library(UTRSoybeanBML)

multiyear_UTR_optim <- function(
    optim_params_short, biocro.fun, ExpData, num_rows, weights, wts2, RootVals){
  cost.avg <- 0
  optim_params <- optim_params_conversion(optim_params_short)
  
  for (i in 1:length(biocro.fun)) {
    gc()
    gro.opt <- match.fun(biocro.fun[[i]])
    result <- biocro.fun[[i]](optim_params)
    if (nrow(result) < num_rows[i] ||                          # 'simulation did not finish'
        max(result$Leaf_substrate_carbon/result$Leaf) > 1.3) { # TNC content exceed 1.3 * 0.3 = 39%
      cost.avg <- 1e10

      break
      
    } else{ # if simulation finished
      TrueValues <- ExpData[[i]]
      RootValues <- RootVals[[i]]
      
      # Predicted values at DOY equal the DOYs in experimental data
      Pred <- data.frame("DOY"=TrueValues$DOY)
      doy_inds <- which(result$fractional_doy %in% TrueValues$DOY) 
      Pred$Stem <- result$Stem[doy_inds]
      Pred$Leaf <- result$Leaf[doy_inds]
      Pred$Pod <- result$Pod[doy_inds]
      
      Pred.Root <- data.frame("DOY"=RootValues$DOY)
      doy_inds.Root <- which(result$fractional_doy %in% RootValues$DOY) 
      Pred.Root$Vals <- result$Root[doy_inds.Root]
      
      cf <- optim_params_short[1]
      Pred$CummulativeLitter <- cf * (result$Leaf_senescence_loss[doy_inds] + 
                                        result$Stem_senescence_loss[doy_inds])
      
      # weights
      wts <- weights[[i]]
      
      # weighted mses
      err.leaf <- sum(wts$Leaf*((Pred$Leaf-TrueValues$Leaf)^2)/(TrueValues$Leaf+0.1))/length(Pred$Leaf)
      err.stem <- sum(wts$Stem*((Pred$Stem-TrueValues$Stem)^2)/(TrueValues$Stem+0.1))/length(Pred$Stem)
      err.pod <- sum(wts$Pod*((Pred$Pod-TrueValues$Pod)^2)/(TrueValues$Pod+0.1))/length(Pred$Pod)
      err.root <- sum(((Pred.Root$Vals-RootValues$Root)^2)/(RootValues$Root+0.1))/length(Pred.Root$Vals)
     
      cost <- wts2$Leaf*err.leaf + wts2$Stem*err.stem + wts2$Pod*err.pod + wts2$Root*err.root
      
      # add litter to the cost function
      err.litter <- sum(wts$CummulativeLitter*(((Pred$CummulativeLitter-TrueValues$CummulativeLitter)^2))/
                          (TrueValues$CummulativeLitter+0.1))/length(Pred$CummulativeLitter)
      
      cost <- cost + wts2$CummulativeLitter * err.litter
      
      # For testing, change the "parallelType" in 
      # "multiyear_utr_params_optimization.R" to 0,
      # so that the following details can be printed.
      
      # print(paste0('err.leaf: ', err.leaf))
      # print(paste0('err.stem: ', err.stem))
      # print(paste0('err.pod: ', err.pod))
      # print(paste0('err.root: ', err.root))
      # print(paste0('err.litter: ', err.litter))
      # print(paste0('cost: ', cost))
      
      if(is.nan(cost)){
        cost.avg <- 1e10
        break
      }
      
      cost <- round(100 * cost,2) / length(biocro.fun)
      cost.avg <- cost.avg + cost
    }
  }

  return(cost.avg)
  
}