library(BioCro)
library(DEoptim)
rm(list=ls())
source('functions_BioCro_run.R') #for my_biocro_run_develop
# Define the RMSE function  
rmse <- function(observed, predicted) {
  # Check if the lengths of the observed and predicted vectors are the same
  if (length(observed) != length(predicted)) {
    stop("The lengths of observed and predicted values must be the same.")
  }      
  
  # Calculate the squared differences
  squared_diff <- (observed - predicted)^2
  
  # Calculate the mean of the squared differences
  mean_squared_diff <- mean(squared_diff)
  
  # Calculate the square root of the mean squared differences (RMSE)
  rmse_value <- sqrt(mean_squared_diff)
  
  return(rmse_value)
} 
#optimize the resistance values
#This function is for the local optimizer
obj<-local({
  iter <- 0
  function(XX,data){
    iter <<- iter + 1
    obs = data[[1]]
    obs_slope = data[[2]]
    years = unique(obs$year)
    output = c()
    for (year in years){
      obs_year = obs[obs$year==year & obs$CO2=="AC",]
      vcmax = 110
      jmax  = 195
      resistance_base    = XX[1]
      para = c(soybean$parameters$b0,soybean$parameters$b1,vcmax,jmax,resistance_base)
      use_varying_sla  = FALSE
      use_pre_opt_vars = TRUE
      results_ctl  = my_biocro_run_develop(para,use_varying_sla,use_pre_opt_vars,year)
      list1_aCO2 = results_ctl[[1]][[1]]

      match_ind = which((list1_aCO2$doy %in% obs_year$DOY) &
                          list1_aCO2$hour==12)

      model_swc = list1_aCO2$soil_water_content[match_ind]
      model_Ci  = list1_aCO2$sunlit_Ci_layer_0[match_ind]
      obs_year = cbind(obs_year,model_swc,model_Ci)
      output = rbind(output,obs_year)

      obs_year = obs[obs$year==year & obs$CO2=="EC",]
      resistance_amplifier = XX[2]
      #set of parameters for eCO2
      vcmax = 110*(1-0.15)  #acclimation
      jmax  = 195*(1-0.05)  #acclimation
      RL25  = soybean$parameters$RL_at_25 * (1+0.37*0.42)
      EleCO2 = 570
      use_varying_sla = TRUE #reduced SLA for eCO2
      #
      para = c(soybean$parameters$b0,soybean$parameters$b1,vcmax,jmax,resistance_base,
               RL25,EleCO2,resistance_amplifier)
      results_ctl = my_biocro_run_develop(para,use_varying_sla,use_pre_opt_vars,year)
      list1_eCO2 = results_ctl[[1]][[1]]
      match_ind = which((list1_eCO2$doy %in% obs_year$DOY) &
                          list1_eCO2$hour==12)

      model_swc = list1_eCO2$soil_water_content[match_ind]
      model_Ci  = list1_eCO2$sunlit_Ci_layer_0[match_ind]
      obs_year = cbind(obs_year,model_swc,model_Ci)
      output = rbind(output,obs_year)
    }
    coefs <- output %>%
      group_by(CO2) %>%
      do({
        fit <- lm(model_Ci ~ model_swc, data = .)
        data.frame(
          slope = coef(fit)[2],
          intercept = coef(fit)[1]
        )
      })
    val = rmse(coefs$slope,obs_slope)
    cat(sprintf(
      "iter %4d | m = %.3f | g0 = %.4f | loss = %.3e\n",
      iter, XX[1], XX[2], val
    ))
    val
  }
})

#this function is for DEoptim
myfun<-function(XX,data){
  obs = data[[1]]
  obs_slope = data[[2]]
  years = unique(obs$year)
  output = c()
  for (year in years){
    obs_year = obs[obs$year==year & obs$CO2=="AC",]
    vcmax = 110
    jmax  = 195
    resistance_base    = XX[1]
    para = c(soybean$parameters$b0,soybean$parameters$b1,vcmax,jmax,resistance_base)
    use_varying_sla  = FALSE
    use_pre_opt_vars = TRUE
    results_ctl  = my_biocro_run_develop(para,use_varying_sla,use_pre_opt_vars,year)
    list1_aCO2 = results_ctl[[1]][[1]]

    match_ind = which((list1_aCO2$doy %in% obs_year$DOY) &
                        list1_aCO2$hour==12)

    model_swc = list1_aCO2$soil_water_content[match_ind]
    model_Ci  = list1_aCO2$sunlit_Ci_layer_0[match_ind]
    obs_year = cbind(obs_year,model_swc,model_Ci)
    output = rbind(output,obs_year)

    obs_year = obs[obs$year==year & obs$CO2=="EC",]
    resistance_amplifier = XX[2]
    vcmax = 110*(1-0.15)  #acclimation
    jmax  = 195*(1-0.05)  #acclimation
    RL25  = soybean$parameters$RL_at_25 * (1+0.37*0.42)
    EleCO2 = 570
    use_varying_sla = TRUE
    para = c(soybean$parameters$b0,soybean$parameters$b1,vcmax,jmax,resistance_base,
             RL25,EleCO2,resistance_amplifier)
    results_ctl = my_biocro_run_develop(para,use_varying_sla,use_pre_opt_vars,year)
    list1_eCO2 = results_ctl[[1]][[1]]
    match_ind = which((list1_eCO2$doy %in% obs_year$DOY) &
                        list1_eCO2$hour==12)

    model_swc = list1_eCO2$soil_water_content[match_ind]
    model_Ci  = list1_eCO2$sunlit_Ci_layer_0[match_ind]
    obs_year = cbind(obs_year,model_swc,model_Ci)
    output = rbind(output,obs_year)
  }
  coefs <- output %>%
    group_by(CO2) %>%
    do({
      fit <- lm(model_Ci ~ model_swc, data = .)
      data.frame(
        slope = coef(fit)[2],
        intercept = coef(fit)[1]
      )
    })
  rmse(coefs$slope,obs_slope)
}

obs_data =
  read.csv("data/doi_10_5061_dryad_g0v62__v20170815/Final_Data_Deposit/Fig_5_ExtendedDataFig8/ABA_soilVWC_middaygasex09_11.csv")

obs_data$xylem_ABA = as.numeric(obs_data$xylem_ABA)
obs_data$Ci  = as.numeric(obs_data$Ci)
obs_data$VWCthru75 = as.numeric(obs_data$VWCthru75)

obs_avg <- obs_data %>%
  group_by(ring,CO2,DOY,year,h2o)%>%
  summarize(
    VWCthru75 = mean(VWCthru75, na.rm = TRUE),
    Ci        = mean(Ci, na.rm = TRUE),
    ABA       = mean(xylem_ABA, na.rm = TRUE),
    .groups = "drop"
  )

obs_avg_CP = obs_avg[obs_avg$h2o=="con",]
obs_avg_CP$VWCthru75 = obs_avg_CP$VWCthru75/100
#for obs's slopes, we still include the replicates of rings
#I think this is more reliable than using just the averages
coefs <- obs_avg_CP %>%
  group_by(CO2) %>%
  do({
    fit <- lm(Ci ~ VWCthru75, data = .)
    data.frame(
      slope = coef(fit)[2],
      intercept = coef(fit)[1]
    )
  })

#average out the rings
#this is for matching DOY only.
obs_avg_CP <- obs_avg_CP %>%
  group_by(CO2,DOY,year) %>%
  summarize(
    VWCthru75 = mean(VWCthru75, na.rm = TRUE),
    Ci        = mean(Ci, na.rm = TRUE),
    ABA       = mean(ABA, na.rm = TRUE),
    .groups = "drop"
  )
obs2use = obs_avg_CP

df = list(obs2use,coefs$slope)
par0  <- c(soybean2$parameters$resistance_base, soybean2$parameters$resistance_amplifier)  
# starting values
lower <- c(0.002, 0.5)
upper <- c(0.2, 10)

## THIS IS FOR DEOPTIM
# maximum number of iterations
max.iter = 1000
stopping_error = 0.5
# Call DEoptim function to run optimization
parVars = c('myfun','df','my_biocro_run_develop','get_CO2','rmse')
cl = makeCluster(8)
# Make sure workers see the same .libPaths() as master (common fix on HPC/mac)
clusterCall(cl, .libPaths, .libPaths())

# Load BioCro + any other packages on every worker
clusterEvalQ(cl, {
  library(BioCro)
  library(dplyr)
  library(tidyr)
  NULL
})

clusterExport(cl, parVars,envir=environment())

de=DEoptim(fn = myfun,
           lower =lower,
           upper = upper,
           data  = df,
           control=list(
                       VTR = stopping_error,itermax=max.iter,
                       parallelType=1,packages=c('BioCro'),parVar=parVars,cl=cl, trace = TRUE)
          )
stopCluster(cl)
saveRDS(de,"optimization_results_Ci_vs_SWC/fit_DEoptim_SLA_reduce.rds")
##
stop()
##LOCAL optimizer
fit <- optim(
  par = par0,
  fn = obj,
  data = df,
  method = "L-BFGS-B",
  lower = lower,
  upper = upper
)
fit_lbfgsb_naive = fit
saveRDS(fit_lbfgsb_naive,"optimization_results_Ci_vs_SWC/fit_lbfgsb_naive_SLA_reduce.rds")
##
# fit$par      # best parameters
# fit$value    # best loss
# fit$convergence

# ggplot(obs_avg_CP, aes(x = VWCthru75, y = Ci, color = CO2)) +
#   geom_point(size = 3) +
#   geom_smooth(method = "lm", se = FALSE) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dotted")+
#   # geom_abline(
#   #   data = coefs,
#   #   aes(slope = slope, intercept = intercept, color = CO2),
#   #   linewidth = 1
#   # ) +
#   theme_minimal(base_size = 16) +
#   # ylim(0,600)+
#   labs(x = "obs",
#        y = "model",
#        color = NULL) +
#   scale_color_manual(
#     values = c("AC" = "steelblue", "EC" = "tomato")
#   )
