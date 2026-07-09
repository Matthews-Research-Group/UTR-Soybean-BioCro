calculate_pod_sensitivity <- function(parameters, solver, yr, co2_opt){
  pod_change_lower <- data.frame()
  pod_change_upper <- data.frame()
  
  for (arg in arg_names){
    lower_params <- parameters
    lower_params[[arg]] <- parameters[[arg]] * 0.9
    
    upper_params <- parameters
    upper_params[[arg]] <- parameters[[arg]] * 1.1
    
    result_lower_param <- solver(lower_params[arg_names]) 
    result_upper_param <- solver(upper_params[arg_names]) 
    
    lower_param_pod_change <- (max(result_lower_param$Pod)-max(result$Pod))/max(result$Pod)
    upper_param_pod_change <- (max(result_upper_param$Pod)-max(result$Pod))/max(result$Pod)
    
    # pod_change_lower[arg, paste0(yr, co2_opt, "-10%")] <- paste0(round(lower_param_pod_change * 100, 2), "%")
    # pod_change_upper[arg, paste0(yr, co2_opt, "+10%")] <- paste0(round(upper_param_pod_change * 100, 2), "%")
    pod_change_lower[arg, paste0(yr, co2_opt, "-10%")] <- lower_param_pod_change
    pod_change_upper[arg, paste0(yr, co2_opt, "+10%")] <- upper_param_pod_change
  }
  return(list(lower = pod_change_lower, upper = pod_change_upper))
}

result_sensitivity <- calculate_pod_sensitivity(parameters, soybean_optsolver[[i]], yr, co2_opt)
# Merge lower (-10%) results
if (is.null(pod_sensitivity_lower)) {
  pod_sensitivity_lower <- result_sensitivity$lower
} else {
  pod_sensitivity_lower <- cbind(pod_sensitivity_lower, result_sensitivity$lower)
}

# Merge upper (+10%) results
if (is.null(pod_sensitivity_upper)) {
  pod_sensitivity_upper <- result_sensitivity$upper
} else {
  pod_sensitivity_upper <- cbind(pod_sensitivity_upper, result_sensitivity$upper)
}
