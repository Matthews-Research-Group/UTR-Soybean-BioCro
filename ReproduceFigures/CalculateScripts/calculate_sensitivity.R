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

# Generate .tex file for the manuscript
generate_latex_table <- function(pod_sensitivity_lower, pod_sensitivity_upper){
  print(pod_sensitivity_lower)
  print(pod_sensitivity_upper)
  
  # ---------------------------------------------------------------------
  # Build sensitivity_table.tex: for each UTR parameter, the min/max
  # percent change in final pod mass across all growing seasons (ambient,
  # elevated CO2, and LD11-2170) when that parameter is perturbed -10%
  # and +10%. Row order follows arg_names.
  # ---------------------------------------------------------------------
  sensitivity_labels <- data.frame(
    arg = arg_names,
    parameter = c(
      'Max utilization rate (leaf)', 
      'Max utilization rate (stem)',
      'Max utilization rate (root)', 
      'Max utilization rate (pod)',
      'Dissociation constant (leaf)', 
      'Dissociation constant (stem)',
      'Dissociation constant (root)', 
      'Dissociation constant (pod)',
      'Respiration factor (stem)', 
      'Respiration factor (root)', 
      'Respiration factor (pod)',
      'Transport conductance from leaf to stem', 
      'Transport conductance from stem to root',
      'Transport conductance stem to pod',
      'Senescence rate max (leaf)', 
      'Senescence rate max (stem)', 
      'Senescence rate max (root)',
      'Senescence change rate with DVI (leaf)', 
      'Senescence change rate with DVI (stem)',
      'Senescence change rate with DVI (root)',
      'Senescence DVI timing (leaf)', 
      'Senescence DVI timing (stem)', 
      'Senescence DVI timing (root)',
      'Retained fraction (leaf)', 
      'Retained fraction (stem)', 
      'Retained fraction (root)',
      'DVI threshold for start of pod growth', 
      'DVI threshold for end of crop growth'
    ),
    symbol = c(
      '$v_{\\max,L}$', 
      '$v_{\\max,S}$', 
      '$v_{\\max,R}$', 
      '$v_{\\max,P}$',
      '$K_{L}$', 
      '$K_{S}$', 
      '$K_{R}$', 
      '$K_{P}$',
      '$k_{res,S}$', 
      '$k_{res,R}$', 
      '$k_{res,P}$',
      '$\\sigma_{L \\to S}$', 
      '$\\sigma_{S \\to R}$', 
      '$\\sigma_{S \\to P}$',
      '$v_{sene,\\max,L}$', 
      '$v_{sene,\\max,S}$', 
      '$v_{sene,\\max,R}$',
      '$\\alpha_{L}$', 
      '$\\alpha_{S}$', 
      '$\\alpha_{R}$',
      '$\\beta_{L}$', 
      '$\\beta_{S}$', 
      '$\\beta_{R}$',
      '$f_{r,L}$', 
      '$f_{r,S}$', 
      '$f_{r,R}$',
      '$\\text{DVI}_\\text{pod start}$', 
      '$\\text{DVI}_\\text{stop growth}$'
    ),
    stringsAsFactors = FALSE
  )
  
  lower_pct <- pod_sensitivity_lower[sensitivity_labels$arg, , drop = FALSE] * 100
  upper_pct <- pod_sensitivity_upper[sensitivity_labels$arg, , drop = FALSE] * 100
  
  sensitivity_labels$lower_min <- apply(lower_pct, 1, min)
  sensitivity_labels$lower_max <- apply(lower_pct, 1, max)
  sensitivity_labels$upper_min <- apply(upper_pct, 1, min)
  sensitivity_labels$upper_max <- apply(upper_pct, 1, max)
  
  fmt <- function(x) sprintf('%.2f', x)
  
  table_body <- paste0(
    '  ', sensitivity_labels$parameter, ' & ', sensitivity_labels$symbol, ' & ',
    fmt(sensitivity_labels$lower_min), ' & ', fmt(sensitivity_labels$lower_max), ' & ',
    fmt(sensitivity_labels$upper_min), ' & ', fmt(sensitivity_labels$upper_max), ' \\\\'
  )
  
  sensitivity_table_tex <- c(
    '\\centering',
    '\\begin{tabular}{llrrrr}',
    '  \\toprule',
    '  & & \\multicolumn{2}{c}{$-10\\%$} & \\multicolumn{2}{c}{$+10\\%$} \\\\',
    '  \\cmidrule(lr){3-4} \\cmidrule(lr){5-6}',
    '  Parameter & Symbol & Min & Max & Min & Max \\\\',
    '  \\midrule',
    table_body,
    '  \\bottomrule',
    '\\end{tabular}'
  )
  
  writeLines(sensitivity_table_tex, 'sensitivity_table.tex')
}
