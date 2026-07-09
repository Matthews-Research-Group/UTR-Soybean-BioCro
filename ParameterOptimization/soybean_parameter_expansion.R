# a function that converts 20 parameters to 28 based on assumptions
optim_params_conversion <- function(optim_params_short){
  optim_params <- c()
  
  # Utilization
  ## utilization rate constant [/hr]
  optim_params[1] = optim_params_short[1] # Leaf
  optim_params[2] = optim_params_short[2] # Stem
  optim_params[3] = optim_params_short[3] # Root
  optim_params[4] = optim_params_short[4] # Pod

  ## Km [dimensionless]
  optim_params[5] = optim_params_short[5]
  optim_params[6] = optim_params_short[6]
  optim_params[7] = optim_params_short[7]
  optim_params[8] = optim_params_short[8]

  ## respiration factor [dimensionless]
  # Leaf respiration is accounted for by the canopy photosynthesis module
  optim_params[c(9, 10, 11)] = optim_params_short[9] # Stem, Root, Pod

  # Transport
  ## substrate conductance [/hr]
  optim_params[12] = optim_params_short[10]  # 'substrate_conductance_Leaf_to_Stem'
  optim_params[13] = optim_params_short[11] # 'substrate_conductance_Stem_to_Root'
  optim_params[14] = optim_params_short[12] # 'substrate_conductance_Stem_to_Pod'

  # Senescence
  ## senescence max rates [/hr]
  optim_params[15] = optim_params_short[13] # 'Leaf_senescence_rate_max'
  optim_params[16] = optim_params_short[14] # 'Stem_senescence_rate_max'
  optim_params[17] = optim_params_short[15] # 'Root_senescence_rate_max'
  ## senescence alphas [dimensionless]
  optim_params[18] = optim_params_short[16] # 'Leaf_senescence_alpha'
  optim_params[19] = optim_params_short[17] # 'Stem_senescence_alpha'
  optim_params[20] = optim_params_short[18] # 'Root_senescence_alpha'
  ## senescence betas [/dvi]
  optim_params[21] = optim_params_short[19] # 'Leaf_senescence_beta'
  optim_params[22] = optim_params_short[20] # 'Stem_senescence_beta'
  optim_params[23] = optim_params_short[21] # 'Root_senescence_beta'
  ## senescence reuse factor
  optim_params[c(24, 25, 26)] = optim_params_short[22]

  # DVI switches
  optim_params[27] = optim_params_short[23] # 'Pod_start_dvi'
  optim_params[28] = optim_params_short[24] # 'stop_growth_dvi'
  
  return(optim_params)
}
