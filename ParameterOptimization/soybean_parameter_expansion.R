# a function that converts 20 parameters to 28 based on assumptions
optim_params_conversion <- function(optim_params_short){
  optim_params <- c()
  
  # carbon to mass factor [(Mg / ha) / (mol / m^2)]
  optim_params[c(1,2,3,4)] = optim_params_short[1]
  
  # Utilization
  ## utilization rate constant [/hr]
  optim_params[5] = optim_params_short[2] # Leaf
  optim_params[6] = optim_params_short[3] # Stem
  optim_params[7] = optim_params_short[4] # Root
  optim_params[8] = optim_params_short[5] # Pod
  
  ## Km [mol / Mg]
  optim_params[9] = optim_params_short[6]
  optim_params[10] = optim_params_short[7]
  optim_params[11] = optim_params_short[8]
  optim_params[12] = optim_params_short[9]
  
  ## respiration factor [dimensionless]
  optim_params[13] = 0 # Leaf respiration is accounted for by the canopy photosynthesis module
  optim_params[c(14, 15, 16)] = optim_params_short[10]
  
  # Transport
  ## substrate conductance [Mg / hr / [Mg / ha]]
  optim_params[17] = optim_params_short[11]  # 'substrate_conductance_Leaf_to_Stem'
  optim_params[18] = optim_params_short[12] # 'substrate_conductance_Stem_to_Root'
  optim_params[19] = optim_params_short[13] # 'substrate_conductance_Stem_to_Pod'
  
  # Senescence
  ## senescence max rates [/hr]
  optim_params[20] = optim_params_short[14] # 'Leaf_senescence_rate_max'
  optim_params[21] = optim_params_short[15] # 'Stem_senescence_rate_max'
  optim_params[22] = optim_params_short[16] # 'Root_senescence_rate_max'
  ## senescence alphas [dimensionless]
  optim_params[23] = optim_params_short[17] # 'Leaf_senescence_alpha'
  optim_params[24] = optim_params_short[18] # 'Stem_senescence_alpha'
  optim_params[25] = optim_params_short[19] # 'Root_senescence_alpha'
  ## senescence betas [/dvi]
  optim_params[26] = optim_params_short[20] # 'Leaf_senescence_beta'
  optim_params[27] = optim_params_short[21] # 'Stem_senescence_beta'
  optim_params[28] = optim_params_short[22] # 'Root_senescence_beta'
  
  return(optim_params)
}
