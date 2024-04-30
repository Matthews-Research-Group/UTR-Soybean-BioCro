optim_params_short <-c(0.489689,    0.095792,    0.090287,    0.087059,    
                       0.981135,    0.443254,    0.126132,    0.077829,    
                       0.104172,    0.491694,    0.472555,    0.461276,    
                       4.793529,    0.001352,    0.000901,    0.002863,    
                       9.413562,    6.123178,    1.895080,    
                       1.507053,    1.501336,    1.872893, 
                       1.0, 2.0)

optim_params_short_SoyFACE <-c(0.385976,    0.005841,    0.030903,    0.023331,    
                               0.662863,    0.314122,    0.425564,    0.366934,    
                               0.394770,    0.119610,    0.072072,    0.375569,    
                               0.637063,    0.020278,    0.002460,    0.000764,    
                               5.746050,    7.207536,    1.001118,    
                               1.961757,    1.985062,    1.891085, 
                               1.0, 2.0)

arg_names <- c('Leaf_carbon_to_mass_factor', 'Stem_carbon_to_mass_factor', # 1, 2 
               'Root_carbon_to_mass_factor', 'Pod_carbon_to_mass_factor',  # 3, 4
               'Leaf_utilization_rate_constant', 'Stem_utilization_rate_constant', # 5, 6
               'Root_utilization_rate_constant', 'Pod_utilization_rate_constant',  # 7, 8
               'Leaf_utilization_km', 'Stem_utilization_km', 'Root_utilization_km', 'Pod_utilization_km', # 9, 10, 11, 12
               'Stem_respiration_factor', # 13
               'Root_respiration_factor', 'Pod_respiration_factor', # 14, 15
               'substrate_conductance_Leaf_to_Stem', 'substrate_conductance_Stem_to_Root', # 16, 17
               'substrate_conductance_Stem_to_Pod', #  'transportation_beta_exponent', # 18
               'Leaf_senescence_rate_max','Stem_senescence_rate_max', 'Root_senescence_rate_max', # 19, 20，21
               'Leaf_senescence_alpha', 'Stem_senescence_alpha', 'Root_senescence_alpha',# 22，23, 24,
               'Leaf_senescence_beta', 'Stem_senescence_beta', 'Root_senescence_beta', # 25, 26, 27
               'Pod_start_dvi', 'stop_growth_dvi')# 28, 29

arg_names_short <- c('carbon_to_mass_factor',# 1,
               'Leaf_utilization_rate_constant', 'Stem_utilization_rate_constant', # 2, 3
               'Root_utilization_rate_constant', 'Pod_utilization_rate_constant',  # 4, 5
               'Leaf_utilization_km', 'Stem_utilization_km', 'Root_utilization_km', 'Pod_utilization_km', # 5,7,8,9
               'respiration_factor', # 10
               'substrate_conductance_Leaf_to_Stem', 'substrate_conductance_Stem_to_Root', # 11, 12
               'substrate_conductance_Stem_to_Pod', #  'transportation_beta_exponent', # 13
               'Leaf_senescence_rate_max','Stem_senescence_rate_max', 'Root_senescence_rate_max', # 14, 15，16
               'Leaf_senescence_alpha', 'Stem_senescence_alpha', 'Root_senescence_alpha',# 17，18, 19,
               'Leaf_senescence_beta', 'Stem_senescence_beta', 'Root_senescence_beta', # 20, 21, 22,
               'Pod_start_dvi', 'stop_growth_dvi')# 23, 24

# utr.model.params <- data.frame(arg_names_short, optim_params_short, optim_params_short_SoyFACE)

# View(utr.model.params)
