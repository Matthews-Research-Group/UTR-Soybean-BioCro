optim_params_short <-c(0.456593,    0.074827,    0.057426,    0.098052,    
                       0.729261,    0.324457,    0.087761,    0.083432,    
                       0.106687,    0.461096,    0.494880,    0.497674,    
                       4.566297,    0.011640,    0.006235,    0.001773,    
                       9.993437,    7.074247,    1.863020,    
                       1.862072,    1.968604,    1.998921,    
                       0.804830,    2.011087)
# optim_params_short <-c(0.455638,    0.014044,    0.071346,    0.089825,    
#                        0.982832,    0.119305,    0.109690,    0.095072,    
#                        0.112983,    0.481014,    0.485386,    0.478083,    
#                        4.457549,    0.013613,    0.004995,    0.003037,    
#                        6.795329,    6.623093,    1.968357,    1.926170,    
#                        1.952613,    1.998567,    0.843298,    2.037519)

optim_params_short_SoyFACE <-c(0.375802,    0.006943,    0.016696,    0.038165,    
                               0.462693,    0.147413,    0.166919,    0.275076,    
                               0.172022,    0.108513,    0.323036,    0.266830,    
                               1.512302,    0.054309,    0.001903,    0.000713,    
                               8.047166,    6.630698,    1.012040,    
                               1.993078,    1.901416,    1.993260,    
                               1.011228,    2.094762)

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
