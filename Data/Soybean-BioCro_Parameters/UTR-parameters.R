optim_params_short <-c(0.388810,    0.066280,    0.089066,    0.021532,
                       0.712018,    0.232411,    0.093201,    0.011088,
                       0.107644,    0.305144,    0.498137,    0.406344,
                       4.426481,    0.015629,    0.007029,    0.007164,
                       6.720799,    6.881872,    1.499185,    1.948120,
                       1.908756,    1.929277,    0.745088,    1.980669)
optim_params_short_SoyFACE <-c(0.375802,    0.006943,    0.016696,    0.038165,    
                               0.462693,    0.147413,    0.166919,    0.275076,    
                               0.172022,    0.108513,    0.323036,    0.266830,    
                               1.512302,    0.054309,    0.001903,    0.000713,    
                               8.047166,    6.630698,    1.012040,    
                               1.993078,    1.901416,    1.993260,    
                               1.011228,    2.094762) # optim_params_short_SoyFACE

# optim_params_short_SoyFACE <-c(0.401895,    0.004262,    0.042255,    0.036734,
#                                0.885907,    0.207854,    0.316474,    0.279429,
#                                0.329332,    0.100413,    0.055206,    0.324411,
#                                2.342311,    0.084638,    0.000346,    0.000765,
#                                9.106215,    9.580776,    1.072162,    1.964001,
#                                1.553672,    1.826175,    0.998824,    2.194033) # with constaints
# 
# 
# optim_params_short_SoyFACE <-c(0.226708,    0.071574,    0.021779,    0.055257,    
#                                0.672280,    0.166042,    0.023193,    0.009908,    
#                                0.386762,    0.106646,    0.496009,    0.470880,    
#                                0.866362,    0.062559,    0.012655,    0.004552,    
#                                9.989777,    8.145810,    1.702138,    1.955685,    
#                                1.996797,    1.928436,    0.919358,    1.949513) # 6-6

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
