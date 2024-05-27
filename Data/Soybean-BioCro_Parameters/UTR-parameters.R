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

optim_params_short <-c(0.483128,    0.009474,    0.093663,    0.053287,    
                       0.992664,    0.060729,    0.096953,    0.006831,    
                       0.118792,    0.379283,    0.479401,    0.307510,    
                       4.905103,    0.010397,    0.007886,    0.005794,    
                       6.015208,    9.331500,    1.523546,    1.833710,    
                       1.919567,    1.790103,    0.805621,    1.972447)


optim_params_short <-c(0.438147,    0.019763,    0.079460,    0.090002,    
                       0.826861,    0.134472,    0.122367,    0.078384,    
                       0.107164,    0.351642,    0.469405,    0.368380,    
                       4.988032,    0.013982,    0.011066,    0.002323,    
                       5.974400,    9.488358,    1.612573,    1.929992,    
                       1.959475,    1.824268,    0.966278,    1.948314)
optim_params_short <-c(0.440478,    0.061241,    0.087037,    0.095199,    
                       0.792503,    0.254855,    0.092947,    0.050497,    
                       0.116680,    0.461427,    0.486983,    0.444954,    
                       2.672736,    0.027970,    0.004377,    0.003015,    
                       8.929532,    7.270175,    1.846408,    1.961543,    
                       1.825278,    1.991865,    0.623041,    1.983941)

optim_params_short <-c(0.454093,    0.081316,    0.072802,    0.033383,    
                       0.936266,    0.305164,    0.073138,    0.014884,    
                       0.102959,    0.431183,    0.495048,    0.499982,    
                       2.725879,    0.019159,    0.007548,    0.006826,    
                       9.470688,    5.767104,    1.830186,    1.918061,    
                       1.999822,    1.999838,    0.602157,    1.996265) #0.642157

optim_params_short_SoyFACE <-c(0.375802,    0.006943,    0.016696,    0.038165,    
                               0.462693,    0.147413,    0.166919,    0.275076,    
                               0.172022,    0.108513,    0.323036,    0.266830,    
                               1.512302,    0.054309,    0.001903,    0.000713,    
                               8.047166,    6.630698,    1.012040,    
                               1.993078,    1.901416,    1.993260,    
                               1.011228,    2.094762) # optim_params_short_SoyFACE

optim_params_short_SoyFACE <-c(0.401895,    0.004262,    0.042255,    0.036734,    
                               0.885907,    0.207854,    0.316474,    0.279429,    
                               0.329332,    0.100413,    0.055206,    0.324411,    
                               2.342311,    0.084638,    0.000346,    0.000765,    
                               9.106215,    9.580776,    1.072162,    1.964001,    
                               1.553672,    1.826175,    0.998824,    2.194033)


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
