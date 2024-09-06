optim_params_short <-c(0.498562,    0.028291,    0.096382,    0.099309,    
                       0.707672,    0.203889,    0.119468,    0.043002,    
                       0.101489,    0.406375,    0.497366,    0.420480,    
                       3.735706,    0.008319,    0.010091,    0.005725,    
                       7.919072,    7.259589,    1.157723,    1.846692,    
                       1.999206,    1.773724,    0.854841,    1.983732) #0.854841
optim_params_short <-c(0.352141,    0.025470,    0.082996,    0.047930,    
                       0.922560,    0.401151,    0.402652,    0.271050,    
                       0.436600,    0.313086,    0.220120,    0.322480,    
                       1.944263,    0.009309,    0.003321,    0.001601,    
                       7.998886,    9.182912,    1.155008,    1.866540,    
                       1.838870,    1.884435,    0.924264,    2.038643) # 0.924264

optim_params_short <-c(0.270905,    0.084463,    0.072071,    0.070452,    
                       0.887315,    0.494679,    0.262958,    0.237797,    
                       0.072732,    0.328827,    0.494764,    0.387722,    
                       2.338996,    0.019494,    0.003122,    0.000927,    
                       9.943553,    9.625340,    1.873967,    1.957134,    
                       1.960444,    1.772205,    1.162893,    2.048275) # 8/17/24 reuse=0

optim_params_short_SoyFACE <-c(0.375802,    0.006943,    0.016696,    0.038165,    
                               0.462693,    0.147413,    0.166919,    0.275076,    
                               0.172022,    0.108513,    0.323036,    0.266830,    
                               1.512302,    0.054309,    0.001903,    0.000713,    
                               8.047166,    6.630698,    1.012040,    
                               1.993078,    1.901416,    1.993260,    
                               1.011228,    2.094762) # 2.094762 optim_params_short_SoyFACE

optim_params_short_SoyFACE <-c(0.401895,    0.004262,    0.042255,    0.036734,
                               0.885907,    0.207854,    0.316474,    0.279429,
                               0.329332,    0.100413,    0.055206,    0.324411,
                               2.342311,    0.084638,    0.000346,    0.000765,
                               9.106215,    9.580776,    1.072162,    1.964001,
                               1.553672,    1.826175,    0.998824,    2.194033) # with constaints

optim_params_short_SoyFACE <-c(0.397978,    0.004262,    0.042255,    0.036734,    
                               0.885907,    0.207854,    0.316474,    0.279429,    
                               0.329332,    0.100413,    0.055206,    0.324411,    
                               2.440188,    0.088113,    0.000346,    0.000765,    
                               9.106215,    9.580776,    1.072162,    1.964001,    
                               1.555754,    1.890465,    1.002523,    2.194033) # optim_params_short_SoyFACE


optim_params_short_SoyFACE <-c(0.262338,    0.032370,    0.072937,    0.096153,    
                               0.917784,    0.289141,    0.139883,    0.119164,    
                               0.046467,    0.539973,    0.366747,    0.472946,    
                               0.926233,    0.051395,    0.001523,    0.001463,    
                               9.107290,    9.807757,    1.709075,    1.989716,    
                               1.948408,    1.981774,    1.081241,    2.017340) # 9/1/24 with 0 reuse 

optim_params_short_SoyFACE <-c(0.262338,    0.032370,    0.072937,    0.096153,    
                               0.917784,    0.289141,    0.139883,    0.119164,    
                               0.046467,    0.539973,    0.366747,    0.472946,    
                               0.926233,    0.051395,    0.001523,    0.001463,    
                               9.107290,    9.807757,    1.709075,    1.989716,    
                               1.948408,    1.981774,    1.081241,    2.017340) # 9/3/24 with 0 reuse and increased Pod 





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

utr.model.params <- data.frame(arg_names_short, optim_params_short, optim_params_short_SoyFACE)
names(utr.model.params) <- c('Parameters', 'EF', 'SF')
utr.model.params$change <- 100*(utr.model.params$EF-utr.model.params$SF)/utr.model.params$SF
# utr.model.params[] <- lapply(utr.model.params, function(x) if(is.numeric(x)&&(x>0.01)) round(x, 2) else x)
View(utr.model.params)
write.csv(utr.model.params, 'UTR_params.csv')
