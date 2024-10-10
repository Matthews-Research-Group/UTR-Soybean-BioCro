optim_params_short <-c(# 0.300294,    
                       0.040210,    0.069855,    0.070787,    0.608021,    
                       0.380683,    0.289376,    0.215143,    0.220713,    0.323682,    
                       0.468640,    0.228597,    2.774941,    0.015044,    0.022064,    
                       0.041733,    9.768698,    6.779807,    1.995181,    1.534624,    
                       1.888500,    1.998240,    0.804800,    1.017007,    2.049954)

optim_params_short <-c(# 0.308663,    
                       0.047340,    0.038147,    0.093637,    0.799633,    
                       0.378279,    0.149303,    0.159667,    0.204075,    0.350841,    
                       0.385811,    0.239362,    2.711464,    0.020813,    0.001882,    
                       0.002422,    9.799749,    7.380839,    1.183634,    1.986812,    
                       1.765750,    1.985843,    0.752978,    0.987921,    2.050715) 
optim_params_short <-c(# 0.345367,    
                       0.002620,    0.034941,    0.058760,    0.312415,    
                       0.143702,    0.289450,    0.317765,    0.134854,    0.103661,    
                       0.048093,    0.319665,    1.506278,    0.024355,    0.001974,    
                       0.006506,    9.108652,    4.088730,    1.659319,    1.843764,    
                       1.959387,    1.984124,    0.999885,    1.153776,    2.010368) # 10/2/2024

optim_params_short_SoyFACE <-c(# 0.316160,    
                               0.005007,    0.013080,    0.095773,    0.750602,    
                               0.187277,    0.091740,    0.245969,    0.118309,    0.144260,    
                               0.082548,    0.487179,    4.378294,    0.027548,    0.000772,    
                               0.003272,    5.868435,    5.329164,    1.049071,    1.993587,    
                               1.584600,    1.896079,    0.996458,    1.157210,    2.062570)


optim_params_short_SoyFACE <-c(# 0.340787,    
                               0.006314,    0.026195,    0.071335,    0.293424,    
                               0.158471,    0.177561,    0.238038,    0.032516,    0.344832,    
                               0.278962,    0.247008,    1.274300,    0.050121,    0.003243,    
                               0.001230,    9.963829,    7.405092,    0.258005,    1.965341,    
                               1.986972,    1.539534,    0.999618,    1.146891,    1.986127)

optim_params_short_SoyFACE <-c(# 0.270165,    
                               0.010248,    0.054891,    0.067387,    0.583982,    
                               0.286407,    0.445138,    0.474320,    0.304828,    0.126300,    
                               0.233157,    0.472530,    1.757361,    0.039268,    0.001498,    
                               0.003244,    9.331028,    5.771867,    1.808724,    1.967422,    
                               1.769574,    1.988078,    0.999927,    1.093697,    1.989707) # 1.093697
optim_params_short_SoyFACE <-c(0.345367,    
                               0.002620,    0.034941,    0.058760,    0.312415,    
                               0.143702,    0.289450,    0.317765,    0.134854,    0.103661,    
                               0.048093,    0.319665,    1.506278,    0.024355,    0.001974,    
                               0.006506,    9.108652,    4.088730,    1.659319,    1.843764,    
                               1.959387,    1.984124,    0.999885,    1.153776,    2.010368) # 10/2/2024 1.153776

optim_params_short_SoyFACE <-c(0.003746,    0.017267,    0.063001,    0.187819,
                               0.377536,    0.363415,    0.040322,    0.169470,    0.102317,
                               0.033388,    0.020354,    2.845359,    0.023626,    0.003685,
                               0.003807,    9.868176,    9.075340,    1.652498,    1.848524,
                               1.991700,    1.984780,    0.998544,    1.172285,    2.000373) # 10/7/2024


arg_names <- c(# 'Leaf_carbon_to_mass_factor', 'Stem_carbon_to_mass_factor', # 1, 2 
               # 'Root_carbon_to_mass_factor', 'Pod_carbon_to_mass_factor',  # 3, 4
               'Leaf_utilization_rate_constant', 'Stem_utilization_rate_constant', # 5, 6
               'Root_utilization_rate_constant', 'Pod_utilization_rate_constant',  # 7, 8
               'Leaf_utilization_km', 'Stem_utilization_km', 'Root_utilization_km', 'Pod_utilization_km', # 9, 10, 11, 12
               'Stem_respiration_factor', # 13
               'Root_respiration_factor', 'Pod_respiration_factor', # 14, 15
               'substrate_conductance_Leaf_to_Stem', 'substrate_conductance_Stem_to_Root', # 16, 17
               'substrate_conductance_Stem_to_Pod', #  'transportation_beta_exponent', # 18
               'Leaf_senescence_fraction_max','Stem_senescence_fraction_max', 'Root_senescence_fraction_max', # 19, 20，21
               'Leaf_senescence_alpha', 'Stem_senescence_alpha', 'Root_senescence_alpha',# 22，23, 24,
               'Leaf_senescence_beta', 'Stem_senescence_beta', 'Root_senescence_beta', # 25, 26, 27
               'Leaf_senescence_reuse_factor', 'Stem_senescence_reuse_factor', 'Root_senescence_reuse_factor', # 28, 29, 30
               'Pod_start_dvi', 'stop_growth_dvi') # 31, 32

arg_names_short <- c(# 'carbon_to_mass_factor',# 1,
               'Leaf_utilization_rate_constant', 'Stem_utilization_rate_constant', # 2, 3
               'Root_utilization_rate_constant', 'Pod_utilization_rate_constant',  # 4, 5
               'Leaf_utilization_km', 'Stem_utilization_km', 'Root_utilization_km', 'Pod_utilization_km', # 5,7,8,9
               'respiration_factor', # 10
               'substrate_conductance_Leaf_to_Stem', 'substrate_conductance_Stem_to_Root', # 11, 12
               'substrate_conductance_Stem_to_Pod', #  'transportation_beta_exponent', # 13
               'Leaf_senescence_fraction_max','Stem_senescence_fraction_max', 'Root_senescence_fraction_max',# 14, 15，16
               'Leaf_senescence_alpha', 'Stem_senescence_alpha', 'Root_senescence_alpha',# 17，18, 19,
               'Leaf_senescence_beta', 'Stem_senescence_beta', 'Root_senescence_beta', # 20, 21, 22,
               'senescence_reuse_factor', # 23
               'Pod_start_dvi', 'stop_growth_dvi')# 24, 25

# utr.model.params <- data.frame(arg_names_short, optim_params_short, optim_params_short_SoyFACE)
# names(utr.model.params) <- c('Parameters', 'EF', 'SF')
# utr.model.params$change <- 100*(utr.model.params$EF-utr.model.params$SF)/utr.model.params$SF
# utr.model.params[] <- lapply(utr.model.params, function(x) {
#   if(is.numeric(x)) { sapply(x, function(y) if(y > 0.01) round(y, 2) else y)
#   } else {x}
# })
# View(utr.model.params)
# write.csv(utr.model.params, 'UTR_params.csv')
