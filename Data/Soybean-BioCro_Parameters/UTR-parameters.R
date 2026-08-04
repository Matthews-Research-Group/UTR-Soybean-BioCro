# Read the file
lines <- readLines("../ParameterOptimization/Optmization_output_2026-08-04.txt") # 07-29

optim_params_short_SoyFACE <- as.numeric(strsplit(trimws(sub(".*bestmemit:", "", lines[1000])), "\\s+")[[1]])

arg_names <- c('Leaf_utilization_rate_constant', 'Stem_utilization_rate_constant', 
               'Root_utilization_rate_constant', 'Pod_utilization_rate_constant',
               'Leaf_utilization_km', 'Stem_utilization_km', 'Root_utilization_km', 'Pod_utilization_km',
               'Stem_respiration_factor',
               'Root_respiration_factor', 'Pod_respiration_factor',
               'substrate_conductance_Leaf_to_Stem', 'substrate_conductance_Stem_to_Root',
               'substrate_conductance_Stem_to_Pod',
               'Leaf_senescence_fraction_max','Stem_senescence_fraction_max', 'Root_senescence_fraction_max', 
               'Leaf_senescence_alpha', 'Stem_senescence_alpha', 'Root_senescence_alpha',
               'Leaf_senescence_beta', 'Stem_senescence_beta', 'Root_senescence_beta', 
               'Leaf_senescence_reuse_factor', 'Stem_senescence_reuse_factor', 'Root_senescence_reuse_factor', 
               'Pod_start_dvi', 'stop_growth_dvi') 

arg_names_short <- c(
               'Leaf_utilization_rate_constant', 'Stem_utilization_rate_constant',
               'Root_utilization_rate_constant', 'Pod_utilization_rate_constant',
               'Leaf_utilization_km', 'Stem_utilization_km', 'Root_utilization_km', 'Pod_utilization_km', 
               'respiration_factor',
               'substrate_conductance_Leaf_to_Stem', 'substrate_conductance_Stem_to_Root', 
               'substrate_conductance_Stem_to_Pod', 
               'Leaf_senescence_fraction_max','Stem_senescence_fraction_max', 'Root_senescence_fraction_max',
               'Leaf_senescence_alpha', 'Stem_senescence_alpha', 'Root_senescence_alpha',
               'Leaf_senescence_beta', 'Stem_senescence_beta', 'Root_senescence_beta', 
               'senescence_reuse_factor',
               'Pod_start_dvi', 'stop_growth_dvi')


