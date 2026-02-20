set_direct_modules<-function(current_direct_modules){
  utr_direct_modules <- c("UTRSoybeanBML:thornley_utilization_calculator_lsrp",
                          "UTRSoybeanBML:thornley_transport_calculator_lsrp",
                          "UTRSoybeanBML:thornley_biomass_calculator_lsrp",
                          "UTRSoybeanBML:lai_from_structural_carbon")
  modules_to_remove <- c("BioCro:maintenance_respiration_calculator",
                         "BioCro:partitioning_coefficient_logistic",
                         "BioCro:partitioning_growth_calculator", 
                         "BioCro:parameter_calculator",
                         "BioCro:sla_linear")
  direct_modules_new <- current_direct_modules[!sapply(current_direct_modules, function(x) x %in% modules_to_remove)]
  direct_modules_new <- c(direct_modules_new, utr_direct_modules)
  return(direct_modules_new)
}

set_differential_modules<-function(current_differential_modules){
  utr_differential_modules <- c("UTRSoybeanBML:thornley_utilization_lsrp",
                                "UTRSoybeanBML:thornley_transport_lsrp")
  modules_to_remove <- c("BioCro:senescence_logistic",
                         "BioCro:maintenance_respiration",
                         "BioCro:partitioning_growth")#,
#                         "BioCro:development_index")
  differential_modules_new <- current_differential_modules[!sapply(current_differential_modules, function(x) x %in% modules_to_remove)]
  differential_modules_new <- c(differential_modules_new, utr_differential_modules)
  return(differential_modules_new)
}

set_init_values<-function(current_initial_values, exp_mass){
  sub_frac <- 0.1           # substrate_fraction
  str_frac <- 1 - sub_frac  # structural_fraction
  seed_mass <- 0.0789
  leaf_frac <- 0.45
  stem_frac <- 0.25
  root_frac <- 0.3
  cf <- 0.3 # optim_params_short_SoyFACE[1]
  
  utr_initial_values <- list(
    Leaf_respiration_loss = 0.0,
    Stem_respiration_loss = 0.0,
    Root_respiration_loss = 0.0,
    Pod_respiration_loss = 0.0,
    # senescence related initial state
    Leaf_senescence_loss = 0.0,
    Stem_senescence_loss = 0.0,
    Root_senescence_loss = 0.0,
    Pod_senescence_loss = 0.0,
    # Other variables
    # TTc = 0.0,
    # cws1=                     0.32,          # dimensionless, current water status, soil layer 1
    # cws2=                     0.32,          # dimensionless, current water status, soil layer 2
    # Rhizome =                 0.0000001,     # Mg / ha
    # RhizomeLitter =           0,               # Mg / ha
    Leaf_substrate_carbon = sub_frac * seed_mass * leaf_frac / cf,
    Leaf_structural_carbon = str_frac * seed_mass * leaf_frac / cf,
    Stem_substrate_carbon = sub_frac * seed_mass * stem_frac / cf,
    Stem_structural_carbon = str_frac * seed_mass * stem_frac / cf,
    Root_substrate_carbon =  sub_frac * seed_mass * root_frac / cf,
    Root_structural_carbon = str_frac * seed_mass * root_frac / cf,
    Pod_substrate_carbon = 1e-4,
    Pod_structural_carbon = 9e-4)
  
  initial_values_to_remove <- c("Leaf", "Stem", "Root", "Grain", "Shell")
  initial_values_new <- current_initial_values[!(names(current_initial_values) %in% initial_values_to_remove)]
  initial_values_new <- c(initial_values_new, utr_initial_values)
  initial_values_new$Sp
  return(initial_values_new)
}

set_parameters<-function(current_parameters){
  
  default_carbon_to_mass_factor <- 0.3
  base_utilization_rate_constant <- 3e-2
  base_utilization_km <- 1
  base_conductance <- 0.5
  
  utr_parameters <- list(
    # Parameters related to the UTR model
    Leaf_carbon_to_mass_factor = default_carbon_to_mass_factor,
    Leaf_utilization_rate_constant = 0.4 * base_utilization_rate_constant,
    Leaf_utilization_km = 1.0 * base_utilization_km,
    Leaf_respiration_factor = 0.0,				# Leaf respiration is accounted for
    # by the canopy photosynthesis module
    Stem_carbon_to_mass_factor = default_carbon_to_mass_factor,
    Stem_utilization_rate_constant = 1.0 * base_utilization_rate_constant,
    Stem_utilization_km = 1.2 * base_utilization_km,
    Stem_respiration_factor = 0.02,
    
    Root_carbon_to_mass_factor = default_carbon_to_mass_factor,
    Root_utilization_rate_constant = 0.1 * base_utilization_rate_constant,
    Root_utilization_km = 0.5 * base_utilization_km,
    Root_respiration_factor = 0.03,
    
    Pod_carbon_to_mass_factor = default_carbon_to_mass_factor,
    Pod_utilization_rate_constant = 1.0 * base_utilization_rate_constant,
    Pod_utilization_km = 1.0 * base_utilization_km,
    Pod_respiration_factor = 0.02,
    
    substrate_conductance_Leaf_to_Stem = base_conductance,
    substrate_conductance_Stem_to_Root = 0.4 * base_conductance,
    substrate_conductance_Stem_to_Pod = 2*base_conductance,
    
    transportation_beta_exponent = 1,
    Pod_start_dvi = 1.0,
    stop_growth_dvi = 2.0,
    
    # senescence_coefficient_logistic module
    Leaf_senescence_fraction_max =           0.003,
    Stem_senescence_fraction_max =           0.003,
    Root_senescence_fraction_max =           0.005,  # senescence of root not simulated in Soybean-BioCro (see Note 1 at end of file)
    Pod_senescence_fraction_max  =           0.0,
    Leaf_senescence_alpha    =               20.0,   
    Stem_senescence_alpha    =               20.0,
    Root_senescence_alpha    =               8.0,    # 3 for 2002,2005; 4 for 2004; 3 for 2006 (not sure)
    Pod_senescence_alpha     =               110.0,
    Leaf_senescence_beta     =               -0.1,
    Stem_senescence_beta     =               -0.1,
    Root_senescence_beta     =               -2.0,
    Pod_senescence_beta      =               -0.4,
    Leaf_senescence_reuse_factor =            0.5, # 0.6,
    Stem_senescence_reuse_factor =            0.5, # 0.6,
    Root_senescence_reuse_factor =            0.5, # 0.6,
    Pod_senescence_reuse_factor =             0 # 0.6,
)
  parameters <- c(current_parameters, utr_parameters)
  return(parameters)
}
