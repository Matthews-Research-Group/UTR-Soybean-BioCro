# Do the calculations inside an empty list so that temporary variables are not created in .Global.
# I have not figured out how to do the calculations inside an empty list while calling elements in a data frame yet...
cf <- 0.25   # [Mg/ha]/[mmol/m2] # from the optimized result
leaf_frac <- 0.45
stem_frac <- 0.25
root_frac <- 0.3
subs_frac <- 0.1
struc_frac <- 1 - subs_frac
# For the initial total seed mass per land area, we use the following equation:
# Number of seeds per meter * weight per seed * 1 / row spacing
#
# Number of seeds per meter = 20 (Morgan et al., 2004, https://doi.org/10.1104/pp.104.043968)
# weight per seed = 0.15 grams / seed (https://www.feedipedia.org/node/42, average of .12 to .18 grams)
# row spacing = 0.38 meters (Morgan et al., 2004)
#
# (20 seeds / meter) * (0.15 grams / seed) * (1 / 0.38 meter) = 7.89 g / m^2 = 0.0789 Mg / ha
# This value is used to determine the initial Leaf, Stem, and Root biomasses
seed_mass <- 0.0789
Leaf <- seed_mass * leaf_frac
Stem <- seed_mass * stem_frac
Root <- seed_mass * root_frac

initial_state <- list(
  
  Leaf_respiration_loss = 0.0,
  Stem_respiration_loss = 0.0,
  Root_respiration_loss = 0.0,
  Pod_respiration_loss  = 0.0,
  
  # senescence related initial state
  Leaf_senescence_loss  = 0.0,
  Stem_senescence_loss  = 0.0,
  Root_senescence_loss  = 0.0,
  Pod_senescence_loss   = 0.0,
  
  # Other variables
  TTc    = 0.0,
  soil_water_content = 0.32,
  
  cws1                  = 0.32,          # dimensionless, current water status, soil layer 1
  cws2                  = 0.32,          # dimensionless, current water status, soil layer 2
  # DVI =                     -1,             # Sowing date: DVI=-1
  # Soybean does not have a rhizome, so these variables will not be used but must be defined
  Rhizome               = 0.0000001,     # Mg / ha
  RhizomeLitter =           0.0,               # Mg / ha
  
  # Biomass
  # Leaf = seed_mass * leaf_frac,
  # Stem = seed_mass * stem_frac,
  # Root = seed_mass * root_frac,
  # Pod = 1e-3, 
  
  # Substrate and structural C
  Leaf_substrate_carbon = subs_frac * Leaf / cf,
  Leaf_structural_carbon = struc_frac * Leaf / cf,
  Stem_substrate_carbon = subs_frac *  Stem / cf,
  Stem_structural_carbon = struc_frac * Stem / cf,
  Root_substrate_carbon =  subs_frac * Root / cf,
  Root_structural_carbon = struc_frac * Root / cf,
  
  Pod_substrate_carbon = 0 ,
  Pod_structural_carbon = 9e-4 / cf   
)

