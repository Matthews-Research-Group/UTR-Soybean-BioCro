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

# =============================================================================
# Build utr-parameters.tex: parameter names, symbols, units, estimation
# ranges, and fitted values (optim_params_short_SoyFACE), in the order of
# arg_names_short.
#
# param_lower / param_upper are the same DEoptim search bounds used in
# ../../ParameterOptimization/multiyear_utr_params_optimization.R
# (kept in sync manually since that script sources this one, so it cannot be
# sourced back here without creating a circular dependency).
# =============================================================================
param_lower <- c(
  0.0, 0.0, 0.0, 0.0,   # utilization rate constant [/hr]
  0.0, 0.0, 0.0, 0.0,   # Km [/]
  0.2,                  # respiration factor [/]
  0.0, 0.0, 0.0,        # substrate conductance [Mg / hr / [Mg / ha]^beta]
  0.0, 0.0, 0.0,        # senescence rate max, LSR
  0.0, 0.0, 0.0,        # senescence alpha, LSR [dimensionless]
  1.5, 1.5, 1.5,        # senescence beta, LSR [DVI^-1]
  0.0,                  # senescence reuse factor
  1.0, 1.8)             # DVI switches
param_upper <- c(
  0.1, 0.1, 0.1, 1.0,   # utilization rate constant [/hr]
  0.5, 0.5, 0.5, 0.5,   # Km [/]
  0.8,                  # respiration factor [/]
  1.0, 1.0, 5.0,        # substrate conductance [Mg / hr / [Mg / ha]^beta]
  0.1, 0.1, 0.01,       # senescence rate max, LSR
  10.0, 10.0, 2.0,      # senescence alpha, LSR [dimensionless]
  2.5, 2.5, 2.5,        # senescence beta, LSR [DVI^-1]
  1.0,                  # senescence reuse factor
  1.3, 2.2)             # DVI switches

param_info <- data.frame(
  section = c(rep('Utilization Parameters', 9),
              rep('Transport Parameters', 3),
              rep('Senescence Parameters', 10),
              rep('DVI switches', 2)),
  parameter = c(
    'Max utilization rate constant (L)', 'Max utilization rate constant (S)',
    'Max utilization rate constant (R)', 'Max utilization rate constant (P)',
    'Half-saturation constant (L)', 'Half-saturation constant (S)',
    'Half-saturation constant (R)', 'Half-saturation constant (P)',
    'Respiration fraction (SRP)',
    'Transport conductance from leaf to stem', 'Transport conductance from stem to root',
    'Transport conductance from stem to pod',
    'Max senescence fraction rate (L)', 'Max senescence fraction rate (S)', 'Max senescence fraction rate (R)',
    'Senescence transition steepness with DVI (L)', 'Senescence transition steepness with DVI (S)',
    'Senescence transition steepness with DVI (R)',
    'Senescence midpoint DVI (L)', 'Senescence midpoint DVI (S)', 'Senescence midpoint DVI (R)',
    'Retained substrate fraction (LSR)',
    'DVI threshold for start of pod growth', 'DVI threshold for end of crop growth'
  ),
  symbol = c(
    '$v_{\\max,L}$', '$v_{\\max,S}$', '$v_{\\max,R}$', '$v_{\\max,P}$',
    '$K_{L}$', '$K_{S}$', '$K_{R}$', '$K_{P}$',
    '$k_{res,O}$',
    '$\\sigma_{L \\to S}$', '$\\sigma_{S \\to R}$', '$\\sigma_{S \\to P}$',
    '$v_{sene,\\max,L}$', '$v_{sene,\\max,S}$', '$v_{sene,\\max,R}$',
    '$\\alpha_{L}$', '$\\alpha_{S}$', '$\\alpha_{R}$',
    '$\\beta_{L}$', '$\\beta_{S}$', '$\\beta_{R}$',
    '$f_{r}$',
    '$\\text{DVI}_\\text{pod start}$', '$\\text{DVI}_\\text{stop growth}$'
  ),
  unit = c(
    'hr$^{-1}$', 'hr$^{-1}$', 'hr$^{-1}$', 'hr$^{-1}$',
    '--', '--', '--', '--',
    '--',
    'hr$^{-1}$', 'hr$^{-1}$', 'hr$^{-1}$',
    'hr$^{-1}$', 'hr$^{-1}$', 'hr$^{-1}$',
    '--', '--', '--',
    'DVI', 'DVI', 'DVI',
    '--',
    '--', '--'
  ),
  min = param_lower,
  max = param_upper,
  eq = c(rep('(\\ref{eq:utilizationrate})', 8),
         '(\\ref{eq:respirationrate})',
         rep('(\\ref{eq:transportrate})', 3),
         rep('(\\ref{eq:senesce_fraction})', 9),
         '(\\ref{eq:substratesenescencerate})(\\ref{eq:retainsenescencerate})',
         '--', '--'),
  stringsAsFactors = FALSE
)

param_info$fit <- sprintf('%.6f', optim_params_short_SoyFACE[seq_along(arg_names_short)])

table_rows <- character(0)
for (i in seq_len(nrow(param_info))) {
  if (i == 1 || param_info$section[i] != param_info$section[i - 1]) {
    table_rows <- c(table_rows,
                     paste0('  \\multicolumn{7}{l}{\\textbf{', param_info$section[i], '}} \\\\'))
  }
  table_rows <- c(table_rows, paste0(
    '  ', param_info$parameter[i], ' & ', param_info$symbol[i], ' & ', param_info$unit[i],
    ' & ', param_info$min[i], ' & ', param_info$max[i], ' & ', param_info$fit[i],
    ' & ', param_info$eq[i], ' \\\\'))
}

utr_parameters_tex <- c(
  '\\begin{table}[H]',
  '\\caption{UTR model parameters, estimation ranges, parameterization results, and corresponding equations. L: Leaf. S: Stem. R: Root. P: Pod}',
  '\\label{tab:utr-param-range-value}',
  '\\centering',
  '\\begin{tabular}{lllrrrc}',
  '  \\toprule',
  '  Parameter & Symbol & Unit & Min & Max & Fit & Eq \\\\',
  '  \\midrule',
  table_rows,
  '  \\bottomrule',
  '\\end{tabular}',
  '\\end{table}'
)

writeLines(utr_parameters_tex, '../Data/Soybean-BioCro_Parameters/utr-parameters.tex')


