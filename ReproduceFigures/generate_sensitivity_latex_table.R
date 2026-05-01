summary_table <- data.frame(
  Lower_Min = round(apply(pod_sensitivity_lower, 1, min, na.rm = TRUE) * 100, 2),
  Lower_Max = round(apply(pod_sensitivity_lower, 1, max, na.rm = TRUE) * 100, 2),
  Upper_Min = round(apply(pod_sensitivity_upper, 1, min, na.rm = TRUE) * 100, 2),
  Upper_Max = round(apply(pod_sensitivity_upper, 1, max, na.rm = TRUE) * 100, 2)
)

colnames(summary_table) <- c("Min (-10%)", "Max (-10%)", "Min (+10%)", "Max (+10%)")

# Export to CSV
write.csv(summary_table, "pod_sensitivity_summary.csv", row.names = TRUE)

# Define parameter names and symbols matching your row order
param_info <- data.frame(
  Parameter = c(
    "Max utilization rate (leaf)",
    "Max utilization rate (stem)", 
    "Max utilization rate (root)",
    "Max utilization rate (pod)",
    "Michaelis-Menten constant (leaf)",
    "Michaelis-Menten constant (stem)",
    "Michaelis-Menten constant (root)",
    "Michaelis-Menten constant (pod)",
    "Respiration factor (stem)",
    "Respiration factor (root)",
    "Respiration factor (pod)",
    "Transport conductance from leaf to stem",
    "Transport conductance from stem to root",
    "Transport conductance stem to pod",
    "Senescence rate max (leaf)",
    "Senescence rate max (stem)",
    "Senescence rate max (root)",
    "Senescence change rate with DVI (leaf)",
    "Senescence change rate with DVI (stem)",
    "Senescence change rate with DVI (root)",
    "Senescence DVI timing (leaf)",
    "Senescence DVI timing (stem)",
    "Senescence DVI timing (root)",
    "Reused fraction (leaf)",
    "Reused fraction (stem)",
    "Reused fraction (root)",
    "DVI threshold for start of pod growth",
    "DVI threshold for end of crop growth"
  ),
  Symbol = c(
    "$r_{\\max,L}$",
    "$r_{\\max,S}$",
    "$r_{\\max,R}$",
    "$r_{\\max,P}$",
    "$K_{L}$",
    "$K_{S}$",
    "$K_{R}$",
    "$K_{P}$",
    "$k_{res,S}$",
    "$k_{res,R}$",
    "$k_{res,P}$",
    "$\\sigma_{L \\to S}$",
    "$\\sigma_{S \\to R}$",
    "$\\sigma_{S \\to P}$",
    "$r_{sene,\\max,L}$",
    "$r_{sene,\\max,S}$",
    "$r_{sene,\\max,R}$",
    "$\\alpha_{L}$",
    "$\\alpha_{S}$",
    "$\\alpha_{R}$",
    "$\\beta_{L}$",
    "$\\beta_{S}$",
    "$\\beta_{R}$",
    "$f_{r,L}$",
    "$f_{r,S}$",
    "$f_{r,R}$",
    "$DVI_{pod~start}$",
    "$DVI_{stop~growth}$"
  ),
  row.names = c(
    "Leaf_utilization_rate_constant",
    "Stem_utilization_rate_constant",
    "Root_utilization_rate_constant",
    "Pod_utilization_rate_constant",
    "Leaf_utilization_km",
    "Stem_utilization_km",
    "Root_utilization_km",
    "Pod_utilization_km",
    "Stem_respiration_factor",
    "Root_respiration_factor",
    "Pod_respiration_factor",
    "substrate_conductance_Leaf_to_Stem",
    "substrate_conductance_Stem_to_Root",
    "substrate_conductance_Stem_to_Pod",
    "Leaf_senescence_fraction_max",
    "Stem_senescence_fraction_max",
    "Root_senescence_fraction_max",
    "Leaf_senescence_alpha",
    "Stem_senescence_alpha",
    "Root_senescence_alpha",
    "Leaf_senescence_beta",
    "Stem_senescence_beta",
    "Root_senescence_beta",
    "Leaf_senescence_reuse_factor",
    "Stem_senescence_reuse_factor",
    "Root_senescence_reuse_factor",
    "Pod_start_dvi",
    "stop_growth_dvi"
  )
)

# Build the table body rows
rows <- mapply(function(param, symbol, lmin, lmax, umin, umax) {
  sprintf("  %s & %s & %.2f & %.2f & %.2f & %.2f \\\\",
          param, symbol, lmin, lmax, umin, umax)
},
param_info[rownames(summary_table), "Parameter"],
param_info[rownames(summary_table), "Symbol"],
summary_table[, "Min (-10%)"],
summary_table[, "Max (-10%)"],
summary_table[, "Min (+10%)"],
summary_table[, "Max (+10%)"]
)

latex_lines <- c(
  "\\begin{tabular}{llrrrr}",
  "  \\toprule",
  # First header row: grouped labels spanning Min/Max pairs
  "  & & \\multicolumn{2}{c}{$-10\\%$} & \\multicolumn{2}{c}{$+10\\%$} \\\\",
  "  \\cmidrule(lr){3-4} \\cmidrule(lr){5-6}",
  # Second header row: Min / Max labels
  "  Parameter & Symbol & Min & Max & Min & Max \\\\",
  "  \\midrule",
  rows,
  "  \\bottomrule",
  "\\end{tabular}"
)

writeLines(latex_lines, "tables/sensitivity_table.tex")
