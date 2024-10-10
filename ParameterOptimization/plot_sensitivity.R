rm(list=ls())
library(ppcor)
sensitivity_result <- read.csv('prcc_ranking_pod_1e+06.csv')
load("lhc_output_1e+06.RData")
source('../Data/Soybean-BioCro_Parameters/UTR-parameters.R')
View(sensitivity_result)
View(prcc_result)

# Function to calculate partial rank correlation and p-value for one parameter
calc_partial_rank_corr <- function(X, y, param_index) {
  # Convert to ranks
  X_ranks <- apply(X, 2, rank)
  y_rank <- rank(y)
  
  # Calculate partial rank correlation
  result <- pcor.test(X_ranks[, param_index], y_rank, X_ranks[, -param_index], method = "spearman")
  
  return(c(correlation = result$estimate, p_value = result$p.value))
}

# Apply the function to all parameters
results <- t(sapply(1:ncol(prcc_result$X), function(i) calc_partial_rank_corr(prcc_result$X, prcc_result$y, i)))

# Create a data frame with results
results_df <- data.frame(
  Parameter = arg_names_short,
  Correlation = results[, 1],
  Correlation_abs = abs(results[, 1]),
  P_value = results[, 2]
)
results_df_sorted <- results_df[order(results_df$Correlation_abs, decreasing = TRUE), ]

# Display results
print(results_df_sorted)

# Save sorted results
write.csv(results_df_sorted, "sorted_results_w_p_values.csv")

# Plot the absolute correlation
# Load required libraries
library(ggplot2)
library(dplyr)
library(viridis)

# Process the data
plot_data <- results_df_sorted %>%
  mutate(log_p_value = log10(P_value)) %>%
  mutate(Parameter = factor(Parameter, levels = Parameter[order(Correlation_abs, decreasing = FALSE)])) %>%
  mutate(significant = ifelse(P_value < 0.05, "*", ""))  # Add this line

parameter_labels <- c(
  'carbon_to_mass_factor' = "f",
  'Leaf_utilization_rate_constant' = "r[max[Leaf]]", 
  'Stem_utilization_rate_constant' = "r[max[Stem]]", 
  'Root_utilization_rate_constant' = "r[max[Root]]", 
  'Pod_utilization_rate_constant'  = "r[max[Pod]]", 
  'Leaf_utilization_km' = "K[Leaf]", 
  'Stem_utilization_km' = "K[Stem]", 
  'Root_utilization_km' = "K[Root]",
  'Pod_utilization_km'  = "K[Pod]", # 5,7,8,9
  'respiration_factor'  = "k[res]", # 10
  'substrate_conductance_Leaf_to_Stem' = "sigma[Leaf~to~Stem]", 
  'substrate_conductance_Stem_to_Root' = "sigma[Stem~to~Root]", 
  'substrate_conductance_Stem_to_Pod'  = "sigma[Stem~to~Pod]", #  'transportation_beta_exponent', # 13
  'Leaf_senescence_fraction_max' = "r[sene~max[Leaf]]",
  'Stem_senescence_fraction_max' = "r[sene~max[Stem]]",
  'Root_senescence_fraction_max' = "r[sene~max[Root]]",
  'Leaf_senescence_alpha' = "alpha[Leaf]", 
  'Stem_senescence_alpha' = "alpha[Stem]", 
  'Root_senescence_alpha' = "alpha[Root]",
  'Leaf_senescence_beta' = "beta[Leaf]", 
  'Stem_senescence_beta' = "beta[Stem]", 
  'Root_senescence_beta' = "beta[Root]", # 20, 21, 22,
  'senescence_reuse_factor' = "f[r]",
  'Pod_start_dvi' = "DVI[pod~start]", 
  'stop_growth_dvi' = "DVI[stop~growth]"
)

# Plot
ggplot(plot_data, aes(x = Correlation, y = Parameter, fill = log_p_value)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_viridis(option = "plasma", name = "log10(P-value)") +
  geom_text(aes(label = significant, x = ifelse(Correlation > 0, max(Correlation), min(Correlation))),
            hjust = ifelse(plot_data$Correlation > 0, -0.5, 1.5),
            size = 5, color = "black") +  
  labs(title = "Parameter Correlations",
       x = "Partial rank correlation Coefficient",
       y = "Parameter") +
  scale_y_discrete(labels = function(x) parse(text = parameter_labels[x])) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right")

