## ================================================================
## Plot Root:Shoot ratios from two models side by side. 
## ================================================================
soybean$parameters$time_zone_offset <- NULL

add_ratio <- function(biocro_res) {
  within(biocro_res, {r_s_ratio = Root / (Leaf + Stem + Grain + Shell)})
}
partition_years <- c('2002', '2004', '2005', '2006')

# Build root_shoot_all for the UTR model 
root_shoot_all <- do.call(rbind, lapply(1:length(years), function(i) {
  data.frame(
    year           = years[i],
    fractional_doy = results[[i]]$fractional_doy,
    DVI            = results[[i]]$DVI,
    ambient_ratio  = with(results[[i]],       Root / (Leaf + Stem + Pod)),
    elevated_ratio = with(results.elevCO2[[i]], Root / (Leaf + Stem + Pod))
  )
}))

# Generate partitioning model results for 2002, 2004, 2005, 2006
root_shoot_all_partitioning <- do.call(rbind, lapply(partition_years, function(yr) {
  pfunc <- with(soybean, {
    partial_run_biocro(
      initial_values,
      parameters,
      soybean_weather[[yr]],
      direct_modules,
      differential_modules,
      ode_solver,
      'Catm'
    )
  })
  
  aCO2 <- add_ratio(pfunc(380))
  eCO2 <- add_ratio(pfunc(550))
  
  data.frame(
    year           = yr,
    fractional_doy = aCO2$fractional_doy,
    DVI            = aCO2$DVI,
    ambient_ratio  = aCO2$r_s_ratio,
    elevated_ratio = eCO2$r_s_ratio
  )
}))

# Combine both models' R/S ratio columns, tagging each with "model" 
root_shoot_all_combined <- rbind(
  within(root_shoot_all,              {model = "UTR"}),
  within(root_shoot_all_partitioning, {model = "Partitioning"})
)

# Pivot to long format
root_shoot_long <- root_shoot_all_combined %>%
  pivot_longer(
    cols      = c(ambient_ratio, elevated_ratio),
    names_to  = "treatment",
    values_to = "ratio"
  ) %>%
  mutate(
    treatment = recode(treatment,
                       ambient_ratio  = "ambient CO2",
                       elevated_ratio = "Elevated CO2"
    ),
    model = factor(model, levels = c("UTR", "Partitioning"))  # UTR left, Partitioning right
  )

# Bin on DVI and summarise mean/min/max per model x treatment
root_shoot_mean <- root_shoot_long %>%
  mutate(dvi_bin = round(DVI / 0.01) * 0.01) %>%
  group_by(dvi_bin, treatment, model) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    min_ratio  = min(ratio,  na.rm = TRUE),
    max_ratio  = max(ratio,  na.rm = TRUE),
    .groups = "drop"
  )

# Plot 
root_shoot_mean_plot <- ggplot(root_shoot_mean,
                               aes(x = dvi_bin, colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = min_ratio, ymax = max_ratio), alpha = 0.2, colour = NA) +
  geom_line(aes(y = mean_ratio), linewidth = 0.8) +
  facet_wrap(~ model, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = c("ambient CO2"  = "steelblue",
                                 "Elevated CO2" = "tomato"),
                      name = "CO2 Level") +
  scale_fill_manual(values = c("ambient CO2"  = "steelblue",
                               "Elevated CO2" = "tomato"),
                    name = "CO2 Level") +
  labs(
    x = "DVI",
    y = "Root:Shoot Ratio"
  ) +
  theme_bw() +
  theme(
    legend.position   = "bottom",
    legend.key.size   = unit(0.4, "cm"),
    strip.background  = element_rect(fill = "grey90"),
    strip.text        = element_text(face = "bold")
  )

print(root_shoot_mean_plot)

ggsave(file.path(FIGURE_DIR, 'Fig-root-shoot-ratio-mean-both-models.png'),
       plot   = root_shoot_mean_plot,
       width  = 7,
       height = 3,
       units  = "in",
       dpi    = 600)
