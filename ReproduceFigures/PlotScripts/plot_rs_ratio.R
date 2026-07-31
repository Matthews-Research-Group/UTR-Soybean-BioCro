library(rlang)
# Build combined data frame across all years
root_shoot_all <- do.call(rbind, lapply(1:length(years), function(i) {
  data.frame(
    year           = years[i],
    fractional_doy = results[[i]]$fractional_doy,
    DVI            = results[[i]]$DVI,
    ambient_ratio  = with(results[[i]],      Root / (Leaf + Stem + Pod)),
    elevated_ratio = with(results.elevCO2[[i]], Root / (Leaf + Stem + Pod))
  )
}))

# Pivot to long format for ggplot
root_shoot_long <- root_shoot_all %>%
  pivot_longer(
    cols      = c(ambient_ratio, elevated_ratio),
    names_to  = "treatment",
    values_to = "ratio"
  ) %>%
  mutate(treatment = recode(treatment,
                            ambient_ratio  = "Ambient CO2",
                            elevated_ratio = "Elevated CO2"
  ))

# Plot
root_shoot_ratio_plot <- ggplot(root_shoot_long, 
                                aes(x = fractional_doy, 
                                    y = ratio, 
                                    colour = treatment)) +
  geom_line() +
  facet_wrap(~ year, ncol = 2) +
  scale_colour_manual(values = c("Ambient CO2" = "steelblue", 
                                 "Elevated CO2" = "tomato")) +
  labs(
    x      = "Day of Year",
    y      = "Root:Shoot Ratio",
    colour = "CO2 Level"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

print(root_shoot_ratio_plot)
ggsave(file.path(FIGURE_DIR, 'Fig-root-shoot-ratio.png'),
       plot = root_shoot_ratio_plot,
       width = 6,
       height = 4,
       units = "in",
       dpi = 600
)


root_shoot_mean <- root_shoot_long %>%
  mutate(dvi_bin = round(DVI / 0.01) * 0.01) %>%   # bin on DVI instead of fractional_doy
  group_by(dvi_bin, treatment) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    min_ratio  = min(ratio,  na.rm = TRUE),
    max_ratio  = max(ratio,  na.rm = TRUE),
    .groups = "drop"
  )

root_shoot_mean_plot <- ggplot(root_shoot_mean,
                               aes(x = dvi_bin, colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = min_ratio, ymax = max_ratio), alpha = 0.2, colour = NA) +
  geom_line(aes(y = mean_ratio), linewidth = 0.8) +
  scale_colour_manual(values = c("Ambient CO2"  = "steelblue",
                                 "Elevated CO2" = "tomato"),
                      name = "CO2 Level") +
  scale_fill_manual(values = c("Ambient CO2"  = "steelblue",
                               "Elevated CO2" = "tomato"),
                    name = "CO2 Level") +
  labs(
    x = "DVI",
    y = "Root:Shoot Ratio"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm")
  )

print(root_shoot_mean_plot)

ggsave(file.path(FIGURE_DIR, 'Fig-root-shoot-ratio-mean.png'),
       plot   = root_shoot_mean_plot,
       width  = 4,
       height = 3,
       units  = "in",
       dpi    = 600)


years_partition <- c("2002", "2004", "2005", "2006")

partition_all <- do.call(rbind, lapply(years_partition, function(yr) {
  
  soybean$parameters$time_zone_offset <- NULL
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
  
  add_ratio <- function(biocro_res) {
    within(biocro_res, {
      r_s_ratio      = Root / (Leaf + Stem + Grain + Shell)
      fractional_doy = doy + hour / 24   # <- adjust if your columns differ
    })
  }
  
  res <- rbind(
    within(add_ratio(pfunc(380)), {co2_treatment = "ambient"}),
    within(add_ratio(pfunc(550)), {co2_treatment = "elevated"})
  )
  
  data.frame(
    year           = as.numeric(yr),
    fractional_doy = res$fractional_doy,
    ratio          = res$r_s_ratio,
    treatment      = recode(res$co2_treatment,
                            ambient  = "Ambient CO2",
                            elevated = "Elevated CO2")
  )
}))

# Plot - partitioning model
partition_ratio_plot <- ggplot(partition_all,
                               aes(x = fractional_doy,
                                   y = ratio,
                                   colour = treatment)) +
  geom_line() +
  facet_wrap(~ year, ncol = 2) +
  scale_colour_manual(values = c("Ambient CO2" = "steelblue",
                                 "Elevated CO2" = "tomato")) +
  labs(
    x      = "Day of Year",
    y      = "Root:Shoot Ratio",
    colour = "CO2 Level"
  ) +
  theme_bw() +
  theme(
    legend.position   = "bottom",
    strip.background  = element_rect(fill = "grey90"),
    strip.text        = element_text(face = "bold")
  )

# Combine UTR (left) and partitioning (right), sharing one legend
library(patchwork)

combined_plot <- (root_shoot_ratio_plot + ggtitle("(A) UTR model")) +
  (partition_ratio_plot + ggtitle("(B) Partitioning model")) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

print(combined_plot)

ggsave(file.path(FIGURE_DIR, 'Fig-root-shoot-ratio-combined.png'),
       plot   = combined_plot,
       width  = 10,
       height = 4,
       units  = "in",
       dpi    = 600
)

