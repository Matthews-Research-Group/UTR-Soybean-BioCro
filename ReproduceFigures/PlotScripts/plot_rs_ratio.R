library(rlang)
library(patchwork)

## ============================================================
## Root:Shoot ratio figures for the UTR and partitioning models.
##
## Writes four figures to FIGURE_DIR:
##   Fig-root-shoot-ratio.png                   UTR, ratio vs DOY by year
##   Fig-root-shoot-ratio-mean.png              UTR, mean ratio vs DVI
##   Fig-root-shoot-ratio-combined.png          UTR + Partitioning by year
##   Fig-root-shoot-ratio-mean-both-models.png  both models, mean ratio vs DVI
##
## Expects the following to already be in scope: `results`,
## `results.elevCO2`, `years`, `soybean`, `soybean_weather`,
## and `FIGURE_DIR`.
## ============================================================

## ---- Shared constants and helpers --------------------------

# Colour for each CO2 treatment
CO2_COLOURS <- c("Ambient CO2" = "steelblue", "Elevated CO2" = "tomato")

# Legend title and key labels, rendering the "2" in CO2 as a subscript.
# CO2_LABELS is positional and must match the order of CO2_COLOURS.
CO2_TITLE  <- expression("CO"[2] * " Level")
CO2_LABELS <- c(expression("Ambient CO"[2]), expression("Elevated CO"[2]))

# CO2 treatments to simulate for the partitioning model, and the years to run.
CO2_PPM <- c(ambient = 372, elevated = 550)
PARTITION_YEARS <- c("2002", "2004", "2005", "2006")

# Theme shared by the "facet by year" panels.
theme_facet <- theme_bw() + theme(
  legend.position  = "bottom",
  strip.background = element_rect(fill = "grey90"),
  strip.text       = element_text(face = "bold")
)

# Pivot a wide (ambient_ratio / elevated_ratio) frame to long format,
# renaming the treatments to their legend labels.
to_long <- function(wide) {
  wide %>%
    pivot_longer(c(ambient_ratio, elevated_ratio),
                 names_to = "co2_level", values_to = "ratio") %>%
    mutate(co2_level = recode(co2_level,
                              ambient_ratio  = "Ambient CO2",
                              elevated_ratio = "Elevated CO2"))
}

# Bin ratios on DVI and summarise mean / min / max per CO2 level
# (optionally split by an extra grouping column such as "model").
summarise_by_dvi <- function(long, extra_group = NULL) {
  long %>%
    mutate(dvi_bin = round(DVI / 0.01) * 0.01) %>%
    group_by(across(all_of(c("dvi_bin", "co2_level", extra_group)))) %>%
    summarise(
      mean_ratio = mean(ratio, na.rm = TRUE),
      min_ratio  = min(ratio,  na.rm = TRUE),
      max_ratio  = max(ratio,  na.rm = TRUE),
      .groups = "drop"
    )
}

# Ratio vs day-of-year, one panel per year.
plot_by_year <- function(long) {
  ggplot(long, aes(fractional_doy, ratio, colour = co2_level)) +
    geom_line() +
    facet_wrap(~ year, ncol = 2) +
    scale_colour_manual(values = CO2_COLOURS, labels = CO2_LABELS) +
    labs(x = "Day of Year", y = "Root:Shoot Ratio", colour = CO2_TITLE) +
    theme_facet
}

# Mean ratio vs DVI with min/max ribbon; optionally facet by model.
plot_mean <- function(mean_df, facet_model = FALSE) {
  p <- ggplot(mean_df, aes(dvi_bin, colour = co2_level, fill = co2_level)) +
    geom_ribbon(aes(ymin = min_ratio, ymax = max_ratio), alpha = 0.2, colour = NA) +
    geom_line(aes(y = mean_ratio), linewidth = 0.8) +
    scale_colour_manual(values = CO2_COLOURS, labels = CO2_LABELS, name = CO2_TITLE) +
    scale_fill_manual(values = CO2_COLOURS, labels = CO2_LABELS, name = CO2_TITLE) +
    labs(x = "DVI", y = "Root:Shoot Ratio") +
    theme_bw() +
    theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"))
  if (facet_model) {
    p <- p + facet_wrap(~ model, ncol = 2, scales = "free_y") +
      theme(strip.background = element_rect(fill = "grey90"),
            strip.text       = element_text(face = "bold"))
  }
  p
}

# Save a figure to FIGURE_DIR at 600 dpi.
save_fig <- function(plot, file, width, height) {
  ggsave(file.path(FIGURE_DIR, file), plot = plot,
         width = width, height = height, units = "in", dpi = 600)
}

## ---- Build datasets (once each) ----------------------------

# UTR model: root:shoot ratio for ambient and elevated CO2, per year.
utr_wide <- do.call(rbind, lapply(seq_along(years), function(i) {
  data.frame(
    year           = years[i],
    fractional_doy = results[[i]]$fractional_doy,
    DVI            = results[[i]]$DVI,
    ambient_ratio  = with(results[[i]],         Root / (Leaf + Stem + Pod)),
    elevated_ratio = with(results.elevCO2[[i]], Root / (Leaf + Stem + Pod))
  )
}))
utr_long <- to_long(utr_wide)

# Partitioning model: run BioCro for each year at both CO2 levels.
soybean$parameters$time_zone_offset <- NULL

# Uses the `fractional_doy` column already present in the BioCro output.
add_ratio <- function(res) within(res, {
  r_s_ratio = Root / (Leaf + Stem + Grain + Shell)
})

partition_wide <- do.call(rbind, lapply(PARTITION_YEARS, function(yr) {
  pfunc <- with(soybean, partial_run_biocro(
    initial_values, parameters, soybean_weather[[yr]],
    direct_modules, differential_modules, ode_solver, 'Catm'
  ))
  aCO2 <- add_ratio(pfunc(CO2_PPM[["ambient"]]))
  eCO2 <- add_ratio(pfunc(CO2_PPM[["elevated"]]))
  data.frame(
    year           = yr,   # character, to match utr_wide$year for bind_rows()
    fractional_doy = aCO2$fractional_doy,
    DVI            = aCO2$DVI,
    ambient_ratio  = aCO2$r_s_ratio,
    elevated_ratio = eCO2$r_s_ratio
  )
}))
partition_long <- to_long(partition_wide)

# Both models tagged and stacked, for the side-by-side mean figure.
both_long <- bind_rows(
  mutate(utr_long,       model = "UTR"),
  mutate(partition_long, model = "Partitioning")
) %>%
  mutate(model = factor(model, levels = c("UTR", "Partitioning")))

## ---- Figures -----------------------------------------------

# UTR: ratio vs DOY by year.
root_shoot_ratio_plot <- plot_by_year(utr_long)
print(root_shoot_ratio_plot)
save_fig(root_shoot_ratio_plot, "Fig-root-shoot-ratio.png", 6, 4)

# UTR: mean ratio vs DVI.
root_shoot_mean_plot <- plot_mean(summarise_by_dvi(utr_long))
print(root_shoot_mean_plot)
save_fig(root_shoot_mean_plot, "Fig-root-shoot-ratio-mean.png", 4, 3)

# UTR (A) and partitioning (B) by year, sharing one legend.
partition_ratio_plot <- plot_by_year(partition_long)
combined_plot <- (root_shoot_ratio_plot + ggtitle("(A) UTR model")) +
  (partition_ratio_plot + ggtitle("(B) Partitioning model")) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")
print(combined_plot)
save_fig(combined_plot, "Fig-root-shoot-ratio-combined.png", 10, 4)

# Both models: mean ratio vs DVI, one panel per model.
root_shoot_mean_both_plot <- plot_mean(
  summarise_by_dvi(both_long, extra_group = "model"),
  facet_model = TRUE
)
print(root_shoot_mean_both_plot)
save_fig(root_shoot_mean_both_plot, "Fig-root-shoot-ratio-mean-both-models.png", 7, 3)
