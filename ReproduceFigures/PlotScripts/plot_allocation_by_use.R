library(patchwork)

# Helper function: cumulative carbon allocation from hourly/daily utilization rates
calculate_cumulative_values <- function(df, i) {
  df %>%
    arrange(row_number()) %>%
    mutate(
      Leaf_growth         = cumsum(`Leaf_utilization_rate`),
      Leaf_respiration    = (cumsum(`canopy_gross_assimilation_rate`) - cumsum(`canopy_assimilation_rate`)) / 0.3,
      Stem_growth         = cumsum(`Stem_utilization_rate`) - Stem_respiration_loss,
      Stem_respiration    = Stem_respiration_loss,
      Root_growth         = cumsum(`Root_utilization_rate`) - Root_respiration_loss,
      Root_respiration    = Root_respiration_loss,
      Pod_growth          = cumsum(`Pod_utilization_rate`)  - Pod_respiration_loss,
      Pod_respiration     = Pod_respiration_loss,
      Leaf                = cumsum(`Leaf_utilization_rate`),
      Stem                = cumsum(`Stem_utilization_rate`),
      Root                = cumsum(`Root_utilization_rate`),
      Pod                 = cumsum(`Pod_utilization_rate`),
      total_usage_w_leaf_respiration =
        Leaf_growth + Leaf_respiration +
        Stem_growth + Stem_respiration +
        Root_growth + Root_respiration +
        Pod_growth  + Pod_respiration,
      total_usage_wo_leaf_respiration =
        Leaf_growth +
        Stem_growth + Stem_respiration +
        Root_growth + Root_respiration +
        Pod_growth  + Pod_respiration,
      cumulative_gross_canopy = cumsum(`canopy_gross_assimilation_rate`) / 0.3,
      cumulative_net_canopy   = cumsum(`canopy_assimilation_rate`) / 0.3,
      year_id                 = i
    )
}

years <- c(2002, 2004, 2005, 2006)
year_labeller <- as_labeller(setNames(as.character(years), seq_along(years)))

# ============================================================================
# SECTION 1 — Ambient CO2 only, Growth and respiration separated per organ
# ============================================================================
carbon_use_types <- c(
  "Pod_growth",   "Pod_respiration",
  "Root_growth",  "Root_respiration",
  "Stem_growth",  "Stem_respiration",
  "Leaf_growth",  "Leaf_respiration"
)
carbon_use_labels <- c(
  "Pod growth",   "Pod respiration",
  "Root growth",  "Root respiration",
  "Stem growth",  "Stem respiration",
  "Leaf growth",  "Leaf respiration"
)

pal <- setNames(
  colorRampPalette(c("#882255", "#ff99cc",
                     "#332288", "#a291f2",
                     "#999933", "#ffff8f",
                     "#117733", "#9dfabc"))(8),
  carbon_use_labels
)

cumulative_all <- map_dfr(seq_along(results), 
                          ~ calculate_cumulative_values(results[[.x]], .x)) %>%
  group_by(year_id) %>%
  mutate(time_step = row_number()) %>%
  ungroup()

cumulative_long <- cumulative_all %>%
  select(year_id, time_step, fractional_doy, cumulative_gross_canopy,
         total_usage_w_leaf_respiration, all_of(carbon_use_types)) %>%
  pivot_longer(cols      = all_of(carbon_use_types),
               names_to  = "component",
               values_to = "cum_value") %>%
  mutate(component = factor(component,
                            levels = carbon_use_types,
                            labels = carbon_use_labels))

# ────── Figure S2: Cumulative carbon use vs gross canopy assimilation ─────────
p_stack <- ggplot(cumulative_long,
                  aes(x = fractional_doy, y = cum_value, fill = component)) +
  geom_area(position = "stack", colour = NA, alpha = 0.85) +
  geom_line(aes(x = fractional_doy, y = cumulative_gross_canopy,
                linetype = "Gross Canopy\nAssimilation",
                colour   = "Gross Canopy\nAssimilation"),
            linewidth = 0.8,
            inherit.aes = FALSE,
            data = cumulative_all) +
  facet_wrap(~ year_id, labeller = year_labeller) +
  scale_fill_manual(values = pal, name = "Carbon Use") +
  scale_colour_manual(
    name   = NULL,
    values = c("Gross Canopy\nAssimilation" = "black")
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c("Gross Canopy\nAssimilation" = "dashed")
  ) +
  labs(
    x = "Day of Year",
    y = expression("Cumulative Carbon Use (mol / m"^2*")")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.4, "cm"),
    strip.background = element_rect(fill = "grey90"),
    legend.key = element_rect(fill = NA)
  ) +
  guides(
    colour   = guide_legend(order = 2, override.aes = list(fill = NA)),
    linetype = guide_legend(order = 2, override.aes = list(fill = NA)),
    fill     = guide_legend(order = 1)
  )

print(p_stack)
ggsave(file.path(FIGURE_DIR, 'FigS-stacked_cumulative_carbon_use.png'),
       plot = p_stack, width = 8, height = 5, units = "in", dpi = 600)

# ───── Figure S3: Percentage of accumulated C allocation based on usage ───────
p_contrib_w_leaf_respiration <- cumulative_long %>%
  mutate(pct = cum_value / total_usage_w_leaf_respiration * 100) %>%
  ggplot(aes(x = fractional_doy, y = pct, fill = component)) +
  geom_col(position = "stack", width = 0.15) +
  facet_wrap(~ year_id, labeller = year_labeller) +
  scale_fill_manual(values = pal, name = "Carbon Use") +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     labels = scales::label_number(suffix = "%")) +
  labs(title = NULL,
       x     = "Day of Year",
       y     = "Proportion of carbon usage (%)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

print(p_contrib_w_leaf_respiration)
ggsave(file.path(FIGURE_DIR, 'FigS-stacked_cumulative_carbon_use_pct.png'),
       plot = p_contrib_w_leaf_respiration, 
       width = 8.5, height = 5, units = "in", dpi = 600)

# ────────── Figure 4: 2002 C Allocation (Cumulative + Percentage) ─────────────
cumulative_long_2002 <- cumulative_long %>% filter(year_id == 1)
cumulative_all_2002  <- cumulative_all  %>% filter(year_id == 1)

p_left <- ggplot(cumulative_long_2002,
                 aes(x = fractional_doy, y = cum_value, fill = component)) +
  geom_area(position = "stack", colour = NA, alpha = 0.85) +
  geom_line(aes(x = fractional_doy, y = cumulative_gross_canopy),
            colour      = "black",
            linewidth   = 0.8,
            linetype    = "dashed",
            inherit.aes = FALSE,
            data        = cumulative_all_2002) +
  scale_fill_manual(values = pal, name = "Carbon Use") +
  labs(
    x = "Day of Year",
    y = expression("Cumulative Carbon Use (mol / m"^2*")")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p_right <- cumulative_long_2002 %>%
  mutate(pct = cum_value / total_usage_w_leaf_respiration * 100) %>%
  ggplot(aes(x = fractional_doy, y = pct, fill = component)) +
  geom_col(width = 0.15, position = "stack", colour = NA, alpha = 0.85) +
  scale_fill_manual(values = pal, name = "Carbon Use") +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     labels = scales::label_number(suffix = "%")) +
  labs(
    x = "Day of Year",
    y = "Proportion of carbon usage (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p_left  <- p_left  + labs(title = "(A)") + theme(plot.title = element_text(hjust = 0))
p_right <- p_right + labs(title = "(B)") + theme(plot.title = element_text(hjust = 0)) + guides(fill = "none")

p_combined_2002 <- (p_left + p_right) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))

print(p_combined_2002)
ggsave(file.path(FIGURE_DIR, "Fig-2002_combined_cumulative_and_pct.png"),
       plot = p_combined_2002, width = 8, height = 4, units = "in", dpi = 600)

# ============================================================================
# SECTION 2 — SoyFACE, Ambient + Elevated CO2, Growth and Respiration Combined
# ============================================================================
carbon_use_types  <- c("Pod", "Root", "Stem", "Leaf")
carbon_use_labels <- carbon_use_types

pal <- setNames(
  colorRampPalette(c("#882255", "#332288", "#999933", "#117733"))(4),
  carbon_use_labels
)

organ_pal <- c(
  Leaf = "#117733",
  Stem = "#999933",
  Root = "#332288",
  Pod  = "#882255"
)

cumulative_all <- bind_rows(
  map_dfr(seq_along(results),         
          ~ calculate_cumulative_values(results[[.x]],         .x)) %>% 
    mutate(treatment = "Ambient CO2"),
  map_dfr(seq_along(results.elevCO2), 
          ~ calculate_cumulative_values(results.elevCO2[[.x]], .x)) %>% 
    mutate(treatment = "Elevated CO2")
) %>%
  group_by(treatment, year_id) %>%
  mutate(time_step = row_number()) %>%
  ungroup() %>%
  mutate(treatment = 
           factor(treatment, 
                  levels = c("Ambient CO2", "Elevated CO2")))

treatment_year_labeller <- labeller(
  year_id   = year_labeller,
  treatment = label_value
)

# ──── Figure S5: Difference between cumulative carbon use and gross assimilation ──
p_diff <- cumulative_all %>%
  mutate(diff = total_usage_wo_leaf_respiration - cumulative_net_canopy) %>%
  ggplot(aes(x = fractional_doy, y = diff)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_grid(treatment ~ year_id, 
             labeller = treatment_year_labeller, 
             scales = "free_y") +
  labs(
    title = NULL,
    x     = "Day of Year",
    y     = expression("Total Carbon Use - 
                 Cumulative Canopy Assimilation Rate (mol C / m"^2*")")
  ) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey90"))

print(p_diff)
ggsave(file.path(FIGURE_DIR, 'FigS-diff_usage_vs_assim_doy.png'),
       plot = p_diff, width = 12, height = 6, units = "in", dpi = 600)

# ── Back-calculate daily growth from cumulative values ──────────────────────
usage_long <- cumulative_all %>%
  group_by(treatment, year_id) %>%
  arrange(fractional_doy, .by_group = TRUE) %>%
  mutate(
    Leaf = c(NA, diff(Leaf)),
    Stem = c(NA, diff(Stem)),
    Root = c(NA, diff(Root)),
    Pod  = c(NA, diff(Pod))
  ) %>%
  ungroup() %>%
  select(year_id, treatment, fractional_doy, Leaf, Stem, Root, Pod) %>%
  pivot_longer(cols      = c(Leaf, Stem, Root, Pod),
               names_to  = "organ",
               values_to = "growth_rate") %>%
  mutate(organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")))

y_max <- max(usage_long$growth_rate, na.rm = TRUE)
y_min <- min(usage_long$growth_rate, na.rm = TRUE)

p_growth <- ggplot(usage_long,
                   aes(x = fractional_doy, y = growth_rate, colour = organ)) +
  geom_line(linewidth = 0.6, na.rm = TRUE, alpha = 0.6) +
  facet_grid(treatment ~ year_id, labeller = treatment_year_labeller, scales = "fixed") +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(title = NULL,
       x     = "Day of Year",
       y     = expression("Organ Growth Rate (mol / m"^2*" / hour)")) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

print(p_growth)
ggsave(file.path(FIGURE_DIR, 'FigS-organ_growth_rate_doy.png'),
       plot = p_growth, width = 8, height = 5, units = "in", dpi = 600)

# ============================================================================
# SECTION 3 — LD11 Organ Growth Rate
# ============================================================================
cumulative_ld11 <- map_dfr(seq_along(ld11.results),
                           ~ calculate_cumulative_values(ld11.results[[.x]], .x)) %>%
  group_by(year_id) %>%
  mutate(time_step = row_number()) %>%
  ungroup()

usage_long_ld11 <- cumulative_ld11 %>%
  group_by(year_id) %>%
  arrange(fractional_doy, .by_group = TRUE) %>%
  mutate(
    Leaf = c(NA, diff(Leaf)),
    Stem = c(NA, diff(Stem)),
    Root = c(NA, diff(Root)),
    Pod  = c(NA, diff(Pod))
  ) %>%
  ungroup() %>%
  select(year_id, fractional_doy, Leaf, Stem, Root, Pod) %>%
  pivot_longer(cols      = c(Leaf, Stem, Root, Pod),
               names_to  = "organ",
               values_to = "growth_rate") %>%
  mutate(organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")))

y_max_ld11 <- max(usage_long_ld11$growth_rate, na.rm = TRUE)
y_min_ld11 <- min(usage_long_ld11$growth_rate, na.rm = TRUE)

ld11_years <- c(2021, 2022, 2023, 2024)
ld11_year_labeller <- as_labeller(setNames(as.character(ld11_years), seq_along(ld11_years)))

p_growth_ld11 <- ggplot(usage_long_ld11,
                        aes(x = fractional_doy, y = growth_rate, colour = organ)) +
  geom_line(linewidth = 0.6, na.rm = TRUE, alpha = 0.6) +
  facet_wrap(~ year_id, labeller = ld11_year_labeller) +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_y_continuous(limits = c(y_min_ld11, y_max_ld11)) +
  labs(title = NULL,
       x     = "Day of Year",
       y     = expression("Organ Growth Rate (mol / m"^2*" / hour)")) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

print(p_growth_ld11)
ggsave(file.path(FIGURE_DIR, 'FigS-organ_growth_rate_doy_ld11.png'),
       plot = p_growth_ld11, width = 8, height = 5, units = "in", dpi = 600)


# ============================================================================
# SECTION 4 — Cumulative C allocation plotted against DVI (developmental stage)
# ============================================================================
cumulative_all_dvi <- cumulative_all %>%
  mutate(DVI_rounded = round(DVI, 2)) %>%
  group_by(treatment, year_id, DVI_rounded) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE),
            .groups = "drop")

# ────────────────────────── Difference plot ───────────────────────────────────
p_diff_dvi <- cumulative_all_dvi %>%
  mutate(diff = total_usage_wo_leaf_respiration - cumulative_net_canopy) %>%
  ggplot(aes(x = DVI_rounded, y = diff)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_grid(treatment ~ year_id, labeller = treatment_year_labeller, scales = "free_y") +
  labs(x = "DVI",
       y = expression("Total Carbon Use - Cumulative Canopy Assimilation Rate (mol C / m"^2*")")) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey90"))

print(p_diff_dvi)
ggsave(file.path(FIGURE_DIR, 'FigS-diff_usage_vs_assim_dvi.png'),
       plot = p_diff_dvi, width = 8, height = 6, units = "in", dpi = 600)

# ── Back-calculate daily growth from DVI-rounded cumulative values ─────────
usage_long_dvi <- cumulative_all_dvi %>%
  group_by(treatment, year_id) %>%
  arrange(DVI_rounded, .by_group = TRUE) %>%
  mutate(
    Leaf = c(NA, diff(Leaf)),
    Stem = c(NA, diff(Stem)),
    Root = c(NA, diff(Root)),
    Pod  = c(NA, diff(Pod))
  ) %>%
  ungroup() %>%
  select(year_id, treatment, DVI_rounded, Leaf, Stem, Root, Pod) %>%
  pivot_longer(cols      = c(Leaf, Stem, Root, Pod),
               names_to  = "organ",
               values_to = "growth_rate") %>%
  mutate(organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")))

organ_pct_long_dvi <- usage_long_dvi %>%
  group_by(treatment, year_id, DVI_rounded) %>%
  mutate(
    total_growth = sum(growth_rate, na.rm = TRUE),
    pct          = ifelse(total_growth == 0, NA, growth_rate / total_growth * 100)
  ) %>%
  ungroup()

organ_summary2_dvi <- organ_pct_long_dvi %>%
  filter(!is.na(pct)) %>%
  group_by(treatment, organ, DVI_rounded) %>%
  summarise(
    mean_pct = mean(pct, na.rm = TRUE),
    min_pct  = min(pct,  na.rm = TRUE),
    max_pct  = max(pct,  na.rm = TRUE),
    .groups  = "drop"
  )

ambient_summary2_dvi  <- organ_summary2_dvi %>%
  filter(treatment == "Ambient CO2",
         DVI_rounded <= parameters$stop_growth_dvi)

elevated_summary2_dvi <- organ_summary2_dvi %>%
  filter(treatment == "Elevated CO2",
         DVI_rounded <= parameters$stop_growth_dvi)

# ─── Figure 3: growth-rate based carbon allocation, ambient vs elevated CO2 ───
p_organ_overlap_growth_dvi <- ggplot(mapping = aes(x = DVI_rounded, 
                                                   colour = organ,
                                                   fill = organ, 
                                                   group = organ)) +
  geom_ribbon(data = ambient_summary2_dvi,
              aes(ymin = min_pct, ymax = max_pct),
              alpha = 0.5, colour = NA) +
  geom_line(data = ambient_summary2_dvi,
            aes(y = mean_pct, linetype = "Ambient CO2"),
            linewidth = 0.9) +
  geom_ribbon(data = elevated_summary2_dvi,
              aes(ymin = min_pct, ymax = max_pct),
              alpha = 0.2, colour = NA) +
  geom_line(data = elevated_summary2_dvi,
            aes(y = mean_pct, linetype = "Elevated CO2"),
            linewidth = 0.9) +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_fill_manual(values = organ_pal, name = "Organ") +
  scale_linetype_manual(name   = "Treatment",
                        values = c("Ambient CO2" = "solid", 
                                   "Elevated CO2" = "dashed")) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     labels = scales::label_number(suffix = "%")) +
  coord_cartesian(xlim = c(-1, parameters$stop_growth_dvi)) +
  labs(x = "DVI",
       y = "Proportion of carbon usage (%)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm")) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 0.9),
                                 keywidth = unit(0.9, "cm")),
         colour   = guide_legend(override.aes = list(linewidth = 1.5)),
         fill     = "none")

print(p_organ_overlap_growth_dvi)
ggsave(file.path(FIGURE_DIR, 'Fig-organ_growth_rate_pct_dvi_overlapped.png'),
       plot = p_organ_overlap_growth_dvi, 
       width = 8, height = 4, units = "in", dpi = 600)