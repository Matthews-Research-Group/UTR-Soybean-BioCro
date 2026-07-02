library(patchwork)

#__________________Plot the Accumulative Allocation_____________
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
      cumulative_gross_canopy   = cumsum(`canopy_gross_assimilation_rate`) / 0.3,
      cumulative_net_canopy   = cumsum(`canopy_assimilation_rate`) / 0.3,
      sim_id              = i
    )
}

component_levels <- c(
  "Pod_growth",   "Pod_respiration",
  "Root_growth",  "Root_respiration",
  "Stem_growth",  "Stem_respiration",
  "Leaf_growth",  "Leaf_respiration"
)
component_labels <- component_levels

pal <- setNames(
  colorRampPalette(c("#882255", "#ff99cc",
                     "#332288", "#a291f2",
                     "#999933", "#ffff8f",
                     "#117733", "#9dfabc"))(8),
  component_labels
)

# Build combined long-format data
cumulative_all <- map_dfr(seq_along(results), ~ calculate_cumulative_values(results[[.x]], .x)) %>%
  group_by(sim_id) %>%
  mutate(time_step = row_number()) %>%
  ungroup()

years <- c(2002, 2004, 2005, 2006)
# Facet label mapper: sim_id -> year
year_labeller <- as_labeller(setNames(as.character(years), seq_along(years)))

# ── Plot 1: Stacked components vs cumulative_gross_canopy line ──────────────────
cumulative_long <- cumulative_all %>%
  select(sim_id, time_step, fractional_doy, cumulative_gross_canopy,
         total_usage_w_leaf_respiration, 
         total_usage_wo_leaf_respiration,
         all_of(component_levels)) %>%
  pivot_longer(cols      = all_of(component_levels),
               names_to  = "component",
               values_to = "cum_value") %>%
  mutate(component = factor(component,
                            levels = component_levels,
                            labels = component_labels))

p_stack <- ggplot(cumulative_long,
                  aes(x = fractional_doy, y = cum_value, fill = component)) +
  geom_area(position = "stack", colour = NA, alpha = 0.85) +
  geom_line(aes(x = fractional_doy, y = cumulative_gross_canopy),
            colour = "black", linewidth = 0.8, linetype = "dashed",
            inherit.aes = FALSE,
            data = cumulative_all) +
  facet_wrap(~ sim_id, labeller = year_labeller, scales = "free_y") +
  scale_fill_manual(values = pal, name = "Carbon Use") +
  labs(title = NULL,
       x = "Day of Year",
       y = expression("Cumulative Carbon Use (mol / m"^2*")")) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

p_stack <- ggplot(cumulative_long,
                  aes(x = fractional_doy, y = cum_value, fill = component)) +
  geom_area(position = "stack", colour = NA, alpha = 0.85) +
  # Map 'linetype' aesthetic to a named string so it appears in legend
  geom_line(aes(x = fractional_doy, y = cumulative_gross_canopy,
                linetype = "Gross Canopy\nAssimilation",
                colour   = "Gross Canopy\nAssimilation"),
            linewidth = 0.8,
            inherit.aes = FALSE,
            data = cumulative_all) +
  facet_wrap(~ sim_id, labeller = year_labeller, scales = "free_y") +
  scale_fill_manual(values = pal, name = "Carbon Use") +
  # Gives the dashed line its colour in the legend
  scale_colour_manual(
    name   = NULL,
    values = c("Gross Canopy\nAssimilation" = "black")
  ) +
  # Gives the dashed line its linetype in the legend
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
    # Ensure the two new legends (colour + linetype) merge into one entry
    legend.key = element_rect(fill = NA)
  ) +
  # Merge colour and linetype guides into a single legend entry
  guides(
    colour   = guide_legend(order = 2, override.aes = list(fill = NA)),
    linetype = guide_legend(order = 2, override.aes = list(fill = NA)),
    fill     = guide_legend(order = 1)
  )

print(p_stack)

ggsave('Fig-stacked_cumulative_carbon_use.png',
       plot   = p_stack,
       width  = 8.5,
       height = 5,
       units  = "in",
       dpi    = 600)

# ── Plot 2: Stacked % of total_usage_w_leaf_respiration (sums to 100%) ───────────────
p_contrib_w_leaf_respiration <- cumulative_long %>%
  mutate(pct = cum_value / total_usage_w_leaf_respiration * 100) %>%
  ggplot(aes(x = fractional_doy, y = pct, fill = component)) +
  # geom_area(position = "stack", colour = NA, alpha = 0.85) + 
  geom_col(position = "stack", width = 0.15) +
  facet_wrap(~ sim_id, labeller = year_labeller) +
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

ggsave('Fig-stacked_cumulative_carbon_use_pct.png',
       plot   = p_contrib_w_leaf_respiration,
       width  = 8.5,
       height = 5,
       units  = "in",
       dpi    = 600)


# Summarise: sum growth + respiration for each organ
organ_long <- cumulative_long %>%
  mutate(
    organ = str_extract(component, "^[A-Za-z]+(?=_)"),  # extract "Leaf", "Stem", etc.
    pct   = cum_value / total_usage_w_leaf_respiration * 100
  ) %>%
  group_by(sim_id, fractional_doy, organ) %>%
  summarise(organ_pct = sum(pct), .groups = "drop")

# Define organ colours (adjust to taste)
organ_pal <- c(
  Leaf = "#117733",
  Stem = "#999933",
  Root = "#332288",
  Pod  = "#882255"
)

p_organ_lines <- organ_long %>%
  ggplot(aes(x = fractional_doy, y = organ_pct,
             colour = organ, group = organ)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ sim_id, labeller = year_labeller) +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = scales::label_number(suffix = "%")
  ) +
  labs(
    x = "Day of Year",
    y = "Proportion of carbon usage of gross canopy assimilation Rate (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.4, "cm"),
    strip.background = element_rect(fill = "grey90")
  )

print(p_organ_lines)

# Summarise organ totals per year first, then compute mean/min/max across years
organ_summary <- organ_long %>%
  group_by(fractional_doy, organ) %>%
  summarise(
    mean_pct = mean(organ_pct),
    min_pct  = min(organ_pct),
    max_pct  = max(organ_pct),
    .groups  = "drop"
  )

p_organ_mean <- organ_summary %>%
  ggplot(aes(x = fractional_doy, colour = organ, fill = organ, group = organ)) +
  # Shaded min–max range
  geom_ribbon(aes(ymin = min_pct, ymax = max_pct), alpha = 0.2, colour = NA) +
  # Mean line on top
  geom_line(aes(y = mean_pct), linewidth = 0.9) +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_fill_manual(values = organ_pal, name = "Organ") +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = scales::label_number(suffix = "%")
  ) +
  labs(
    x = "Day of Year",
    y = "Proportion of carbon usage of gross canopy assimilation Rate (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.4, "cm")
  )

print(p_organ_mean)


p_diff <- cumulative_all %>%
  mutate(diff = total_usage_w_leaf_respiration - cumulative_gross_canopy) %>%
  ggplot(aes(x = fractional_doy, y = diff)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ sim_id, labeller = year_labeller, scales = "free_y") +
  labs(
    title = NULL,
    x     = "Day of Year",
    y     = expression("Total Caron Use - Cumulative Canopy Assimilation Rate (mol C / m"^2*")")
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90")
  )

print(p_diff)
ggsave('FigS-diff_usage_vs_assim_dvi.png',
       plot   = p_diff,
       width  = 8,
       height = 6,
       units  = "in",
       dpi    = 600)

# ── Back-calculate daily growth from cumulative values ────────────────────
growth_long <- cumulative_all %>%
  group_by(sim_id) %>%
  arrange(fractional_doy, .by_group = TRUE) %>%
  mutate(
    Leaf = c(NA, diff(Leaf_growth)),
    Stem = c(NA, diff(Stem_growth)),
    Root = c(NA, diff(Root_growth)),
    Pod  = c(NA, diff(Pod_growth))
  ) %>%
  ungroup() %>%
  select(sim_id, fractional_doy, Leaf, Stem, Root, Pod) %>%
  pivot_longer(cols      = c(Leaf, Stem, Root, Pod),
               names_to  = "organ",
               values_to = "growth_rate") %>%
  mutate(organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")))

use_long <- cumulative_all %>%
  group_by(sim_id) %>%
  arrange(fractional_doy, .by_group = TRUE) %>%
  mutate(
    Leaf = c(NA, diff(Leaf)),
    Stem = c(NA, diff(Stem)),
    Root = c(NA, diff(Root)),
    Pod  = c(NA, diff(Pod))
  ) %>%
  ungroup() %>%
  select(sim_id, fractional_doy, Leaf, Stem, Root, Pod) %>%
  pivot_longer(cols      = c(Leaf, Stem, Root, Pod),
               names_to  = "organ",
               values_to = "growth_rate") %>%
  mutate(organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")))



# ── Plot ──────────────────────────────────────────────────────────────────
# Compute shared y limits across all simulations
y_max <- max(use_long$growth_rate, na.rm = TRUE)
y_min <- min(use_long$growth_rate, na.rm = TRUE)

p_use <- ggplot(use_long,
                   aes(x = fractional_doy, y = growth_rate,
                       colour = organ)) +
  geom_line(linewidth = 0.6, na.rm = TRUE, alpha = 0.6) +
  facet_wrap(~ sim_id, labeller = year_labeller,
             scales = "fixed") +          # <-- shared axes
  scale_colour_manual(
    values = c("Leaf" = "#117733",
               "Stem" = "#999933",
               "Root" = "#332288",
               "Pod"  = "#882255"),
    name = "Organ"
  ) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(title = NULL,
       x     = "Day of Year",
       y     = expression("Organ Carbon Growth Rate (mol / m"^2*" / hour)")) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

print(p_use)
ggsave('FigS-organ_use_rate_hourly.png',
       plot   = p_use,
       width  = 8,
       height = 5,
       units  = "in",
       dpi    = 600)


# ──────────────── plot 2002 graph only ────────────────

# ── Filter to 2002 only (sim_id == 1) ────────────────────────────────────────
cumulative_long_2002 <- cumulative_long %>% filter(sim_id == 1)
cumulative_all_2002  <- cumulative_all  %>% filter(sim_id == 1)


library(gridExtra)

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
    x     = "Day of Year",
    y     = expression("Cumulative Carbon Use (mol / m"^2*")")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5)
  )

p_right <- cumulative_long_2002 %>%
  mutate(pct = cum_value / total_usage_w_leaf_respiration * 100) %>%
  ggplot(aes(x = fractional_doy, y = pct, fill = component)) +
  geom_col(width = 0.15, position = "stack", colour = NA, alpha = 0.85) +
  scale_fill_manual(values = pal, name = "Carbon Use") +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     labels = scales::label_number(suffix = "%")) +
  labs(
    x     = "Day of Year",
    y     = "Proportion of carbon usage (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5)
  )


p_left  <- p_left  + labs(title = "(A)")
p_right <- p_right + labs(title = "(B)") + guides(fill = "none")   # <- the fix

p_combined_2002 <- (p_left + p_right) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))

ggsave("Fig-2002_combined_cumulative_and_pct.png",
       plot   = p_combined_2002,
       width  = 8, height = 4, units = "in", dpi = 600)
