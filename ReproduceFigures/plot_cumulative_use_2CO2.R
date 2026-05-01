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
         total_usage_w_leaf_respiration, all_of(component_levels)) %>%
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
       width  = 7,
       height = 5,
       units  = "in",
       dpi    = 600)

# -------- Plot with net assimilation rate ----------- # 
component_levels <- c(
  "Pod",
  "Root",
  "Stem",
  "Leaf"
)
component_labels <- component_levels

pal <- setNames(
  colorRampPalette(c("#882255",
                     "#332288", 
                     "#999933", 
                     "#117733"))(4),
  component_labels
)


# Build combined long-format data
# ── Build combined data from both treatments ──────────────────────────────────
cumulative_all <- bind_rows(
  map_dfr(seq_along(results),          ~ calculate_cumulative_values(results[[.x]],          .x)) %>% mutate(treatment = "Ambient CO2"),
  map_dfr(seq_along(results.elevCO2),  ~ calculate_cumulative_values(results.elevCO2[[.x]],  .x)) %>% mutate(treatment = "Elevated CO2")
) %>%
  group_by(treatment, sim_id) %>%
  mutate(time_step = row_number()) %>%
  ungroup() %>%
  mutate(treatment = factor(treatment, levels = c("Ambient CO2", "Elevated CO2")))
years <- c(2002, 2004, 2005, 2006)
# Facet label mapper: sim_id -> year
year_labeller <- as_labeller(setNames(as.character(years), seq_along(years)))

treatment_year_labeller <- labeller(
  sim_id    = year_labeller,
  treatment = label_value          # just prints "Ambient CO2" / "Elevated CO2"
)


# ── Plot 1: Stacked components vs cumulative_net_canopy line ──────────────────
cumulative_long <- cumulative_all %>%
  select(sim_id, treatment, time_step, fractional_doy, cumulative_net_canopy,
         total_usage_wo_leaf_respiration, all_of(component_levels)) %>%
  pivot_longer(cols      = all_of(component_levels),
               names_to  = "component",
               values_to = "cum_value") %>%
  mutate(component = factor(component,
                            levels = component_levels,
                            labels = component_labels))

p_stack <- ggplot(cumulative_long,
                  aes(x = fractional_doy, y = cum_value, fill = component)) +
  # geom_area(position = "stack", colour = NA, alpha = 0.85) +
  geom_area(aes(group = interaction(component, treatment)),
            position = "stack", colour = NA, alpha = 0.85) +
  geom_line(aes(x = fractional_doy, y = cumulative_net_canopy),
            colour = "black", linewidth = 0.8, linetype = "dashed",
            inherit.aes = FALSE,
            data = cumulative_all) +
  facet_grid(treatment ~ sim_id, labeller = treatment_year_labeller, scales = "free_y") +
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
  geom_line(aes(x = fractional_doy, y = cumulative_net_canopy,
                linetype = "Net Canopy\nAssimilation",
                colour   = "Net Canopy\nAssimilation"),
            linewidth = 0.8,
            inherit.aes = FALSE,
            data = cumulative_all) +
  facet_grid(treatment ~ sim_id, labeller = treatment_year_labeller, scales = "free_y") +
  scale_fill_manual(values = pal, name = "Carbon Use") +
  # Gives the dashed line its colour in the legend
  scale_colour_manual(
    name   = NULL,
    values = c("Net Canopy\nAssimilation" = "black")
  ) +
  # Gives the dashed line its linetype in the legend
  scale_linetype_manual(
    name   = NULL,
    values = c("Net Canopy\nAssimilation" = "dashed")
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
       width  = 7,
       height = 5,
       units  = "in",
       dpi    = 600)

# ── Plot 2: Stacked % of total_usage_wo_leaf_respiration (sums to 100%) ───────────────
p_contrib <- cumulative_long %>%
  mutate(pct = cum_value / total_usage_wo_leaf_respiration * 100) %>%
  ggplot(aes(x = fractional_doy, y = pct, fill = component)) +
  geom_area(position = "stack", colour = NA, alpha = 0.85) +
  facet_grid(treatment ~ sim_id, labeller = treatment_year_labeller) +
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

print(p_contrib)



# Change to line plots
organ_long <- cumulative_all %>%
  select(sim_id, treatment, fractional_doy, total_usage_wo_leaf_respiration,
         Leaf, Stem, Root, Pod) %>%
  pivot_longer(
    cols      = c(Leaf, Stem, Root, Pod),
    names_to  = "organ",
    values_to = "cum_value"
  ) %>%
  mutate(
    organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")),
    pct   = cum_value / total_usage_wo_leaf_respiration * 100
  ) %>%
  group_by(sim_id, treatment, fractional_doy, organ) %>%
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
  facet_grid(treatment ~ sim_id, labeller = treatment_year_labeller)+
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = scales::label_number(suffix = "%")
  ) +
  labs(
    x = "Day of Year",
    y = "Proportion of cumulative carbon usage (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.4, "cm"),
    strip.background = element_rect(fill = "grey90")
  )

print(p_organ_lines)

organ_summary <- organ_long %>%
  group_by(treatment, fractional_doy, organ) %>%
  summarise(
    mean_pct = mean(organ_pct),
    min_pct  = min(organ_pct),
    max_pct  = max(organ_pct),
    .groups  = "drop"
  )

p_organ_mean <- organ_summary %>%
  ggplot(aes(x = fractional_doy, colour = organ, fill = organ, group = organ)) +
  geom_ribbon(aes(ymin = min_pct, ymax = max_pct), alpha = 0.2, colour = NA) +
  geom_line(aes(y = mean_pct), linewidth = 0.9) +
  facet_wrap(~ treatment) +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_fill_manual(values = organ_pal, name = "Organ") +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = scales::label_number(suffix = "%")
  ) +
  labs(
    x = "Day of Year",
    y = "Proportion of cumulative carbon usage (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.4, "cm"),
    strip.background = element_rect(fill = "grey90")
  )

print(p_organ_mean)



# -----Plot Difference Between Assimilation and Usage ---- #
p_diff <- cumulative_all %>%
  mutate(diff = total_usage_wo_leaf_respiration - cumulative_net_canopy) %>%
  ggplot(aes(x = fractional_doy, y = diff)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_grid(treatment ~ sim_id, labeller = treatment_year_labeller, scales = "free_y") +
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
usage_long <- cumulative_all %>%
  group_by(treatment, sim_id) %>%
  arrange(fractional_doy, .by_group = TRUE) %>%
  mutate(
    Leaf = c(NA, diff(Leaf)),
    Stem = c(NA, diff(Stem)),
    Root = c(NA, diff(Root)),
    Pod  = c(NA, diff(Pod))
  ) %>%
  ungroup() %>%
  select(sim_id, treatment, fractional_doy, Leaf, Stem, Root, Pod) %>%
  pivot_longer(cols      = c(Leaf, Stem, Root, Pod),
               names_to  = "organ",
               values_to = "growth_rate") %>%
  mutate(organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")))

# ── Plot ──────────────────────────────────────────────────────────────────
# Compute shared y limits across all simulations
y_max <- max(usage_long$growth_rate, na.rm = TRUE)
y_min <- min(usage_long$growth_rate, na.rm = TRUE)

p_growth <- ggplot(usage_long,
                   aes(x = fractional_doy, y = growth_rate,
                       colour = organ)) +
  geom_line(linewidth = 0.6, na.rm = TRUE, alpha = 0.6) +
  facet_grid(treatment ~ sim_id, labeller = treatment_year_labeller,
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
       y     = expression("Organ Growth Rate (mol / m"^2*" / hour)")) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

print(p_growth)
ggsave('FigS-organ_growth_rate_dvi.png',
       plot   = p_growth,
       width  = 8,
       height = 5,
       units  = "in",
       dpi    = 600)



ambient_summary  <- organ_summary %>% filter(treatment == "Ambient CO2")
elevated_summary <- organ_summary %>% filter(treatment == "Elevated CO2")

p_organ_overlap <- ggplot(mapping = aes(x = fractional_doy, colour = organ,
                                        fill = organ, group = organ)) +
  # Ambient: filled ribbon (higher alpha) + solid line
  geom_ribbon(data = ambient_summary,
              aes(ymin = min_pct, ymax = max_pct),
              alpha = 0.3, colour = NA) +
  geom_line(data = ambient_summary,
            aes(y = mean_pct, linetype = "Ambient CO2"),
            linewidth = 0.9) +
  # Elevated: filled ribbon (lower alpha) + dashed line
  geom_ribbon(data = elevated_summary,
              aes(ymin = min_pct, ymax = max_pct),
              alpha = 0.1, colour = NA) +
  geom_line(data = elevated_summary,
            aes(y = mean_pct, linetype = "Elevated CO2"),
            linewidth = 0.9) +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_fill_manual(values = organ_pal, name = "Organ") +
  scale_linetype_manual(
    name   = "Treatment",
    values = c("Ambient CO2" = "solid", "Elevated CO2" = "dashed")
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = scales::label_number(suffix = "%")
  ) +
  labs(
    x = "Day of Year",
    y = "Proportion of cumulative carbon usage (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.4, "cm")
  )+
  guides(
    linetype = guide_legend(
      override.aes = list(linewidth = 0.9),
      keywidth = unit(0.9, "cm")   # wider key so solid vs dashed is obvious
    ),
    colour = guide_legend(
      override.aes = list(linewidth = 1.5)  # thicker lines in organ legend
    ),
    fill = "none"  # hide the redundant fill legend
  )

print(p_organ_overlap)

# ── Compute percentages from usage_long ──────────────────────────────────
organ_pct_long <- usage_long %>%
  group_by(treatment, sim_id, fractional_doy) %>%
  mutate(
    total_growth = sum(growth_rate, na.rm = TRUE),
    pct          = ifelse(total_growth == 0, NA, growth_rate / total_growth * 100)
  ) %>%
  ungroup()

# ── Summarise across sim_ids (mean/min/max ribbon) ────────────────────────
organ_summary2 <- organ_pct_long %>%
  group_by(treatment, organ, fractional_doy) %>%
  summarise(
    mean_pct = mean(pct, na.rm = TRUE),
    min_pct  = min(pct,  na.rm = TRUE),
    max_pct  = max(pct,  na.rm = TRUE),
    .groups  = "drop"
  )

ambient_summary2  <- organ_summary2 %>% filter(treatment == "Ambient CO2")
elevated_summary2 <- organ_summary2 %>% filter(treatment == "Elevated CO2")

p_organ_overlap <- ggplot(mapping = aes(x = fractional_doy, colour = organ,
                                        fill = organ, group = organ)) +
  geom_ribbon(data = ambient_summary2,
              aes(ymin = min_pct, ymax = max_pct),
              alpha = 0.3, colour = NA) +
  geom_line(data = ambient_summary2,
            aes(y = mean_pct, linetype = "Ambient CO2"),
            linewidth = 0.9) +
  geom_ribbon(data = elevated_summary2,
              aes(ymin = min_pct, ymax = max_pct),
              alpha = 0.1, colour = NA) +
  geom_line(data = elevated_summary2,
            aes(y = mean_pct, linetype = "Elevated CO2"),
            linewidth = 0.9) +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_fill_manual(values = organ_pal, name = "Organ") +
  scale_linetype_manual(
    name   = "Treatment",
    values = c("Ambient CO2" = "solid", "Elevated CO2" = "dashed")
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = scales::label_number(suffix = "%")
  ) +
  labs(
    x = "Day of Year",
    y = "Proportion of carbon usage (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.4, "cm")
  )+
  guides(
    linetype = guide_legend(
      override.aes = list(linewidth = 0.9),
      keywidth = unit(0.9, "cm")   # wider key so solid vs dashed is obvious
    ),
    colour = guide_legend(
      override.aes = list(linewidth = 1.5)  # thicker lines in organ legend
    ),
    fill = "none"  # hide the redundant fill legend
  )

print(p_organ_overlap)

#########──────────────── Plot against DVI-rounded ────────────────#############
#─── Build DVI-rounded summary data frame ──────────────────────────────────────
cumulative_all_dvi <- cumulative_all %>%
  mutate(DVI_rounded = round(DVI, 2)) %>%
  group_by(treatment, sim_id, DVI_rounded) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), 
            .groups = "drop")

# ── Rebuild long-format data using DVI_rounded ────────────────────────────────
cumulative_long_dvi <- cumulative_all_dvi %>%
  select(sim_id, treatment, DVI_rounded, cumulative_net_canopy,
         total_usage_wo_leaf_respiration, all_of(component_levels)) %>%
  pivot_longer(cols      = all_of(component_levels),
               names_to  = "component",
               values_to = "cum_value") %>%
  mutate(component = factor(component,
                            levels = component_levels,
                            labels = component_labels))

# ── Plot 1 (DVI): Stacked components vs cumulative_net_canopy line ────────────
p_stack_dvi <- ggplot(cumulative_long_dvi,
                      aes(x = DVI_rounded, y = cum_value, fill = component)) +
  geom_area(position = "stack", colour = NA, alpha = 0.85) +
  geom_line(aes(x = DVI_rounded, y = cumulative_net_canopy,
                linetype = "Net Canopy\nAssimilation",
                colour   = "Net Canopy\nAssimilation"),
            linewidth = 0.8,
            inherit.aes = FALSE,
            data = cumulative_all_dvi) +
  facet_grid(treatment ~ sim_id, labeller = treatment_year_labeller, scales = "free_y") +
  scale_fill_manual(values = pal, name = "Carbon Use") +
  scale_colour_manual(name   = NULL,
                      values = c("Net Canopy\nAssimilation" = "black")) +
  scale_linetype_manual(name   = NULL,
                        values = c("Net Canopy\nAssimilation" = "dashed")) +
  labs(x = "DVI",
       y = expression("Cumulative Carbon Use (mol / m"^2*")")) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"),
        legend.key       = element_rect(fill = NA)) +
  guides(colour   = guide_legend(order = 2, override.aes = list(fill = NA)),
         linetype = guide_legend(order = 2, override.aes = list(fill = NA)),
         fill     = guide_legend(order = 1))

print(p_stack_dvi)

ggsave('Fig-stacked_cumulative_carbon_use_CO2_dvi.png',
       plot = p_stack_dvi, width = 8, height = 5, units = "in", dpi = 600)

# ── Plot 2 (DVI): Proportional contribution ───────────────────────────────────
p_contrib_dvi <- cumulative_long_dvi %>%
  mutate(pct = cum_value / total_usage_wo_leaf_respiration * 100) %>%
  ggplot(aes(x = DVI_rounded, y = pct, fill = component)) +
  geom_area(position = "stack", colour = NA, alpha = 0.85) +
  facet_grid(treatment ~ sim_id, labeller = treatment_year_labeller) +
  scale_fill_manual(values = pal, name = "Carbon Use") +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     labels = scales::label_number(suffix = "%")) +
  labs(x = "DVI",
       y = "Proportion of carbon usage (%)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

print(p_contrib_dvi)

# ── Organ long-format (DVI) ───────────────────────────────────────────────────
organ_long_dvi <- cumulative_all_dvi %>%
  select(sim_id, treatment, DVI_rounded, total_usage_wo_leaf_respiration,
         Leaf, Stem, Root, Pod) %>%
  pivot_longer(cols      = c(Leaf, Stem, Root, Pod),
               names_to  = "organ",
               values_to = "cum_value") %>%
  mutate(organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")),
         pct   = cum_value / total_usage_wo_leaf_respiration * 100) %>%
  group_by(sim_id, treatment, DVI_rounded, organ) %>%
  summarise(organ_pct = sum(pct), .groups = "drop")

# ── Per-sim line plot (DVI) ───────────────────────────────────────────────────
p_organ_lines_dvi <- organ_long_dvi %>%
  ggplot(aes(x = DVI_rounded, y = organ_pct, colour = organ, group = organ)) +
  geom_line(linewidth = 0.8) +
  facet_grid(treatment ~ sim_id, labeller = treatment_year_labeller) +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     labels = scales::label_number(suffix = "%")) +
  labs(x = "DVI",
       y = "Proportion of carbon usage (%)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

print(p_organ_lines_dvi)

# ── Mean ± range ribbon plot (DVI) ────────────────────────────────────────────
organ_summary_dvi <- organ_long_dvi %>%
  group_by(treatment, DVI_rounded, organ) %>%
  summarise(mean_pct = mean(organ_pct),
            min_pct  = min(organ_pct),
            max_pct  = max(organ_pct),
            .groups  = "drop")

p_organ_mean_dvi <- organ_summary_dvi %>%
  ggplot(aes(x = DVI_rounded, colour = organ, fill = organ, group = organ)) +
  geom_ribbon(aes(ymin = min_pct, ymax = max_pct), alpha = 0.2, colour = NA) +
  geom_line(aes(y = mean_pct), linewidth = 0.9) +
  facet_wrap(~ treatment) +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_fill_manual(values = organ_pal, name = "Organ") +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     labels = scales::label_number(suffix = "%")) +
  labs(x = "DVI",
       y = "Proportion of carbon usage (%)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

print(p_organ_mean_dvi)

# ── Difference plot (DVI) ─────────────────────────────────────────────────────
p_diff_dvi <- cumulative_all_dvi %>%
  mutate(diff = total_usage_wo_leaf_respiration - cumulative_net_canopy) %>%
  ggplot(aes(x = DVI_rounded, y = diff)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_grid(treatment ~ sim_id, labeller = treatment_year_labeller, scales = "free_y") +
  labs(x = "DVI",
       y = expression("Total Carbon Use - Cumulative Canopy Assimilation Rate (mol C / m"^2*")")) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey90"))

print(p_diff_dvi)

ggsave('FigS-diff_usage_vs_assim_dvi_rounded.png',
       plot = p_diff_dvi, width = 8, height = 6, units = "in", dpi = 600)

# ── Overlapping treatment comparison (DVI) ────────────────────────────────────
ambient_summary_dvi  <- organ_summary_dvi %>% filter(treatment == "Ambient CO2")
elevated_summary_dvi <- organ_summary_dvi %>% filter(treatment == "Elevated CO2")

p_organ_overlap_dvi <- ggplot(mapping = aes(x = DVI_rounded, colour = organ,
                                            fill = organ, group = organ)) +
  geom_ribbon(data = ambient_summary_dvi,
              aes(ymin = min_pct, ymax = max_pct),
              alpha = 0.3, colour = NA) +
  geom_line(data = ambient_summary_dvi,
            aes(y = mean_pct, linetype = "Ambient CO2"),
            linewidth = 0.9) +
  geom_ribbon(data = elevated_summary_dvi,
              aes(ymin = min_pct, ymax = max_pct),
              alpha = 0.1, colour = NA) +
  geom_line(data = elevated_summary_dvi,
            aes(y = mean_pct, linetype = "Elevated CO2"),
            linewidth = 0.9) +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_fill_manual(values = organ_pal, name = "Organ") +
  scale_linetype_manual(name   = "Treatment",
                        values = c("Ambient CO2" = "solid", "Elevated CO2" = "dashed")) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     labels = scales::label_number(suffix = "%")) +
  labs(x = "DVI",
       y = "Proportion of cumulative carbon usage (%)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm")) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 0.9),
                                 keywidth = unit(0.9, "cm")),
         colour   = guide_legend(override.aes = list(linewidth = 1.5)),
         fill     = "none")

print(p_organ_overlap_dvi)


# ── Back-calculate daily growth from DVI-rounded cumulative values ────────
usage_long_dvi <- cumulative_all_dvi %>%
  group_by(treatment, sim_id) %>%
  arrange(DVI_rounded, .by_group = TRUE) %>%
  mutate(
    Leaf = c(NA, diff(Leaf)),
    Stem = c(NA, diff(Stem)),
    Root = c(NA, diff(Root)),
    Pod  = c(NA, diff(Pod))
  ) %>%
  ungroup() %>%
  select(sim_id, treatment, DVI_rounded, Leaf, Stem, Root, Pod) %>%
  pivot_longer(cols      = c(Leaf, Stem, Root, Pod),
               names_to  = "organ",
               values_to = "growth_rate") %>%
  mutate(organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")))

# ── Compute percentages from usage_long_dvi ─────────────────────────────
organ_pct_long_dvi <- usage_long_dvi %>%
  group_by(treatment, sim_id, DVI_rounded) %>%
  mutate(
    total_growth = sum(growth_rate, na.rm = TRUE),
    pct          = ifelse(total_growth == 0, NA, growth_rate / total_growth * 100)
  ) %>%
  ungroup()

# ── Summarise across sim_ids ──────────────────────────────────────────────
organ_summary2_dvi <- organ_pct_long_dvi %>%
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


# ── p_organ_overlap_dvi (growth-rate based) ───────────────────────────────
p_organ_overlap_growth_dvi <- ggplot(mapping = aes(x = DVI_rounded, colour = organ,
                                                   fill = organ, group = organ)) +
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
                        values = c("Ambient CO2" = "solid", "Elevated CO2" = "dashed")) +
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
ggsave('Fig-organ_growth_rate_pct_dvi_overlapped.png',
       plot = p_organ_overlap_growth_dvi, width = 8, height = 4, units = "in", dpi = 600)



