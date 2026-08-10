library(patchwork)

calculate_cumulative_values <- function(df, i) {
  df %>%
    arrange(row_number()) %>%
    mutate(
      year_id = i,
      # Growth and respiration separated
      Leaf_growth         = Leaf_cumulative_growth,
      Leaf_respiration    = cumulative_gross_assimilation - Leaf_cumulative_net_assimilation,
      Stem_growth         = Stem_cumulative_growth,
      Stem_respiration    = Stem_respiration_loss,
      Root_growth         = Root_cumulative_growth,
      Root_respiration    = Root_respiration_loss,
      Pod_growth          = Pod_cumulative_growth,
      Pod_respiration     = Pod_respiration_loss,
      # Utilization (= growth + respiration)
      Leaf                = Leaf_cumulative_utilization,
      Stem                = Stem_cumulative_utilization,
      Root                = Root_cumulative_utilization,
      Pod                 = Pod_cumulative_utilization,
      # Total C use including leaf respiration
      total_usage_w_leaf_respiration =
        Leaf_growth + Leaf_respiration +
        Stem_growth + Stem_respiration +
        Root_growth + Root_respiration +
        Pod_growth  + Pod_respiration,
      # Total C use including leaf respiration
      total_usage_wo_leaf_respiration =
        Leaf_growth +
        Stem_growth + Stem_respiration +
        Root_growth + Root_respiration +
        Pod_growth  + Pod_respiration,
      # Cumulative gross canopy assimilation
      cumulative_gross_canopy = cumulative_gross_assimilation,
      cumulative_net_canopy = Leaf_cumulative_net_assimilation,
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

cumulative_ambient_all <- 
  map(seq_along(results), 
      function(i) calculate_cumulative_values(results[[i]], i)) %>%
  list_rbind() %>%
  group_by(year_id) %>%
  ungroup()

cumulative_ambient_long <- cumulative_ambient_all %>%
  select(year_id, fractional_doy, cumulative_gross_canopy,
         total_usage_w_leaf_respiration, all_of(carbon_use_types)) %>%
  pivot_longer(cols      = all_of(carbon_use_types),
               names_to  = "component",
               values_to = "cum_value") %>%
  mutate(component = factor(component,
                            levels = carbon_use_types,
                            labels = carbon_use_labels))

# Figure S2: Cumulative carbon use vs gross canopy assimilation
p_ambient_stack <- ggplot(cumulative_ambient_long,
                  aes(x = fractional_doy, y = cum_value, fill = component)) +
  geom_area(position = "stack", colour = NA, alpha = 0.85) +
  geom_line(aes(x = fractional_doy, y = cumulative_gross_canopy,
                linetype = "Gross Canopy\nAssimilation",
                colour   = "Gross Canopy\nAssimilation"),
            linewidth = 0.8,
            inherit.aes = FALSE,
            data = cumulative_ambient_all) +
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
    y = expression(paste("Cumulative Carbon Use (mol", " m"^{-2},")"))
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

print(p_ambient_stack)
ggsave(file.path(FIGURE_DIR, 'FigS-stacked_cumulative_carbon_use.png'),
       plot = p_ambient_stack, width = 8, height = 5, units = "in", dpi = 600)

#  Figure S3: Percentage of accumulated C allocation based on usage 
p_contrib_w_leaf_respiration <- cumulative_ambient_long %>%
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

#  Final (end-of-season) pct of accumulated C allocation per component, by year 
final_pct_by_year <- cumulative_ambient_long %>%
  mutate(pct = cum_value / total_usage_w_leaf_respiration * 100) %>%
  group_by(year_id, component) %>%
  filter(fractional_doy == max(fractional_doy)) %>%
  ungroup() %>%
  mutate(year = years[year_id]) %>%
  select(year, component, pct) %>%
  arrange(year, component)

#  Mean/min/max final pct per component, across years 
final_pct_summary <- final_pct_by_year %>%
  group_by(component) %>%
  summarise(
    mean_pct = mean(pct, na.rm = TRUE),
    min_pct  = min(pct,  na.rm = TRUE),
    max_pct  = max(pct,  na.rm = TRUE),
    .groups  = "drop"
  )

print(final_pct_summary, n = Inf)

#  Figure 4: 2002 C Allocation (Cumulative + Percentage) 
cumulative_ambient_long_2002 <- cumulative_ambient_long %>% filter(year_id == 1)
cumulative_ambient_all_2002  <- cumulative_ambient_all  %>% filter(year_id == 1)

p_left <- ggplot(cumulative_ambient_long_2002,
                 aes(x = fractional_doy, y = cum_value, fill = component)) +
  geom_area(position = "stack", colour = NA, alpha = 0.85) +
  geom_line(aes(x = fractional_doy, y = cumulative_gross_canopy),
            colour      = "black",
            linewidth   = 0.8,
            linetype    = "dashed",
            inherit.aes = FALSE,
            data        = cumulative_ambient_all_2002) +
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
    y = expression(paste("Cumulative Carbon Use (mol", " m"^{-2},")"))
  ) + 
  theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p_right <- cumulative_ambient_long_2002 %>%
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
organ_pal <- c(
  Leaf = "#117733",
  Stem = "#999933",
  Root = "#332288",
  Pod  = "#882255"
)

cumulative_all <- bind_rows(
  map(seq_along(results),
        function(i) calculate_cumulative_values(results[[i]], i)) %>%
  list_rbind() %>%
  mutate(co2_level = "Ambient CO2"),
  map(seq_along(results.elevCO2),
      function(i) calculate_cumulative_values(results.elevCO2[[i]], i)) %>%
    list_rbind() %>%
    mutate(co2_level = "Elevated CO2")
) %>%
  group_by(co2_level, year_id) %>%
  ungroup() %>%
  mutate(co2_level =
           factor(co2_level,
                  levels = c("Ambient CO2", "Elevated CO2")))

co2_level_year_labeller <- labeller(
  year_id   = year_labeller,
  co2_level = label_value
)

#  Figure S5: Difference between cumulative carbon use and net assimilation 
p_diff <- cumulative_all %>%
  mutate(diff = cumulative_net_canopy - total_usage_wo_leaf_respiration) %>%
  ggplot(aes(x = fractional_doy, y = diff)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_grid(co2_level ~ year_id, 
             labeller = co2_level_year_labeller, 
             scales = "free_y") +
  labs(
    title = NULL,
    x     = "Day of Year",
    y     = expression(
      paste("Cumulative Canopy Assimilation - Total Carbon Use", 
            " (mol", " m"^{-2}, ")"))
  ) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey90"))

print(p_diff)
ggsave(file.path(FIGURE_DIR, 'FigS-diff_usage_vs_assim_doy.png'),
       plot = p_diff, width = 12, height = 6, units = "in", dpi = 600)

#  Plot use rate by timestep 
usage_long <- bind_rows(
  imap(results, \(df, i) {
    df %>%
      mutate(year_id = i) %>%
      select(year_id, fractional_doy, doy, DVI,
             Leaf = Leaf_utilization_rate,
             Stem = Stem_utilization_rate,
             Root = Root_utilization_rate,
             Pod  = Pod_utilization_rate)
  }) %>%
    list_rbind() %>%
    mutate(co2_level = "Ambient CO2"),
  imap(results.elevCO2, \(df, i) {
    df %>%
      mutate(year_id = i) %>%
      select(year_id, fractional_doy, doy, DVI,
             Leaf = Leaf_utilization_rate,
             Stem = Stem_utilization_rate,
             Root = Root_utilization_rate,
             Pod  = Pod_utilization_rate)
  }) %>%
    list_rbind() %>%
    mutate(co2_level = "Elevated CO2")
) %>%
  mutate(co2_level = factor(co2_level, levels = c("Ambient CO2", "Elevated CO2"))) %>%
  pivot_longer(cols      = c(Leaf, Stem, Root, Pod),
               names_to  = "organ",
               values_to = "utilization_rate") %>%
  mutate(organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")))

y_max <- max(usage_long$utilization_rate, na.rm = TRUE)
y_min <- min(usage_long$utilization_rate, na.rm = TRUE)

p_use <- ggplot(usage_long,
                   aes(x = fractional_doy, y = utilization_rate, colour = organ)) +
  geom_line(linewidth = 0.6, na.rm = TRUE, alpha = 0.6) +
  facet_grid(co2_level ~ year_id, labeller = co2_level_year_labeller, scales = "fixed") +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(title = NULL,
       x     = "Day of Year",
       y     = expression(paste(
                          "Organ Utilization Rate (mol", 
                          " m"^{-2}, " hour"^{-1}, ")"))) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

print(p_use)
ggsave(file.path(FIGURE_DIR, 'FigS-organ_utilization_rate_doy.png'),
       plot = p_use, width = 8, height = 5, units = "in", dpi = 600)

# # ============================================================================
# # SECTION 3 — LD11 Organ Utilization Rate
# # ============================================================================
cumulative_ld11 <-
  map(seq_along(ld11.results),
      function(i) calculate_cumulative_values(ld11.results[[i]], i)) %>%
  list_rbind() %>%
  group_by(year_id) %>%
  ungroup()

usage_long_ld11 <- imap(ld11.results, \(df, i) {
  df %>%
    mutate(year_id = i) %>%
    select(year_id, fractional_doy,
           Leaf = Leaf_utilization_rate,
           Stem = Stem_utilization_rate,
           Root = Root_utilization_rate,
           Pod  = Pod_utilization_rate)
}) %>%
  list_rbind() %>%
  pivot_longer(cols      = c(Leaf, Stem, Root, Pod),
               names_to  = "organ",
               values_to = "utilization_rate") %>%
  mutate(organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")))

y_max_ld11 <- max(usage_long_ld11$utilization_rate, na.rm = TRUE)
y_min_ld11 <- min(usage_long_ld11$utilization_rate, na.rm = TRUE)

ld11_years <- c(2021, 2022, 2023, 2024)
ld11_year_labeller <- as_labeller(setNames(as.character(ld11_years), seq_along(ld11_years)))

p_use_ld11 <- ggplot(usage_long_ld11,
                        aes(x = fractional_doy, y = utilization_rate, colour = organ)) +
  geom_line(linewidth = 0.6, na.rm = TRUE, alpha = 0.6) +
  facet_wrap(~ year_id, labeller = ld11_year_labeller) +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_y_continuous(limits = c(y_min_ld11, y_max_ld11)) +
  labs(title = NULL,
       x     = "Day of Year",
       y     = expression("Organ Use Rate (mol m"^{-2}*" hour"^{-1}*")")) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

print(p_use_ld11)
ggsave(file.path(FIGURE_DIR, 'FigS-organ_utilization_rate_doy_ld11.png'),
       plot = p_use_ld11, width = 8, height = 5, units = "in", dpi = 600)

#  Sum organ usage by day of year 
usage_long_doy <- usage_long %>%
  group_by(co2_level, year_id, organ, doy) %>%
  summarise(
    daily_use = sum(utilization_rate, na.rm = TRUE), 
    DVI       = mean(DVI, na.rm = TRUE),
    .groups   = "drop",)

#  Convert to % of total use on each day 
organ_pct_long_doy <- usage_long_doy %>%
  group_by(co2_level, year_id, doy) %>%
  mutate(
    total_use = sum(daily_use, na.rm = TRUE),
    pct       = ifelse(total_use == 0, NA, daily_use / total_use * 100)
  ) %>%
  ungroup()

#  Summarise across years: mean/min/max % per organ per day 
organ_use_summary <- organ_pct_long_doy %>%
  filter(!is.na(pct)) %>%
  group_by(co2_level, organ, doy) %>%
  summarise(
    mean_pct = mean(pct, na.rm = TRUE),
    min_pct  = min(pct,  na.rm = TRUE),
    max_pct  = max(pct,  na.rm = TRUE),
    DVI      = mean(DVI, na.rm = TRUE),
    .groups  = "drop"
  )

ambient_summary  <- organ_use_summary %>% filter(co2_level == "Ambient CO2")
elevated_summary <- organ_use_summary %>% filter(co2_level == "Elevated CO2")

#  Figure: growth-rate based carbon allocation, ambient vs elevated CO2, by DVI 
p_organ_overlap_use_dvi <- ggplot(mapping = aes(x = DVI, 
                                                   colour = organ,
                                                   fill = organ, 
                                                   group = organ)) +
  geom_ribbon(data = ambient_summary,
              aes(ymin = min_pct, ymax = max_pct),
              alpha = 0.5, colour = NA) +
  geom_line(data = ambient_summary,
            aes(y = mean_pct, linetype = "Ambient CO2"),
            linewidth = 0.9) +
  geom_ribbon(data = elevated_summary,
              aes(ymin = min_pct, ymax = max_pct),
              alpha = 0.2, colour = NA) +
  geom_line(data = elevated_summary,
            aes(y = mean_pct, linetype = "Elevated CO2"),
            linewidth = 0.9) +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_fill_manual(values = organ_pal, name = "Organ") +
  scale_linetype_manual(name   = expression("CO"[2]*" Level"),
                        values = c("Ambient CO2" = "solid",
                                   "Elevated CO2" = "dashed"),
                        labels = c(expression("Ambient CO"[2]),
                                   expression("Elevated CO"[2]))) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     labels = scales::label_number(suffix = "%")) +
  labs(x = "Development Index (DVI)",
       y = "Proportion of carbon usage (%)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.2, "cm")) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 0.9),
                                 keywidth = unit(0.9, "cm")),
         colour   = guide_legend(override.aes = list(linewidth = 1.5)),
         fill     = "none")

print(p_organ_overlap_use_dvi)
ggsave(file.path(FIGURE_DIR, 'Fig-organ_utilization_rate_pct_doy_overlapped.png'),
       plot = p_organ_overlap_use_dvi,
       width = 7, height = 4, units = "in", dpi = 600)

# ============================================================================
# SECTION 4 — Combined organ carbon-allocation figure: usage (A) vs flux (B), by DVI
# ============================================================================

#  Flux-based daily organ allocation percentage (method from plot_partitioning.R)
calculate_flux_allocation_daily <- function(result) {
  canopy_assim_daily <- aggregate(result$Leaf_substrate_carbon_source_rate_updated,
                                  list(result$doy), FUN = sum)
  avg_dvi_daily <- aggregate(result$DVI, list(result$doy), FUN = mean)

  net_subC_input <- data.frame(Group.1 = avg_dvi_daily$x, x = canopy_assim_daily$x)

  leaf_export_daily <- aggregate(result$substrate_transport_Leaf_to_Stem,
                                 list(result$doy), FUN = sum)
  pod_allocation_daily <- aggregate(result$substrate_transport_Stem_to_Pod,
                                    list(result$doy), FUN = sum)
  root_allocation_daily <- aggregate(result$substrate_transport_Stem_to_Root,
                                     list(result$doy), FUN = sum)
  stem_allocation_daily <- leaf_export_daily - pod_allocation_daily - root_allocation_daily

  leaf_allocation_daily <- canopy_assim_daily - leaf_export_daily

  data.frame(
    doy  = canopy_assim_daily$Group.1,
    DVI  = avg_dvi_daily$x,
    Leaf = 100 * leaf_allocation_daily$x / net_subC_input$x,
    Stem = 100 * stem_allocation_daily$x / net_subC_input$x,
    Root = 100 * root_allocation_daily$x / net_subC_input$x,
    Pod  = 100 * pod_allocation_daily$x / net_subC_input$x
  )
}

flux_pct_long <- bind_rows(
  imap(results, \(result, i) {
    df <- calculate_flux_allocation_daily(result)
    df <- df[1:which.min(abs(df$DVI - parameters$stop_growth_dvi)), ]
    df$year_id <- i
    df
  }) %>%
    list_rbind() %>%
    mutate(co2_level = "Ambient CO2"),
  imap(results.elevCO2, \(result, i) {
    df <- calculate_flux_allocation_daily(result)
    df <- df[1:which.min(abs(df$DVI - parameters$stop_growth_dvi)), ]
    df$year_id <- i
    df
  }) %>%
    list_rbind() %>%
    mutate(co2_level = "Elevated CO2")
) %>%
  mutate(co2_level = factor(co2_level, levels = c("Ambient CO2", "Elevated CO2"))) %>%
  pivot_longer(cols      = c(Leaf, Stem, Root, Pod),
               names_to  = "organ",
               values_to = "pct") %>%
  mutate(organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")))

#  Summarise across years: mean/min/max flux % per organ per day 
flux_summary <- flux_pct_long %>%
  filter(!is.na(pct)) %>%
  group_by(co2_level, organ, doy) %>%
  summarise(
    mean_pct = mean(pct, na.rm = TRUE),
    min_pct  = min(pct,  na.rm = TRUE),
    max_pct  = max(pct,  na.rm = TRUE),
    DVI      = mean(DVI, na.rm = TRUE),
    .groups  = "drop"
  )

flux_ambient_summary  <- flux_summary %>% filter(co2_level == "Ambient CO2")
flux_elevated_summary <- flux_summary %>% filter(co2_level == "Elevated CO2")

#  Panel B: flux-based carbon allocation, ambient vs elevated CO2, by DVI 
p_organ_overlap_flux_dvi <- ggplot(mapping = aes(x = DVI,
                                                   colour = organ,
                                                   fill = organ,
                                                   group = organ)) +
  geom_ribbon(data = flux_ambient_summary,
              aes(ymin = min_pct, ymax = max_pct),
              alpha = 0.5, colour = NA) +
  geom_line(data = flux_ambient_summary,
            aes(y = mean_pct, linetype = "Ambient CO2"),
            linewidth = 0.9) +
  geom_ribbon(data = flux_elevated_summary,
              aes(ymin = min_pct, ymax = max_pct),
              alpha = 0.2, colour = NA) +
  geom_line(data = flux_elevated_summary,
            aes(y = mean_pct, linetype = "Elevated CO2"),
            linewidth = 0.9) +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_fill_manual(values = organ_pal, name = "Organ") +
  scale_linetype_manual(name   = expression("CO"[2]*" Level"),
                        values = c("Ambient CO2" = "solid",
                                   "Elevated CO2" = "dashed"),
                        labels = c(expression("Ambient CO"[2]),
                                   expression("Elevated CO"[2]))) +
  scale_y_continuous(limits = c(-20, 140),
                     breaks = seq(-20, 140, 20),
                     labels = scales::label_number(suffix = "%")) +
  labs(x = "Development Index (DVI)",
       y = "C allocation based on net import (%)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.2, "cm")) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 0.9),
                                 keywidth = unit(0.9, "cm")),
         colour   = guide_legend(override.aes = list(linewidth = 1.5)),
         fill     = "none")

print(p_organ_overlap_flux_dvi)

#  Combine (A) usage-based and (B) flux-based panels 
p_use_labeled  <- p_organ_overlap_use_dvi +
  labs(title = "(A)", y = "C allocation based on usage (%)") +
  theme(plot.title = element_text(hjust = 0))
p_flux_labeled <- p_organ_overlap_flux_dvi +
  labs(title = "(B)") + theme(plot.title = element_text(hjust = 0))

p_organ_use_and_flux_combined <- (p_use_labeled + p_flux_labeled) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))

print(p_organ_use_and_flux_combined)
ggsave(file.path(FIGURE_DIR, 'Fig-organ_use_and_flux_pct_dvi_combined.png'),
       plot = p_organ_use_and_flux_combined,
       width = 10, height = 4.5, units = "in", dpi = 600)
