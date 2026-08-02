library(patchwork)

# Helper function: cumulative carbon allocation from hourly/daily utilization rates
# calculate_cumulative_values <- function(df, i) {
#   df %>%
#     arrange(row_number()) %>%
#     mutate(
#       year_id = i,
#       # Growth and respiration separated
#       Leaf_growth         = cumsum(`Leaf_utilization_rate`),
#       # Leaf_respiration    = (cumsum(`canopy_gross_assimilation_rate`) - cumsum(`canopy_assimilation_rate`)) * 0.6/180.156E-3,
#       # Leaf_respiration    = (cumsum(`canopy_photorespiration_rate`) + cumsum(`canopy_non_photorespiratory_CO2_release_rate`)) / 0.3, #  0.6/180.156E-3,
#       Leaf_respiration    = cumsum(`canopy_gross_assimilation_rate`) / 0.3 - cumsum(`Leaf_substrate_carbon_source_rate_updated`),
#       Stem_growth         = cumsum(`Stem_utilization_rate`) - Stem_respiration_loss,
#       Stem_respiration    = Stem_respiration_loss,
#       Root_growth         = cumsum(`Root_utilization_rate`) - Root_respiration_loss,
#       Root_respiration    = Root_respiration_loss,
#       Pod_growth          = cumsum(`Pod_utilization_rate`)  - Pod_respiration_loss,
#       Pod_respiration     = Pod_respiration_loss,
#       # Utilization (= growth + respiration)
#       Leaf                = cumsum(`Leaf_utilization_rate`),
#       Stem                = cumsum(`Stem_utilization_rate`),
#       Root                = cumsum(`Root_utilization_rate`),
#       Pod                 = cumsum(`Pod_utilization_rate`),
#       # Total C use including leaf respiration
#       total_usage_w_leaf_respiration =
#         Leaf_growth + Leaf_respiration +
#         Stem_growth + Stem_respiration +
#         Root_growth + Root_respiration +
#         Pod_growth  + Pod_respiration,
#       # Total C use including leaf respiration
#       total_usage_wo_leaf_respiration =
#         Leaf_growth +
#         Stem_growth + Stem_respiration +
#         Root_growth + Root_respiration +
#         Pod_growth  + Pod_respiration,
#       # Cumulative gross canopy assimilation
#       cumulative_gross_canopy = cumsum(`canopy_gross_assimilation_rate`) / 0.3,
#       # Cumulative net canopy assimilation
#       # cumulative_net_canopy   = cumsum(`canopy_assimilation_rate`) / 0.3
#       cumulative_net_canopy = cumsum(`Leaf_substrate_carbon_source_rate_updated`),
#     )
# }


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

cumulative_all <- 
  map(seq_along(results), 
      function(i) calculate_cumulative_values(results[[i]], i)) %>%
  list_rbind() %>%
  group_by(year_id) %>%
  ungroup()

cumulative_long <- cumulative_all %>%
  select(year_id, fractional_doy, cumulative_gross_canopy,
         total_usage_w_leaf_respiration, all_of(carbon_use_types)) %>%
  pivot_longer(cols      = all_of(carbon_use_types),
               names_to  = "component",
               values_to = "cum_value") %>%
  mutate(component = factor(component,
                            levels = carbon_use_types,
                            labels = carbon_use_labels))

# Figure S2: Cumulative carbon use vs gross canopy assimilation
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

# ─ Figure S3: Percentage of accumulated C allocation based on usage ─
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

#  Figure 4: 2002 C Allocation (Cumulative + Percentage) ─
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
  mutate(treatment = "Ambient CO2"),
  map(seq_along(results.elevCO2),
      function(i) calculate_cumulative_values(results.elevCO2[[i]], i)) %>%
    list_rbind() %>%
    mutate(treatment = "Elevated CO2")
) %>%
  group_by(treatment, year_id) %>%
  ungroup() %>%
  mutate(treatment =
           factor(treatment,
                  levels = c("Ambient CO2", "Elevated CO2")))

treatment_year_labeller <- labeller(
  year_id   = year_labeller,
  treatment = label_value
)

#  Figure S5: Difference between cumulative carbon use and net assimilation 
p_diff <- cumulative_all %>%
  mutate(diff = cumulative_net_canopy - total_usage_wo_leaf_respiration) %>%
  ggplot(aes(x = fractional_doy, y = diff)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_grid(treatment ~ year_id, 
             labeller = treatment_year_labeller, 
             scales = "free_y") +
  labs(
    title = NULL,
    x     = "Day of Year",
    y     = expression(" 
            Cumulative Canopy Assimilation - Total Carbon Use (mol C / m"^2*")")
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
    mutate(treatment = "Ambient CO2"),
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
    mutate(treatment = "Elevated CO2")
) %>%
  mutate(treatment = factor(treatment, levels = c("Ambient CO2", "Elevated CO2"))) %>%
  pivot_longer(cols      = c(Leaf, Stem, Root, Pod),
               names_to  = "organ",
               values_to = "utilization_rate") %>%
  mutate(organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")))

y_max <- max(usage_long$utilization_rate, na.rm = TRUE)
y_min <- min(usage_long$utilization_rate, na.rm = TRUE)

p_use <- ggplot(usage_long,
                   aes(x = fractional_doy, y = utilization_rate, colour = organ)) +
  geom_line(linewidth = 0.6, na.rm = TRUE, alpha = 0.6) +
  facet_grid(treatment ~ year_id, labeller = treatment_year_labeller, scales = "fixed") +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(title = NULL,
       x     = "Day of Year",
       y     = expression("Organ Utilization Rate (mol / m"^2*" / hour)")) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

print(p_use)
ggsave(file.path(FIGURE_DIR, 'FigS-organ_utilization_rate_doy.png'),
       plot = p_use, width = 8, height = 5, units = "in", dpi = 600)

# ============================================================================
# SECTION 3 — LD11 Organ Utilization Rate
# ============================================================================
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
       y     = expression("Organ Use Rate (mol / m"^2*" / hour)")) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

print(p_use_ld11)
ggsave(file.path(FIGURE_DIR, 'FigS-organ_utilization_rate_doy_ld11.png'),
       plot = p_use_ld11, width = 8, height = 5, units = "in", dpi = 600)


# ============================================================================
# SECTION 4 — Cumulative C allocation plotted against DVI (developmental stage)
# ============================================================================
# cumulative_all_dvi <- cumulative_all %>%
#   mutate(DVI_rounded = round(DVI, 2)) %>%
#   group_by(treatment, year_id, DVI_rounded) %>%
#   summarise(across(where(is.numeric), 
#                    mean, 
#                    na.rm = TRUE),
#             .groups = "drop")
# 
# #  Difference plot by DVI ─
# p_diff_dvi <- cumulative_all_dvi %>%
#   mutate(diff = cumulative_gross_canopy - total_usage_w_leaf_respiration) %>%
#   ggplot(aes(x = DVI_rounded, y = diff)) +
#   geom_line(linewidth = 0.5) +
#   geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
#   facet_grid(treatment ~ year_id, labeller = treatment_year_labeller, scales = "free_y") +
#   labs(x = "DVI",
#        y = expression("Cumulative Canopy Assimilation - Total Carbon Use (mol C / m"^2*")")) +
#   theme_bw(base_size = 12) +
#   theme(strip.background = element_rect(fill = "grey90"))
# 
# print(p_diff_dvi)
# 
# #  Difference plot by DOY ─
# cumulative_all_doy <- cumulative_all %>%
#   group_by(treatment, year_id, doy) %>%
#   summarise(across(where(is.numeric), 
#                    mean),
#             .groups = "drop")
# p_diff_doy <- cumulative_all_doy %>%
#   mutate(diff = cumulative_gross_canopy - total_usage_w_leaf_respiration) %>%
#   ggplot(aes(x = doy, y = diff)) +
#   geom_line(linewidth = 0.5) +
#   geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
#   facet_grid(treatment ~ year_id, labeller = treatment_year_labeller, scales = "free_y") +
#   labs(x = "DOY",
#        y = expression("Cumulative Canopy Assimilation - Total Carbon Use (mol C / m"^2*")")) +
#   theme_bw(base_size = 12) +
#   theme(strip.background = element_rect(fill = "grey90"))
# 
# print(p_diff_doy)
# 
# ggsave(file.path(FIGURE_DIR, 'FigS-diff_usage_vs_assim_dvi.png'),
#        plot = p_diff_dvi, width = 8, height = 6, units = "in", dpi = 600)
# 
#  Sum organ usage by day of year 
usage_long_doy <- usage_long %>%
  group_by(treatment, year_id, organ, doy) %>%
  summarise(
    daily_use = sum(utilization_rate, na.rm = TRUE), 
    DVI       = mean(DVI, na.rm = TRUE),
    .groups   = "drop",)

#  Convert to % of total use on each day 
organ_pct_long_doy <- usage_long_doy %>%
  group_by(treatment, year_id, doy) %>%
  mutate(
    total_use = sum(daily_use, na.rm = TRUE),
    pct       = ifelse(total_use == 0, NA, daily_use / total_use * 100)
  ) %>%
  ungroup()

#  Summarise across years: mean/min/max % per organ per day ─
organ_use_summary <- organ_pct_long_doy %>%
  filter(!is.na(pct)) %>%
  group_by(treatment, organ, doy) %>%
  summarise(
    mean_pct = mean(pct, na.rm = TRUE),
    min_pct  = min(pct,  na.rm = TRUE),
    max_pct  = max(pct,  na.rm = TRUE),
    DVI      = mean(DVI, na.rm = TRUE),
    .groups  = "drop"
  )

ambient_summary  <- organ_use_summary %>% filter(treatment == "Ambient CO2")
elevated_summary <- organ_use_summary %>% filter(treatment == "Elevated CO2")

# ─ Figure: growth-rate based carbon allocation, ambient vs elevated CO2, by DVI ─
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
  scale_linetype_manual(name   = "Treatment",
                        values = c("Ambient CO2" = "solid", 
                                   "Elevated CO2" = "dashed")) +
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


###############################################################  
# Testing why total use is greater than total assimilation rate
###############################################################
walk(seq_along(results), function(i) {
  df <- results[[i]]
  cat("Year", i, "\n")
  print(df %>% summarise(
    min_leaf = min(Leaf_substrate_carbon, na.rm = TRUE),
    min_stem = min(Stem_substrate_carbon, na.rm = TRUE),
    min_root = min(Root_substrate_carbon, na.rm = TRUE),
    min_pod  = min(Pod_substrate_carbon,  na.rm = TRUE)
  ))
})


#  overdraw at each timestep ─
compute_organ_overdraw <- function(res) {
  res %>%
    mutate(
      leaf_available = lag(Leaf_substrate_carbon, 
                           default = first(Leaf_substrate_carbon)) +
                       Leaf_substrate_carbon_source_rate_updated - 
                       substrate_transport_Leaf_to_Stem,
      stem_available = lag(Stem_substrate_carbon, 
                           default = first(Stem_substrate_carbon)) +
                       substrate_transport_Leaf_to_Stem -
                       substrate_transport_Stem_to_Root - 
                       substrate_transport_Stem_to_Pod,
      root_available = lag(Root_substrate_carbon, 
                           default = first(Root_substrate_carbon)) +
                       substrate_transport_Stem_to_Root,
      pod_available  = lag(Pod_substrate_carbon, 
                           default = first(Pod_substrate_carbon)) +
                       substrate_transport_Stem_to_Pod,
      Leaf = pmax(Leaf_utilization_rate - leaf_available, 0),
      Stem = pmax(Stem_utilization_rate - stem_available, 0),
      Root = pmax(Root_utilization_rate - root_available, 0),
      Pod  = pmax(Pod_utilization_rate  - pod_available,  0)
    ) %>%
    select(fractional_doy, Leaf, Stem, Root, Pod)
}

overdraw_long <-   bind_rows(
  imap(results, \(res, i) compute_organ_overdraw(res) %>% mutate(year_id = i)) %>%
    list_rbind() %>%
    mutate(treatment = "Ambient CO2"),
  imap(results.elevCO2, \(res, i) compute_organ_overdraw(res) %>% mutate(year_id = i)) %>%
    list_rbind() %>%
    mutate(treatment = "Elevated CO2")
) %>%
  mutate(treatment = factor(treatment, levels = c("Ambient CO2", "Elevated CO2"))) %>%
  pivot_longer(cols      = c(Leaf, Stem, Root, Pod),
               names_to  = "organ",
               values_to = "overdraw") %>%
  mutate(organ = factor(organ, levels = c("Leaf", "Stem", "Root", "Pod")))

#  Summary stats derived from the same long data frame (no duplicate logic) ─
overdraw_summary <- overdraw_long %>%
  group_by(treatment, year_id, organ) %>%
  summarise(
    total_overdraw   = sum(overdraw, na.rm = TRUE),
    n_overdraw_steps = sum(overdraw > 0, na.rm = TRUE),
    .groups = "drop"
  )

overdraw_summary_by_year <- overdraw_long %>%
  group_by(treatment, year_id) %>%
  summarise(
    total_overdraw   = sum(overdraw, na.rm = TRUE),
    n_overdraw_steps = sum(overdraw > 0, na.rm = TRUE),
    .groups = "drop"
  )

print(overdraw_summary_by_year)

#  Plot overdraw vs day of year 
p_overdraw <- ggplot(overdraw_long, aes(x = fractional_doy, y = overdraw, colour = organ)) +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  facet_grid(treatment ~ year_id, labeller = treatment_year_labeller, scales = "free_y") +
  scale_colour_manual(values = organ_pal, name = "Organ") +
  labs(
    x = "Day of Year",
    y = expression("Carbon Overdraw (mol / m"^2*" / hour)")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey90"))

print(p_overdraw)
ggsave(file.path(FIGURE_DIR, 'FigS-organ_overdraw_doy.png'),
       plot = p_overdraw, width = 12, height = 6, units = "in", dpi = 600)


# Compare overdraw with negative canopy assimilation rate
for (i in seq_along(results)) {
  r <- results[[i]]
  
  neg_assim_when_leaf_substrate_is_neg <- sum(
    r$Leaf_substrate_carbon_source_rate_updated[which(r$Leaf_substrate_carbon < 0 & r$Leaf_substrate_carbon_source_rate_updated < 0)]
  )
  
  final_diff <- cumulative_all %>%
    filter(treatment == "Ambient CO2", year_id == i) %>%
    mutate(diff = cumulative_net_canopy - total_usage_wo_leaf_respiration) %>%
    arrange(fractional_doy) %>%
    slice_tail(n = 1) %>%
    pull(diff)
}

# compare overdraw with negative canopy assimilation rate — Ambient + Elevated CO2
co2_level_data <- list(
  "Ambient CO2"  = results,
  "Elevated CO2" = results.elevCO2
)

overdraw_vs_diff_summary <- map_dfr(names(co2_level_data), function(co2_level_name) {
  res_list <- co2_level_data[[co2_level_name]]
  
  map_dfr(seq_along(res_list), function(i) {
    r <- res_list[[i]]
    
    neg_assim_when_leaf_substrate_is_neg <- sum(
      r$Leaf_substrate_carbon_source_rate_updated[which(r$Leaf_substrate_carbon < 0 & r$Leaf_substrate_carbon_source_rate_updated < 0)]
    )
    
    final_diff <- cumulative_all %>%
      filter(treatment == co2_level_name, year_id == i) %>%
      mutate(diff = cumulative_net_canopy - total_usage_wo_leaf_respiration) %>%
      arrange(fractional_doy) %>%
      slice_tail(n = 1) %>%
      pull(diff)
    
    tibble(
      co2_level    = co2_level_name,
      year_id      = i,
      neg_assim    = neg_assim_when_leaf_substrate_is_neg,
      final_diff   = final_diff,
      discrepancy  = neg_assim_when_leaf_substrate_is_neg - final_diff
    )
  })
})

print(overdraw_vs_diff_summary)

r <- results[[1]]
r_wdiff <- r %>%
 mutate(
  Leaf_diff = lag(Leaf_substrate_carbon, 
                            default = first(Leaf_substrate_carbon)) - 
    Leaf_utilization_rate +
    canopy_assimilation_rate/0.3 - # Leaf_substrate_carbon_source_rate_updated - 
    substrate_transport_Leaf_to_Stem - 
    Leaf_substrate_carbon - 
    Leaf_substrate_senescence_rate,
  Stem_diff = lag(Stem_substrate_carbon, 
                  default = first(Stem_substrate_carbon)) -
    Stem_utilization_rate +
    substrate_transport_Leaf_to_Stem -
    substrate_transport_Stem_to_Root - 
    substrate_transport_Stem_to_Pod - 
    Stem_substrate_carbon - 
    Stem_substrate_senescence_rate,
  Root_diff = lag(Root_substrate_carbon, 
                       default = first(Root_substrate_carbon)) -
    Root_utilization_rate +
    substrate_transport_Stem_to_Root - 
    Root_substrate_carbon - 
    Root_substrate_senescence_rate,
  Pod_diff  = lag(Pod_substrate_carbon, 
                       default = first(Pod_substrate_carbon))  -
    Pod_utilization_rate + substrate_transport_Stem_to_Pod - 
    Pod_substrate_carbon - 
    Pod_substrate_senescence_rate
  )

# library(lattice)
# xyplot(data=r_wdiff, 
#          Leaf_diff+
#          Stem_diff+
#          Root_diff+
#          Pod_diff~
#          fractional_doy,
#        auto.key=TRUE)

l_diff_sum <- sum(r_wdiff$Leaf_diff)
s_diff_sum <- sum(r_wdiff$Stem_diff)
r_diff_sum <- sum(r_wdiff$Root_diff)
p_diff_sum <- sum(r_wdiff$Pod_diff)
all_diff <- l_diff_sum + s_diff_sum + r_diff_sum + p_diff_sum
print(all_diff)

# plot adjustment to the negative assimilation rate
for(i in 1:4){
  print(xyplot(data=results[[i]], (Leaf_substrate_carbon_source_rate_updated-canopy_assimilation_rate/0.3)~fractional_doy))
  print(xyplot(data=results.elevCO2[[i]], (Leaf_substrate_carbon_source_rate_updated-canopy_assimilation_rate/0.3)~fractional_doy))
}

# ggplot version (combined)
combine_list <- function(lst, co2_level) {
  purrr::imap_dfr(lst, function(df, i) {
    df %>% mutate(yr = years[i], CO2_level = co2_level)
  })
}

plot_data <- bind_rows(
  combine_list(results, "ambient"),
  combine_list(results.elevCO2, "elevated")
) %>%
  mutate(assim_adj = Leaf_substrate_carbon_source_rate_updated - canopy_assimilation_rate * 0.6/180.156E-3)

ggplot(plot_data, aes(x = fractional_doy, y = assim_adj)) +
  geom_point(alpha = 0.4, size = 1) +
  facet_grid(CO2_level ~ yr) +
  labs(
    x = "DOY",
    y = "Canopy Assimilation Adjustment (mol/m^2/hr)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    panel.grid.minor = element_blank()
  )

adjustment_df <- plot_data %>%
  filter(assim_adj > 0) %>%
  select(i, CO2_level, fractional_doy, hour, solar, assim_adj) %>%
  arrange(desc(assim_adj))

print(unique(adjustment_df$hour))
View(adjustment_df)
