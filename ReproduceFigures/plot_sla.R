# Combine all years into one data frame
plot_sla <- function(ExpBiomass, years){
  sla_df <- do.call(rbind, lapply(seq_along(ExpBiomass), function(i) {
    data.frame(
      time = ExpBiomass[[i]]$time,
      SLA  = ExpBiomass[[i]]$SLA,
      year = as.factor(years[i])
    )
  }))
  
  # Plot
  ggplot(sla_df, aes(x = (time/24 + 1), y = SLA, color = year)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    labs(
      x     = "DOY",
      y     = "Specific leaf area (ha / Mg)",
      color = "Year",
      title = "Observed SLA over Time"
    ) +
    theme_bw()
}
