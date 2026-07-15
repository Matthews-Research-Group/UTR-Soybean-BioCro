library(BioCro)
library(lattice)
library(ggplot2)

soybean$parameters$time_zone_offset <- NULL

pfunc <- with(soybean, {partial_run_biocro(
  initial_values,
  parameters,
  soybean_weather[['2002']],
  direct_modules,
  differential_modules,
  ode_solver,
  'Catm'
)})

add_ratio <- function(biocro_res) {
  within(biocro_res, {r_s_ratio = Root / (Leaf + Stem + Grain + Shell)})
}

full_res <- rbind(
  within(add_ratio(pfunc(380)), {co2_treatment = 'ambient'}),
  within(add_ratio(pfunc(550)), {co2_treatment = 'elevated'})
)

ggplot(
  full_res,
  aes(x = DVI, y = r_s_ratio, color = co2_treatment, group = co2_treatment)
) +
  geom_line() +
  scale_y_continuous(limits = c(0, NA)) +
  labs(y = 'Root:Shoot Ratio from the Partitioning Model', x = 'DVI', color = 'CO2 treatment') +
  theme_bw()

xyplot(
  r_s_ratio ~ DVI,
  group = co2_treatment,
  data = full_res,
  type = 'l',
  auto = TRUE,
  ylab = 'Root / (Leaf + Stem + Pod)',
  ylim = c(0, NA)
)

print(full_res$DVI[which.max(full_res$r_s_ratio)]) 
print(full_res$kRoot[which.max(full_res$r_s_ratio)]) 

xyplot(
  Root ~ DVI,
  group = co2_treatment,
  data = full_res,
  type = 'l',
  auto = TRUE,
  ylab = 'Root (Mg / ha)',
  ylim = c(0, 2)
)



xyplot(
  DVI ~ fractional_doy,
  group = co2_treatment,
  data = full_res,
  type = 'l',
  auto = TRUE,
  ylab = 'Root (Mg / ha)',
  ylim = c(0, 2)
)