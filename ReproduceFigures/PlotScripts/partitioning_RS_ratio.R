library(BioCro)
library(lattice)

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

xyplot(
  r_s_ratio ~ DVI,
  group = co2_treatment,
  data = full_res,
  type = 'l',
  auto = TRUE,
  ylab = 'Root / (Leaf + Stem + Pod)',
  ylim = c(0, 2)
)
