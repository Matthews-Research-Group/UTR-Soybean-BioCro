# To run this script, first install BioCro from the `master` branch of
# `biocro-dev`. It may not work with other versions of the BioCro R package.

# Clear workspace
rm(list=ls())

# Load libraries
library(BioCro)
library(lattice)

# Here we want to simulate defoliation by running the simulation until a certain
# time, at which point we suddenly decrease the leaf mass. Let's write an R
# function that does this. This function could be improved in several ways, such
# as allowing the user to specify more than one defoliation event, or more than
yr <- 2003
weather <- read.csv(file = paste0('../Data/Weather_data/', yr,'_Bondville_IL_daylength_wDVI.csv'))
# one tissue types from which to remove mass.
run_biocro_with_defoliation <- function(
  model = soybean,
  weather = weather,
  defoliation_time = 198.5,
  defoliation_fraction = 0.6,
  leaf_name = 'Leaf'
)
{
  # Make sure the weather data includes a `time` column
  weather <- add_time_to_weather_data(weather)

  # If the defoliation time is too large, just run the model as usual.
  # Otherwise, we will run the model in stages.
  if (defoliation_time >= max(weather$time)) {
    run_biocro(
      model$initial_values,
      model$parameters,
      weather,
      model$direct_modules,
      model$differential_modules,
      model$ode_solver
    )
  } else {
    # Truncate the weather data so it stops at the defoliation time
    weather_1 <- weather[weather$time <= defoliation_time, ]

    # Run the model until that time
    result_1 <- run_biocro(
      model$initial_values,
      model$parameters,
      weather_1,
      model$direct_modules,
      model$differential_modules,
      model$ode_solver
    )

    # Get the final values of the differential quantities; these will be the
    # values just before defoliation
    differential_quantities_just_before_defoliation <-
      as.list(result_1[nrow(result_1), names(model$initial_values)])

    # Now reduce the leaf mass
    differential_quantities_just_after_defoliation <-
      differential_quantities_just_before_defoliation

    differential_quantities_just_after_defoliation[[leaf_name]] <-
      differential_quantities_just_before_defoliation[[leaf_name]] * defoliation_fraction


    # Truncate the weather data so it starts at the defoliation time
    weather_2 <- weather[weather$time >= defoliation_time, ]

    # Run the model starting from the defoliation time
    result_2 <- run_biocro(
      differential_quantities_just_after_defoliation,
      model$parameters,
      weather_2,
      model$direct_modules,
      model$differential_modules,
      model$ode_solver
    )

    # Combine both results and return them, making sure to not double count the
    # defoliation time
    rbind(result_1[seq_len(nrow(result_1) - 1), ], result_2)
  }
}

# Run the model with and without defoliation
defoliation_run <- run_biocro_with_defoliation()
regular_run <- run_biocro_with_defoliation(defoliation_time = 1000)

# Plot biomass values
full_result <- rbind(
  within(defoliation_run, {type = 'defoliation'}),
  within(regular_run, {type = 'no defoliation'})
)

print(xyplot(
  Leaf + Stem + Root + Pod ~ time,
  group = type,
  data = full_result,
  type = 'l',
  auto = TRUE,
  grid = TRUE
))
