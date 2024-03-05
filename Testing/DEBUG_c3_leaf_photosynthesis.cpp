#include "BioCro.h" // for c3EvapoTrans
#include "c3_leaf_photosynthesis.h"
#include "c3photo.h"    // for c3photoC
#include <Rinternals.h> // for debugging
#include <cmath>        // for std::isnan

using standardBML::c3_leaf_photosynthesis;

string_vector c3_leaf_photosynthesis::get_inputs() {
  return {
      "absorbed_ppfd",               // micromol / (m^2 leaf) / s
      "temp",                        // deg. C
      "rh",                          // dimensionless
      "vmax1",                       // micromole / m^2 / s
      "jmax",                        // micromole / m^2 / s
      "tpu_rate_max",                // micromole / m^2 / s
      "Rd",                          // micromole / m^2 / s
      "b0",                          // mol / m^2 / s
      "b1",                          // dimensionless
      "Gs_min",                      // mol / m^2 / s
      "Catm",                        // micromole / mol
      "atmospheric_pressure",        // Pa
      "O2",                          // mmol / mol
      "theta",                       // dimensionless
      "StomataWS",                   // dimensionless
      "electrons_per_carboxylation", // electron / carboxylation
      "electrons_per_oxygenation",   // electron / oxygenation
      "average_absorbed_shortwave",  // J / (m^2 leaf) / s
      "windspeed",                   // m / s
      "height",                      // m
      "specific_heat_of_air",        // J / kg / K
      "minimum_gbw",                 // mol / m^2 / s
      "windspeed_height",            // m
      "beta_PSII" // dimensionless (fraction of absorbed light that reaches
                  // photosystem II)
  };
}

string_vector c3_leaf_photosynthesis::get_outputs() {
  return {
      "Assim",            // micromol / m^2 /s
      "GrossAssim",       // micromol / m^2 /s
      "Rp",               // micromol / m^2 / s
      "Ci",               // micromole / mol
      "Gs",               // mmol / m^2 / s
      "Cs",               // micromol / m^2 / s
      "RHs",              // dimensionless from Pa / Pa
      "TransR",           // mmol / m^2 / s
      "EPenman",          // mmol / m^2 / s
      "EPriestly",        // mmol / m^2 / s
      "leaf_temperature", // deg. C
      "gbw"               // mol / m^2 / s
  };
}

void c3_leaf_photosynthesis::do_operation() const {
  // Make an initial guess for boundary layer conductance
  double const gbw_guess{1.2}; // mol / m^2 / s

  // Get an initial estimate of stomatal conductance, assuming the leaf is at
  // air temperature
  double const initial_stomatal_conductance =
      c3photoC(absorbed_ppfd, ambient_temperature, ambient_temperature, rh,
               vmax1, jmax, tpu_rate_max, Rd, b0, b1, Gs_min, Catm,
               atmospheric_pressure, O2, theta, StomataWS,
               electrons_per_carboxylation, electrons_per_oxygenation,
               beta_PSII, gbw_guess)
          .Gs; // mmol / m^2 / s

  // Calculate a new value for leaf temperature using the estimate for
  // stomatal conductance
  const ET_Str et =
      c3EvapoTrans(average_absorbed_shortwave, ambient_temperature, rh,
                   windspeed, height, specific_heat_of_air,
                   initial_stomatal_conductance, minimum_gbw, windspeed_height);

  double const leaf_temperature = ambient_temperature + et.Deltat; // deg. C

  // Calculate final values for assimilation, stomatal conductance, and Ci
  // using the new leaf temperature
  const photosynthesis_outputs photo = c3photoC(
      absorbed_ppfd, leaf_temperature, ambient_temperature, rh, vmax1, jmax,
      tpu_rate_max, Rd, b0, b1, Gs_min, Catm, atmospheric_pressure, O2, theta,
      StomataWS, electrons_per_carboxylation, electrons_per_oxygenation,
      beta_PSII, et.boundary_layer_conductance);

  std::string message = "c3_leaf_photosynthesis:\n";
  char buff[128];
  bool print_msg = true;

  if (print_msg && std::isnan(photo.Assim)) {
    sprintf(buff, "  absorbed_ppfd = %e\n", absorbed_ppfd);
    message += std::string(buff);
    sprintf(buff, "  temp = %e\n", ambient_temperature);
    message += std::string(buff);
    sprintf(buff, "  rh = %e\n", rh);
    message += std::string(buff);
    sprintf(buff, "  vmax1 = %e\n", vmax1);
    message += std::string(buff);
    sprintf(buff, "  jmax = %e\n", jmax);
    message += std::string(buff);
    sprintf(buff, "  tpu_rate_max = %e\n", tpu_rate_max);
    message += std::string(buff);
    sprintf(buff, "  Rd = %e\n", Rd);
    message += std::string(buff);
    sprintf(buff, "  b0 = %e\n", b0);
    message += std::string(buff);
    sprintf(buff, "  b1 = %e\n", b1);
    message += std::string(buff);
    sprintf(buff, "  Gs_min = %e\n", Gs_min);
    message += std::string(buff);
    sprintf(buff, "  Catm = %e\n", Catm);
    message += std::string(buff);
    sprintf(buff, "  atmospheric_pressure = %e\n", atmospheric_pressure);
    message += std::string(buff);
    sprintf(buff, "  O2 = %e\n", O2);
    message += std::string(buff);
    sprintf(buff, "  theta = %e\n", theta);
    message += std::string(buff);
    sprintf(buff, "  StomataWS = %e\n", StomataWS);
    message += std::string(buff);
    sprintf(buff, "  electrons_per_carboxylation = %e\n",
            electrons_per_carboxylation);
    message += std::string(buff);
    sprintf(buff, "  electrons_per_oxygenation = %e\n",
            electrons_per_oxygenation);
    message += std::string(buff);
    sprintf(buff, "  average_absorbed_shortwave = %e\n",
            average_absorbed_shortwave);
    message += std::string(buff);
    sprintf(buff, "  windspeed = %e\n", windspeed);
    message += std::string(buff);
    sprintf(buff, "  height = %e\n", height);
    message += std::string(buff);
    sprintf(buff, "  specific_heat_of_air = %e\n", specific_heat_of_air);
    message += std::string(buff);
    sprintf(buff, "  minimum_gbw = %e\n", minimum_gbw);
    message += std::string(buff);
    sprintf(buff, "  windspeed_height = %e\n", windspeed_height);
    message += std::string(buff);
    sprintf(buff, "  beta_PSII = %e\n", beta_PSII);
    message += std::string(buff);

    message += "\n";
    Rprintf(message.c_str());
    throw std::range_error("Found NaN in c3_leaf_photosynthesis");
  }

  // Update the outputs
  update(Assim_op, photo.Assim);
  update(GrossAssim_op, photo.GrossAssim);
  update(Rp_op, photo.Rp);
  update(Ci_op, photo.Ci);
  update(Gs_op, photo.Gs);
  update(Cs_op, photo.Cs);
  update(RHs_op, photo.RHs);
  update(TransR_op, et.TransR);
  update(EPenman_op, et.EPenman);
  update(EPriestly_op, et.EPriestly);
  update(leaf_temperature_op, leaf_temperature);
  update(gbw_op, et.boundary_layer_conductance);
}
