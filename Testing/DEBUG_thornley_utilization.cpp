#include <vector>
#include <string>
#include <Rinternals.h>                            // for debugging
#include <cmath>                                   // for std::isnan
#include "../framework/module.h"                   // for differential_module and update
#include "../framework/module_helper_functions.h"  // for get_ip
#include "../framework/constants.h"                // for getting physical constants
#include "thornley_nutrient_dynamics.h"            // for the thornley_nutrient_dynamics namespace
#include "thornley_utilization.h"

using thornley_nutrient_dynamics::generate_multi_organ_quantity_names;
using thornley_nutrient_dynamics::get_external_substrate_ips;
using thornley_nutrient_dynamics::get_external_substrate_quantity_names;
using thornley_nutrient_dynamics::get_multi_organ_ips;
using thornley_nutrient_dynamics::get_multi_organ_ops;
using thornley_nutrient_dynamics::organ;

thornley_utilization::thornley_utilization(
    std::vector<organ> const& organs,
    state_map const& input_quantities,
    state_map* output_quantities)
    :  // Define basic module properties by passing its name to its parent class
      differential_module(),

      // Store the organ names
      organs(organs),

      // Get pointers to input parameters
      substrate_carbon_source_rate_ips(get_external_substrate_ips(input_quantities, organs)),
      mass_ips(get_multi_organ_ips(input_quantities, organs, "")),
      utilization_rate_ips(get_multi_organ_ips(input_quantities, organs, "utilization_rate")),
      senescence_rate_ips(get_multi_organ_ips(input_quantities, organs, "senescence_rate")),
      respiration_factor_ips(get_multi_organ_ips(input_quantities, organs, "respiration_factor")),
      senescence_reuse_factor_ips(get_multi_organ_ips(input_quantities, organs, "senescence_reuse_factor")),
      // Get reference to input parameters
      stop_growth_dvi(get_input(input_quantities, "stop_growth_dvi")),
      DVI(get_input(input_quantities, "DVI")),  // any difference from {}?

      // Get pointers to output parameters
      structural_carbon_ops(get_multi_organ_ops(output_quantities, organs, "structural_carbon")),
      substrate_carbon_ops(get_multi_organ_ops(output_quantities, organs, "substrate_carbon")),
      respiration_loss_ops(get_multi_organ_ops(output_quantities, organs, "respiration_loss")),
      senescence_loss_ops(get_multi_organ_ops(output_quantities, organs, "senescence_loss"))
{
}

std::vector<std::string> thornley_utilization::get_inputs(std::vector<organ> const& organs)
{
    // List the quantity names that are guaranteed to exist for each organ
    std::vector<std::string> quantities_for_each_organ = {
        "",                         // Mg / ha
        "utilization_rate",         // mol / m^2 / hr
        "senescence_rate",          // mol / m^2 / hr
        "respiration_factor",       // dimensionless
        "senescence_reuse_factor",  // dimensionless
    };

    // Append the organ names as prefixes
    std::vector<std::string> inputs = generate_multi_organ_quantity_names(organs, quantities_for_each_organ);

    // Add the external substrate sources
    std::vector<std::string> external_substrate_source_names = get_external_substrate_quantity_names(organs);  // Mg / ha / hr
    inputs.insert(inputs.end(), external_substrate_source_names.begin(), external_substrate_source_names.end());
    inputs.push_back("stop_growth_dvi");
    inputs.push_back("DVI");
    return inputs;
}

std::vector<std::string> thornley_utilization::get_outputs(std::vector<organ> const& organs)
{
    // List the quantity names that exist for each organ
    std::vector<std::string> quantities_for_each_organ = {
        "structural_carbon",  // mol / m^2
        "substrate_carbon",   // mol / m^2
        "respiration_loss",   // mol / m^2 / hr
        "senescence_loss"     // mol / m^2 / hr
    };

    // Append the organ names as prefixes
    std::vector<std::string> outputs = generate_multi_organ_quantity_names(organs, quantities_for_each_organ);

    return outputs;
}

void thornley_utilization::do_multi_organ_operation() const
{
    // Each substrate pool loses carbon to utilization and may gain carbon from an external
    // source (such as photosynthesis). Each structural carbon supply gains carbon due to
    // utilization, but some of the utilized substrate is lost due to respiration.

    std::string message = "thornley_utilization:\n";
    char buff[128];
    bool found_NaN = false;
    bool print_msg = true;

    // For Pod, utilization is switched off before a threshold DOY
    for (size_t i = 0; i < organs.size(); ++i) {
        double change_in_substrate_pool_per_m2 = -*utilization_rate_ips[i] +  // mol / m^2 / hr
                                                 (*senescence_rate_ips[i]) * (*senescence_reuse_factor_ips[i]);

        sprintf(buff, " organ = %i\n", i);
        message += std::string(buff);

        sprintf(buff, " utilization_rate = %e\n", *utilization_rate_ips[i]);
        message += std::string(buff);

        sprintf(buff, " senescence_rate = %e\n", *senescence_rate_ips[i]);
        message += std::string(buff);

        sprintf(buff, " senescence_reuse_factor = %e\n", *senescence_reuse_factor_ips[i]);
        message += std::string(buff);

        sprintf(buff, " change_in_substrate_pool_per_m2 = %e\n", change_in_substrate_pool_per_m2);
        message += std::string(buff);

        if (substrate_carbon_source_rate_ips[i] && DVI < stop_growth_dvi) {
            // TEMPORARY: We must convert the substrate carbon source rate from
            // Mg / ha / hr to mol / m^2 / hr, assuming that all carbon was
            // converted into biomass in the form of glucose (C6H12O6), i.e.,
            // six assimilated CO2 molecules contribute one glucose molecule.
            // Using the molar mass of glucose in kg / mol, the conversion can
            // be accomplished with the following factor:
            //   (1000 kg / 1 Mg) * (6 C / 1 glucose) * (1 ha / 1e4 m^2)
            //    = 0.6 kg * C * ha / (Mg * glucose * m^2)
            // Essentially, we are partially undoing the conversion performed at
            // the end of `CanAC` and `c3CanAC` that switches the assimilation
            // rate from a molar flux density to a mass flux density. In the
            // long term, it would be nice for BioCro to use molar flux
            // densities everywhere, and this conversion will no longer be
            // required in that case.
            double const cf = 0.6 / physical_constants::molar_mass_of_glucose;             // (mol C / m^2) / (Mg glucose / hr)
            change_in_substrate_pool_per_m2 += *substrate_carbon_source_rate_ips[i] * cf;  // mol C / m^2

            sprintf(buff, " substrate_carbon_source_rate = %e\n", *substrate_carbon_source_rate_ips[i]);
            message += std::string(buff);

            sprintf(buff, " change_in_substrate_pool_per_m2 = %e\n", change_in_substrate_pool_per_m2);
            message += std::string(buff);
        }

        if (std::isnan(change_in_substrate_pool_per_m2)) {
            found_NaN = true;
        }

        // senescence rate
        double litter = *senescence_rate_ips[i] * (1 - *senescence_reuse_factor_ips[i]);                                  // mol C / m^2 / hr
        double respiration_loss = *respiration_factor_ips[i] * *utilization_rate_ips[i];                                  // mol C / m^2 / hr
        double growth_from_utilization_per_m2 = (*utilization_rate_ips[i] - respiration_loss - *senescence_rate_ips[i]);  // mol / m^2 / hr

        update(substrate_carbon_ops[i], change_in_substrate_pool_per_m2);  // mol / m^2
        update(structural_carbon_ops[i], growth_from_utilization_per_m2);  // mol / m^2
        update(respiration_loss_ops[i], respiration_loss);                 // mol C / m^2
        update(senescence_loss_ops[i], litter);                            // mol C / m^2
    }

    if (print_msg && found_NaN) {
        message += "\n";
        Rprintf(message.c_str());
        throw std::range_error("Found NaN in thornley_utilization");
    }
}
